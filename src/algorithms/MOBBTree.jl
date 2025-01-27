using DataStructures # for queue

@enum PrunedType NONE INFEASIBILITY INTEGRITY DOMINANCE


mutable struct QCRcoefficients
    Q::Vector{Matrix{Float64}}
    c::Vector{Vector{Float64}}
    constant::Vector{Float64}
end

function QCRcoefficients()
    return QCRcoefficients(Vector{Matrix{Float64}}(), Vector{Vector{Float64}}(), Vector{Float64}())
end

# ----------------------------------
# ---- SupportedSolutionPoint ------
# ----------------------------------
mutable struct SupportedSolutionPoint
    x::Set{Vector{Float64}}
    y::Vector{Float64}
    λ :: Vector{Float64}
    is_integer :: Bool 
end

# todo : tolerance
function dominates(a::SupportedSolutionPoint, b::SupportedSolutionPoint)
    if a.y == b.y
        return false
    else 
        return all(a.y .<= b.y)
    end
end

function Base.isapprox(a::SupportedSolutionPoint, b::SupportedSolutionPoint)
    tol = MOI.get(MultiObjectiveBranchBound(), Tolerance())
    return all(abs.(a.y - b.y) .≤ tol)
end

function Base.:show(io::IO, sol::SupportedSolutionPoint)
    println(io, "(y=", sol.y,  
    "; λ=", sol.λ,
    "; is_integer=", sol.is_integer,
    "; x = ", sol.x,
    " )"
    )
end

function Base.:show(io::IO, sols::Vector{SupportedSolutionPoint})
    print(io, "[")
    for sol in sols
        println(io, sol, " , ")
    end
    println(io, "] size=", length(sols))
end

"""
    Return `true` if the given solution point is near to integer under a tolerance.
"""
function _is_integer(x::Vector{Float64})::Bool
    tol = MOI.get(MultiObjectiveBranchBound(), Tolerance())

    for val in x
        if !(abs(val - floor(Int64, val)) < tol || abs(ceil(Int64, val) - val ) < tol)
            return false
        end
    end

    return true
end

# ----------------------------------
# ---------- Node ------------------
# ----------------------------------
mutable struct Node 
    num::Int64                    
    depth::Int64                # depth in tree
    pred::Node                  # predecessor
    succs::Vector{Node}         # successors
    var_idx::Union{Nothing,MOI.VariableIndex}   # index of the chosen variable to be split
    var_bound::Float64            # variable bound
    bound_type :: Int64         # 1 : <= ; 2 : >= 
    bound_ctr :: Union{Nothing,MOI.ConstraintIndex}    # the branching constraint
    activated::Bool             # if the node is active
    pruned::Bool                # if the node is pruned
    pruned_type::PrunedType      # if the node is fathomed, restore pruned type
    deleted::Bool               # if the node is supposed to be deleted
    lower_bound_set::Vector{SupportedSolutionPoint}        # local lower bound set    
    assignment::Dict{MOI.VariableIndex, Float64}  # (varidex, varbound, boundtype)
    qcr_coeff::QCRcoefficients

    Node() = new()

    function Node(num::Int64, depth::Int64 ;
        pred::Node=Node(), succs::Vector{Node}=Vector{Node}(), var_idx=nothing, var_bound::Float64=0.0, bound_type::Int64=0, bound_ctr=nothing
   )
        n = new()
        n.num = num
        n.depth = depth
        n.pred = pred
        n.succs = succs
        n.var_idx = var_idx
        n.var_bound = var_bound
        n.bound_type = bound_type
        n.bound_ctr = bound_ctr

        n.activated = true 
        n.pruned = false
        n.pruned_type = NONE
        n.deleted = false
        n.lower_bound_set = Vector{SupportedSolutionPoint}()
        n.assignment = Dict{MOI.VariableIndex, Float64}()
        n.qcr_coeff = QCRcoefficients()

        f(t) = nothing 
        finalizer(f, n)
    end
end

function Base.:show(io::IO, n::Node)
    println(io, "\n\n # ----------- node $(n.num) : \n", 
    "depth = $(n.depth) \n",
    "pred = $(n.pred.num) \n",
    "var[ $(n.var_idx) ] = $(n.var_bound) \n",
    "activated = $(n.activated) \n",
    "pruned = $(n.pruned) \n",
    "pruned_type = $(n.pruned_type)"
    )
    print(io, "succs = [ ")
    for s in n.succs print(io, "$(s.num), ") end
    println(io, " ]")

    println(io, "LBS = ", n.lower_bound_set)
end


"""
Return `true` if the given node is the root of a branch-and-bound tree.
"""
function isRoot(node::Node)
    return node.depth == 0 # !isdefined(node, :pred) 
end

"""
Return `true` if the given `node` has activated/non-explored successor(s).
"""
function hasNonExploredChild(node::Node)
    for c in node.succs
        if c.activated return true end 
    end
    return false
end

"""
Delete the given node in B&B tree. (should be a private function)
"""
function Base.delete!(node::Node)           # todo : check
    node.deleted = true ; node = nothing               # remove from the memory
end

"""
Prune the given node in a B&B tree and delete all successors of the pruned node.
"""
function prune!(node::Node, reason::PrunedType)     # todo : check 
    node.pruned = true
    node.pruned_type = reason

    to_delete = node.succs[:]
    node.succs = Vector{Node}()

    while length(to_delete) > 0
        n = pop!(to_delete)
        to_delete = vcat(to_delete, n.succs[:])
        delete!(n)
    end
end


"""
From the actual node, go up to the root to get the partial assignment of variables.
"""
function getPartialAssign(actual::Node)::Dict{MOI.VariableIndex, Float64} 
    assignment = Dict{MOI.VariableIndex, Float64}()
    if isRoot(actual) # the actual node is the root 
        return assignment
    end
    predecessor = actual.pred
    assignment[actual.var_idx] = actual.var_bound

    while !isRoot(predecessor)     
        actual = predecessor ; predecessor = actual.pred
        if actual.bound_ctr !== nothing
            assignment[actual.var_idx] = actual.var_bound
        end
    end
    return assignment
end


"""
Going through all the predecessors until the root, add variables or objective bounds branched in the predecessors.

Return a list of objective bounds (symbolic constraint).
"""
function setVarBounds(actual::Node, model, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}})
    if isRoot(actual) # the actual node is the root 
        return 
    end
    predecessor = actual.pred

    # set actual objective/variable bound
    addOneBound(actual, model, Bounds)

    # set actual objective/variable bounds in predecessors
    while !isRoot(predecessor)    
        actual = predecessor ; predecessor = actual.pred 
        addOneBound(actual, model, Bounds)
    end
end


"""
Remove variables or objective bounds set in the predecessors.
"""
function removeVarBounds(actual::Node, model, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}})
    if isRoot(actual) # the actual node is the root 
        return 
    end
    predecessor = actual.pred
    removeOneBound(actual, model, Bounds)

    while !isRoot(predecessor)     
        actual = predecessor ; predecessor = actual.pred
        removeOneBound(actual, model, Bounds)
    end
end

function removeOneBound(actual::Node, model, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}})
    MOI.delete(model, actual.bound_ctr)

    if actual.bound_type == 1 
        Bounds[2][actual.var_idx] = MOI.add_constraint(model, actual.var_idx, MOI.LessThan(1.0))
    elseif actual.bound_type == 2
        Bounds[1][actual.var_idx] = MOI.add_constraint(model, actual.var_idx, MOI.GreaterThan(0.0))
    else
        error("bound_type unknown in removeVarBounds() .")
    end
end


"""
Given a partial assignment on variables values, add the corresponding bounds.
"""
function addOneBound(actual::Node, model, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}})
    lower_bounds, upper_bounds = Bounds # >= and <= 
    if actual.bound_type == 1 
        MOI.delete(model, upper_bounds[actual.var_idx])
        actual.bound_ctr = MOI.add_constraint(model, actual.var_idx, MOI.LessThan(actual.var_bound))
    elseif actual.bound_type == 2
        MOI.delete(model, lower_bounds[actual.var_idx])
        actual.bound_ctr = MOI.add_constraint(model, actual.var_idx, MOI.GreaterThan(actual.var_bound))
    else
        error("bound_type unknown in addOneBound() .")
    end
end



# ----------------------------------
# ---------- branching -------------
# ----------------------------------
"""
Return an initialized todo list according to the fixed parameter.
"""
function initTree(algorithm)
    if MOI.get(algorithm, TraverseOrder()) == :bfs
        return Queue{Base.RefValue{Node}}()
    elseif MOI.get(algorithm, TraverseOrder()) == :dfs
        return Stack{Base.RefValue{Node}}() 
    elseif MOI.get(algorithm, TraverseOrder()) == :arbitrary
        return Vector{Base.RefValue{Node}}()
    else
        @error "Unknown traverse parameter $(MOI.get(algorithm, TraverseOrder()))\n please set attribute with :bfs, :dfs or :arbitrary ."
    end
end


"""
Add a node identify in the todo list.
"""
function addTree(todo, algorithm, node::Node)
    if MOI.get(algorithm, TraverseOrder()) == :bfs
        enqueue!(todo, Ref(node))
    elseif MOI.get(algorithm, TraverseOrder()) == :dfs
        push!(todo, Ref(node)) 
    elseif MOI.get(algorithm, TraverseOrder()) == :arbitrary
        push!(todo, Ref(node))
    else
        @error "Unknown traverse parameter $(MOI.get(algorithm, TraverseOrder()))\n please set attribute with :bfs, :dfs or :arbitrary ."
    end
end

"""
Return the next element in the todo list.
"""
function nextNodeTree(todo, algorithm)
    if MOI.get(algorithm, TraverseOrder()) == :bfs
        return dequeue!(todo)
    elseif MOI.get(algorithm, TraverseOrder()) == :dfs
        return pop!(todo) 
    elseif MOI.get(algorithm, TraverseOrder()) == :arbitrary
        i = rand(1:length(todo))
        next = todo[i]
        deleteat!(todo, i)
        return next
    else
        @error "Unknown traverse parameter $(MOI.get(algorithm, TraverseOrder()))\n please set attribute with :bfs, :dfs or :arbitrary ."
    end
end

"""
Pick up a free variable to be split according to the prefiexd strategy.
"""
# todo : add other strategies ...
function pickUpAFreeVar(assignment::Dict{MOI.VariableIndex, Float64}, model) :: Union{Nothing, MOI.VariableIndex}
    vars_idx = MOI.get(model, MOI.ListOfVariableIndices())

    if length(assignment) == length(vars_idx) return nothing end # todo : for binary vars only !!

    for v in vars_idx
        if !haskey(assignment, v) return v end 
    end
    return nothing
end


# ----------------------------------
# ---------- Bound Sets ------------
# ----------------------------------
function _fix_λ(λ_count, p, Λ)
    λ = Λ[end]
    # in case of LowerBoundsLimit() > p+1
    for i in 1:p
        if λ_count <= 0 return end
        λ_ = (λ .+ Λ[i]) ./2
        push!(Λ, λ_) ; λ_count -= 1
    end

    while λ_count > 0
        i = rand(1:length(Λ)) ; j = rand(1:length(Λ))
        if i != j
            λ_ = (Λ[j] + Λ[i]) ./2
            push!(Λ, λ_) ; λ_count -= 1
        end
    end
end

# todo : add equivalent solutions
function push_avoiding_duplicate(vec::Vector{SupportedSolutionPoint}, candidate::SupportedSolutionPoint) :: Bool
    for sol in vec
        if sol ≈ candidate return false end 
    end
    push!(vec, candidate) ; return true
end

"""
Stop looking for lower bounds if duplicate is encounterd
"""
function MOLP(algorithm, 
                model::Optimizer, 
                node::Node;
    )

    Λ = [] ; p = MOI.output_dimension(model.f)

    # convexifier each p objective function 
    for i in 1:p
        λ = zeros(p) ; λ[i] = 1.0
        push!(Λ, λ) 

        if MOI.get(algorithm, ConvexQCR())
            is_solved = QCR_csdp(algorithm.Qs[i], zeros(algorithm.nb_vars), 0.0, 
                                    model, algorithm, node.qcr_coeff
                        )
            is_solved ? nothing : return
        end
    end

    λ = [1/p for _ in 1:p] ; push!(Λ, λ)
    λ_count = MOI.get(algorithm, LowerBoundsLimit()) - length(Λ)
    _fix_λ(λ_count, p, Λ)

    iter = 0
    for λ in Λ
        iter += 1
        # status, x, y = solve_weighted_sum(model, λ, MOI.get(algorithm, ConvexQCR()), node.qcr_coeff)

        # if status==OPTIMAL println("λ = ", λ , "\t qcr value = ", λ'*y ) end 
        println("λ ", λ)

        status, x, y = column_generation_algorithm(model, algorithm, sum( λ[i].* node.qcr_coeff.Q[i] for i in 1:p), 
                                        sum( λ[i].* node.qcr_coeff.c[i] for i in 1:p ), λ'* node.qcr_coeff.constant, 
                                        sum( λ[i].* algorithm.Qs[i] for i in 1:p)
                                    )

        # if status==OPTIMAL println("DWR value = ", λ'*y) end 

        if _is_scalar_status_optimal(status)
           
            sol = SupportedSolutionPoint(Set( [ x ] ), y, λ, _is_integer( x ) ) 
            
            if sol.is_integer
                x_val = round.(Int64, first(sol.x)) .*1.0
                sol.x = Set( [x_val])
                sol.y = [ x_val'* algorithm.Qs[i] *x_val for i in 1:p]
            end
            # if any(test -> test ≈ sol, node.lower_bound_set)
            #     nothing
            # else
            #     is_new_point = push_avoiding_duplicate(node.lower_bound_set, sol)
            #     if !is_new_point return end
            # end

            push_filtering_dominance(node.lower_bound_set, sol)
        end
        if iter==2 && length(node.lower_bound_set) < iter return end  
    end

end

"""
Compute and stock the relaxed bound set (i.e. the LP relaxation) of the (sub)problem defined by the given node.
Return `true` if the node is pruned by infeasibility.
"""
function computeLBS(node::Node, model::Optimizer, algorithm, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}})::Bool
    setVarBounds(node, model, Bounds)
    
    MOLP(algorithm, model, node)

    removeVarBounds(node, model, Bounds)

    return length(node.lower_bound_set) == 0
end


function push_filtering_dominance(vec::Vector{SupportedSolutionPoint}, candidate::SupportedSolutionPoint)
    i = 0 ; to_delete = []

    for sol in vec
        i += 1

        if sol ≈ candidate
            # Point already added to nondominated solutions. Don't add
            for equiv in candidate.x push!(sol.x, equiv) end 
            return
        elseif dominates(sol, candidate)
            # Point is dominated. Don't add
            return
        elseif dominates(candidate, sol)
            # new dominating point 
            push!(to_delete, i)
        end
    end

    deleteat!(vec, to_delete) ; push!(vec, candidate)
    sort!(vec; by = sol ->  sol.y) ; 
end

"""
At the given node, update (filtered by dominance) the global upper bound set.
Return `true` if the node is pruned by integrity.
"""
function updateUBS(node::Node, UBS::Vector{SupportedSolutionPoint})::Bool
    integrity = false
    if length(node.lower_bound_set) ==1 && node.lower_bound_set[1].is_integer
        s = node.lower_bound_set[1] ; push_filtering_dominance(UBS, s)
        integrity = true
    else
        for i in 1:length(node.lower_bound_set) 
            if node.lower_bound_set[i].is_integer
                s = node.lower_bound_set[i] ; push_filtering_dominance(UBS, s)
            end
        end
    end
    println("UBS = ", UBS)

    return integrity
end

# ----------------------------------
# ---------- fathoming -------------
# ----------------------------------
"""
Return local nadir points (so-called corner points) of the given UBS.
"""
# todo : p>3 ?!!!!!!!!!!
function getNadirPoints(UBS::Vector{SupportedSolutionPoint}, model) :: Vector{SupportedSolutionPoint}
    p = MOI.output_dimension(model.f)
    nadir_pts = Vector{SupportedSolutionPoint}()

    if length(UBS) == 1 return UBS end 

    if p == 2
        for i in 1:length(UBS)-1
            push!(nadir_pts, SupportedSolutionPoint(Set{Vector{Float64}}(), 
                                                    [UBS[i+1].y[1], UBS[i].y[2]], 
                                                    Vector{Float64}(), false
                                                )
            )
        end
    else
        nothing     # todo p > 3
    end
    return nadir_pts
end

"""
A fully explicit dominance test, and prune the given node if it's fathomed by dominance.
(i.e. ∀ l∈L: ∃ u∈U s.t. λu ≤ λl )
Return `true` if the given node is fathomed by dominance.
"""
function fullyExplicitDominanceTest(node::Node, UBS::Vector{SupportedSolutionPoint}, model)
    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(UBS) == 0 || length(node.lower_bound_set)==0 return false end

    p = MOI.output_dimension(model.f) ; nadir_pts = getNadirPoints(UBS, model)

    # ------------------------------------------
    # if the LBS consists of a single point
    # ------------------------------------------
    if length(node.lower_bound_set) == 1
        for u ∈ nadir_pts                   # if there exists an upper bound u s.t. u≦l
            if dominates(u, node.lower_bound_set[1])
                return true
            end
        end
        return false
    end

    UBS_ideal = UBS[1].y[:] ; LBS_ideal = node.lower_bound_set[1].y[:]

    for i in 2:length(UBS)
        for z in 1:p
            if UBS[i].y[z] < UBS_ideal[z] UBS_ideal[z] = UBS[i].y[z] end 
        end
    end
    for i in 2:length(node.lower_bound_set)
        for z in 1:p
            if node.lower_bound_set[i].y[z] < LBS_ideal[z] LBS_ideal[z] = node.lower_bound_set[i].y[z] end 
        end
    end
    UBS_ideal_sp = SupportedSolutionPoint(Set{Vector{Float64}}(), UBS_ideal, Vector{Float64}(), false)
    LBS_ideal_sp = SupportedSolutionPoint(Set{Vector{Float64}}(), LBS_ideal, Vector{Float64}(), false)

    # ----------------------------------------------
    # if the LBS consists of hyperplanes
    # ----------------------------------------------

    # if only one feasible point in UBS 
    if length(UBS) == 1 
        return dominates(UBS_ideal_sp, LBS_ideal_sp)
    end

    # test range condition necessary 1 : LBS ⊆ UBS (i.e. UBS includes the LP lexico-optimum)
    if !dominates( UBS_ideal_sp, LBS_ideal_sp)  return false end

    # test condition necessary 2 : UBS dominates LBS 
    # iterate of all local nadir points
    for u ∈ nadir_pts
        existence = false

        # case 1 : if u is dominates the ideal point of LBS 
        if dominates(u, LBS_ideal_sp)
            return true
        end

        # case 3 : complete pairwise comparison
        for l in node.lower_bound_set             # ∀ l ∈ LBS 
            if l.λ'*u.y < l.λ'*l.y         # todo : add TOL ? 
                existence = true ; break
            end
        end
        
        if !existence return false end
    end

    return true
end

