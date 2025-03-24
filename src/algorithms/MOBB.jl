include("MOBBTree.jl")

"""
   MultiObjectiveBranchBound()

`MultiObjectiveBranchBound` implements the multi-objective branch&bound framework.

## Supported optimizer attributes

* `MOA.LowerBoundsLimit()`: the maximum number of lower bounds calculated at each B&B node.

 ## Hypothesis :

 * only consider BINARY LINEAR programs for now (but not limited to) # todo : change branching strategy 

 * no not deal with objective with type `FEASIBILITY_SENSE`

"""
mutable struct MultiObjectiveBranchBound <: AbstractAlgorithm
    lowerbounds_limit::Union{Nothing,Int}                   # the number of lower bounds solved at each node 
    traverse_order :: Union{Nothing, Symbol}                # the traversing order of B&B tree
    tolerance :: Union{Nothing, Float64}                    # numerical tolerance
    convex_qcr :: Union{Nothing, Bool}                      # QCR convexification parameter 
    heuristic :: Union{Nothing, Bool}
    preproc :: Union{Nothing, Int64}                        # preprocessing choice -1, 0, 1, 2
    tight_root:: Union{Nothing, Int64}                       # calculate the tight LBS at root (0 : false, 1: QCR, 2 : DW)

    # --------------- informations for getting attributes 
    pruned_nodes :: Union{Nothing, Int64}
    pruned_dominance_nodes :: Union{Nothing, Int64}
    heuristic_time :: Union{Nothing, Float64}

    nb_vars :: Union{Nothing, Int64}
    Qs :: Union{Nothing, Vector{Matrix{Float64}}}
    Ls :: Union{Nothing, Vector{Vector{Float64}}}
    Cs :: Union{Nothing, Vector{Float64}}
    A_eq :: Union{Nothing, Matrix{Float64}}
    A_iq :: Union{Nothing, Matrix{Float64}}
    b_eq :: Union{Nothing, Vector{Float64}}
    b_iq :: Union{Nothing, Vector{Float64}}
    variables:: Union{Nothing, Vector{MOI.VariableIndex}}
    variableIndex :: Union{Nothing, Dict{MOI.VariableIndex, Int64}}
    indexVariable :: Union{Nothing, Dict{Int64, MOI.VariableIndex}}

    preproc_μ :: Union{Nothing, Vector{Vector{Vector{Float64}}}} # ∀ p obj, ∀ levels, μ

    MultiObjectiveBranchBound() = new(nothing, nothing, nothing, true, false, 0, 0,
                                      nothing, nothing, nothing,
                                      nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                                      nothing
                                )
end



# -------------------------------------
# ----------- parameters --------------
# -------------------------------------
MOI.supports(::MultiObjectiveBranchBound, ::LowerBoundsLimit) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::LowerBoundsLimit, value)
    alg.lowerbounds_limit = value ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::LowerBoundsLimit)
    return something(alg.lowerbounds_limit, default(alg, attr))
end


MOI.supports(::MultiObjectiveBranchBound, ::Heuristic) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::Heuristic, value)
    alg.heuristic = value ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::Heuristic)
    return something(alg.heuristic, default(alg, attr))
end


MOI.supports(::MultiObjectiveBranchBound, ::TraverseOrder) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::TraverseOrder, order)
    alg.traverse_order = order ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::TraverseOrder)
    return something(alg.traverse_order, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::Tolerance) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::Tolerance, tol)
    alg.tolerance = tol ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::Tolerance)
    return something(alg.tolerance, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::ConvexQCR) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::ConvexQCR, value)
    alg.convex_qcr = value ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::ConvexQCR)
    return something(alg.convex_qcr, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::Preproc) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::Preproc, value)
    alg.preproc = value ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::Preproc)
    return something(alg.preproc, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::TightRoot) = true

function MOI.set(alg::MultiObjectiveBranchBound, ::TightRoot, value)
    alg.tight_root = value ; return
end

function MOI.get(alg::MultiObjectiveBranchBound, attr::TightRoot)
    return something(alg.tight_root, default(alg, attr))
end

# --------- attributes only for getting 
MOI.supports(::MultiObjectiveBranchBound, ::PrunedNodeCount) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::PrunedNodeCount)
    return something(alg.pruned_nodes, default(alg, attr))
end


# --------- attributes only for getting 
MOI.supports(::MultiObjectiveBranchBound, ::PrunedDominanceNodeCount) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::PrunedDominanceNodeCount)
    return something(alg.pruned_dominance_nodes, default(alg, attr))
end



MOI.supports(::MultiObjectiveBranchBound, ::HeuristicTime) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::HeuristicTime)
    return something(alg.heuristic_time, default(alg, attr))
end


MOI.supports(::MultiObjectiveBranchBound, ::NBvars) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::NBvars)
    return something(alg.nb_vars, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::QObj) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::QObj)
    return something(alg.Qs, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::LObj) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::LObj)
    return something(alg.Ls, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::CObj) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::CObj)
    return something(alg.Cs, default(alg, attr))
end


MOI.supports(::MultiObjectiveBranchBound, ::Aeq) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::Aeq)
    return something(alg.A_eq, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::Aiq) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::Aiq)
    return something(alg.A_iq, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::Beq) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::Beq)
    return something(alg.b_eq, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::Biq) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::Biq)
    return something(alg.b_iq, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::PreprocMu) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::PreprocMu)
    return something(alg.preproc_μ, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::Variables) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::Variables)
    return something(alg.variables, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::VariablesIndex) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::VariablesIndex)
    return something(alg.variableIndex, default(alg, attr))
end

MOI.supports(::MultiObjectiveBranchBound, ::IndexVariables) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::IndexVariables)
    return something(alg.indexVariable, default(alg, attr))
end

"""
    Load coefficients matrices of the model. 
"""
function _loadMatrices(algorithm::MultiObjectiveBranchBound, model::Optimizer)
    N = algorithm.nb_vars
    varIndex = algorithm.variableIndex

    # read objective Qs # todo : only for multiple quadratic functions !! 
    if typeof(model.f) == MOI.VectorQuadraticFunction{Float64}
        algorithm.Qs = Vector{Matrix{Float64}}() ; algorithm.Ls = Vector{Vector{Float64}}() ; algorithm.Cs = Vector{Float64}()
        for _ in 1:MOI.output_dimension(model.f) 
            push!(algorithm.Qs, zeros(N, N) ) ; push!(algorithm.Ls, zeros(N)) ; push!(algorithm.Cs, 0.0)
        end

        for term in model.f.quadratic_terms
            i = varIndex[term.scalar_term.variable_1 ]; j = varIndex[term.scalar_term.variable_2 ]
            algorithm.Qs[term.output_index][i, j] = i==j ? term.scalar_term.coefficient/2 : term.scalar_term.coefficient
        end

        for term in model.f.affine_terms
            i = varIndex[term.scalar_term.variable ] 
            algorithm.Ls[term.output_index][i] = term.scalar_term.coefficient
        end

        algorithm.Cs = model.f.constants
          
    end

    # read constraints ax=b     # todo :  VectorAffineFunction case 
    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}())
    algorithm.A_eq = zeros(length(ctr_t), N) ; algorithm.b_eq = zeros(length(ctr_t))
    i = 0
    for ci in ctr_t
        i +=1
        f_ = MOI.get(model.inner, MOI.ConstraintFunction(), ci) 

        for term in f_.terms
            j = varIndex[term.variable]
            algorithm.A_eq[i, j] = term.coefficient
        end
        algorithm.b_eq[i] = MOI.get(model.inner, MOI.ConstraintSet(), ci).value
    end

    # read constraint ax <= b   # todo :  VectorAffineFunction
    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}())
    algorithm.A_iq = zeros(length(ctr_t), N) ; algorithm.b_iq = zeros(length(ctr_t))
    i = 0
    for ci in ctr_t
        i +=1
        f_ = MOI.get(model.inner, MOI.ConstraintFunction(), ci) 

        for term in f_.terms
            j = varIndex[term.variable]
            algorithm.A_iq[i, j] = term.coefficient
        end
        algorithm.b_iq[i] = MOI.get(model.inner, MOI.ConstraintSet(), ci).upper
    end
    
end

"""
    Relax binary variables to continous between 0.0 and 1.0.
"""
function relaxVariables(model::Optimizer) :: Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}}
    vars_idx = MOI.get(model, MOI.ListOfVariableIndices())

    # todo : to complete contraints type see  https://jump.dev/MathOptInterface.jl/stable/manual/constraints/
    for (t1, t2) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        ctr_t =  MOI.get(model, MOI.ListOfConstraintIndices{t1,t2}())

        if t1 == MOI.VariableIndex
            for ci in ctr_t
                MOI.delete(model, ci)
            end
        end
    end

    lower_bounds = Dict{MOI.VariableIndex, MOI.ConstraintIndex}(
        var => MOI.add_constraint(model, var, MOI.GreaterThan(0.0)) for var in vars_idx
    )
    
    upper_bounds = Dict{MOI.VariableIndex, MOI.ConstraintIndex}(
        var => MOI.add_constraint(model, var, MOI.LessThan(1.0)) for var in vars_idx
    )

    return[lower_bounds, upper_bounds]
end


function _preprocessing_UQCR(model, algorithm::MultiObjectiveBranchBound)
    p = MOI.output_dimension(model.f) ; algorithm.preproc_μ = Vector{Vector{Vector{Float64}}}()

    # level 0 (all vars are free)
    for p_ in 1:p
        μ = klevel_UQCR_csdp(algorithm.nb_vars, algorithm.Qs[p_], 0, algorithm)
        if μ === nothing @error("UQCR error in preprocessing at level 0 for objective $p_ ! ") end 
        push!(algorithm.preproc_μ, [μ])
    end

    for k in 1:algorithm.nb_vars-1
        for p_ in 1:p
            μ = klevel_UQCR_csdp(algorithm.nb_vars - k, algorithm.Qs[p_][k+1:end, k+1:end], k, algorithm)
            if μ === nothing @error("UQCR error in preprocessing at level $k for objective $p_ ! ") end 
            push!(algorithm.preproc_μ[p_], μ)
        end
    end
end

function MOBB(algorithm::MultiObjectiveBranchBound, model::Optimizer, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}},
            tree, node::Node, UBS::Vector{SupportedSolutionPoint}, LBS::Vector{SupportedSolutionPoint}
)

    # println("\n\n -------------- node $(node.num) ")


    # get the actual node
    @assert node.activated == true "the actual node is not activated "
    node.activated = false ;  node.assignment = getPartialAssign(node)

    # calculate the lower bound set 
    if computeLBS(node, model, algorithm, Bounds)
        prune!(node, INFEASIBILITY) ; algorithm.pruned_nodes += 1
        # println(node)
        return
    end

    # if node.depth <= algorithm.nb_vars/2 # !isRoot(node) && node.var_bound == 1.0 # && 
    #     # before = length(UBS)
    #     # varArray = MOI.get(model, MOI.ListOfVariableIndices())
    #     # varidx = Dict(varArray[i] => i for i in 1:algorithm.nb_vars)
    #     # algorithm.heuristic_time += heuristic_local(model, UBS, algorithm, node.assignment, varidx)
    #     # # println("UBS changed ", length(UBS) - before)

    #     start = time()
    #     p = MOI.output_dimension(model.f)
    #     for sol in node.lower_bound_set
    #         if !sol.is_integer
    #             x = [ v>0.9 ? 1 : 0 for v in first(sol.x)]
    #             if sum(x) == algorithm.b_eq[1] && x'* algorithm.A_iq[1, :] < algorithm.b_iq[1]
    #                 push_filtering_dominance(UBS, SupportedSolutionPoint(Set([x]),
    #                                                     [x'*algorithm.Qs[k]'*x + x'*algorithm.Ls[k] + algorithm.Cs[k] for k in 1:p],
    #                                                     sol.λ, true
    #                                     ) 
    #                 )
    #             end
    #         end
    #     end
    #     algorithm.heuristic_time += time() - start

    # end
    
    # update the upper bound set 
    if updateUBS(node, UBS)
        prune!(node, INTEGRITY) ; algorithm.pruned_nodes += 1 
        # println(node)
        return 
    end 
   
    newLBS = Vector{SupportedSolutionPoint}()

    # todo : global LBS intersection 
    if MOI.get(algorithm, TightRoot()) > 0 
        if isRoot(node)
            for p in node.lower_bound_set push!(LBS, p) end 
            newLBS = node.lower_bound_set
        else
            
            dominated = [false for _ in node.lower_bound_set] ; i = 0
            for p in node.lower_bound_set
                i += 1  
                for l in LBS
                    if p.y'*l.λ < l.y'*l.λ
                        dominated[i] = true ; break
                    end
                end
            end
            for i in 1:length(node.lower_bound_set)
                if !dominated[i]
                    push!(newLBS, node.lower_bound_set[i])
                end
            end

            dominated = [false for _ in LBS] ; i = 0
            for l in LBS
                i += 1  
                for p in node.lower_bound_set
                    if l.y'*p.λ < p.y'*p.λ
                        dominated[i] = true ; break
                    end
                end
            end
            for i in 1:length(LBS)
                if !dominated[i]
                    push!(newLBS, LBS[i])
                end
            end

        end
    else
        newLBS = node.lower_bound_set
    end

    # println("newLBS = ", newLBS)
    # test dominance 
    if fullyExplicitDominanceTest(newLBS, UBS, model)
        prune!(node, DOMINANCE) ; algorithm.pruned_nodes += 1 ; algorithm.pruned_dominance_nodes += 1
        # println(node)
        return
    end


    # todo :: for k-item knapsack instance only 
    if length(algorithm.b_eq)>0 && sum(collect(values(node.assignment))) >= algorithm.b_eq[1] return end 


    # otherwise this node is not fathomed, continue to branch on free variable
    var = pickUpAFreeVar(node.assignment, model)
    if var === nothing return end

    children =[ Node(model.total_nodes + 1, node.depth + 1, pred = node, var_idx = var, var_bound = 1.0, bound_type = 2),
                Node(model.total_nodes + 2, node.depth + 1, pred = node, var_idx = var, var_bound = 0.0, bound_type = 1)
    ]
    for child in children
        addTree(tree, algorithm, child) ; model.total_nodes += 1
        if MOI.get(algorithm, Preproc()) == 0 child.qcr_coeff = node.qcr_coeff end 
        push!(node.succs, child)
    end

    # println(node)
end

# -------------------------------------
# ----------- main program ------------
# -------------------------------------

function optimize_multiobjective!(
    algorithm::MultiObjectiveBranchBound,
    model::Optimizer,
    # verbose :: Bool = false,
)
    model.total_nodes = 0 ; algorithm.pruned_nodes = 0 ; algorithm.pruned_dominance_nodes = 0 
    start_time = time()

    # step1 - set tolerance to inner model 
    if MOI.get(algorithm, Tolerance()) != default(algorithm, Tolerance())
        MOI.set(model, MOI.RawOptimizerAttribute("tol_inconsistent"), MOI.get(algorithm, Tolerance()))
    end

    # step2 - check lower bounds limit 
    if MOI.get(algorithm, LowerBoundsLimit()) <= MOI.output_dimension(model.f)
        # at least p lower bounds optimized on each objective 
        MOI.set(algorithm, LowerBoundsLimit(), MOI.output_dimension(model.f) +1 )
    end

    # step3 - LP relaxation 
    Bounds = relaxVariables(model)
    algorithm.variables = MOI.get(model, MOI.ListOfVariableIndices())
    algorithm.nb_vars = length(algorithm.variables)
    algorithm.variableIndex = Dict(algorithm.variables[i] => i for i=1:algorithm.nb_vars)
    algorithm.indexVariable = Dict(i => algorithm.variables[i] for i=1:algorithm.nb_vars)
    
    _loadMatrices(algorithm, model)

    # preprocessing phase
    if MOI.get(algorithm, Preproc()) > 0
        _preprocessing_UQCR(model, algorithm)
    end

    # step4 - initialization UBS
    UBS = Vector{SupportedSolutionPoint}()
    # global heuristic 
    if MOI.get(algorithm, Heuristic())
        algorithm.heuristic_time = heuristic(model, UBS, algorithm)
    end

    LBS = Vector{SupportedSolutionPoint}()

     
    tree = initTree(algorithm)
    model.total_nodes += 1 ; root = Node(model.total_nodes, 0)
    addTree(tree, algorithm, root)

    status = MOI.OPTIMAL

    # step5 - study every node in tree 
    while length(tree) > 0
        if _time_limit_exceeded(model, start_time)
            status = MOI.TIME_LIMIT
            break
        end

        node_ref = nextNodeTree(tree, algorithm)

        MOBB(algorithm, model, Bounds, tree, node_ref[], UBS, LBS)
        
        if node_ref[].deleted
            finalize(node_ref[])
        end
    end
    
    vars_idx = MOI.get(model, MOI.ListOfVariableIndices())

    return status, [SolutionPoint(
                                    Dict(vars_idx[i] => first(sol.x)[i] for i in 1:length(vars_idx) ) , sol.y
                    ) 
                    for sol in UBS
                ]
end
