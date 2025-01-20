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

    # --------------- informations for getting attributes 
    pruned_nodes :: Union{Nothing, Int64}
    heuristic_time :: Union{Nothing, Float64}

    nb_vars :: Union{Nothing, Int64}
    Qs :: Union{Nothing, Vector{Matrix{Float64}}}
    A_eq :: Union{Nothing, Matrix{Float64}}
    A_iq :: Union{Nothing, Matrix{Float64}}
    b_eq :: Union{Nothing, Vector{Float64}}
    b_iq :: Union{Nothing, Vector{Float64}}

    MultiObjectiveBranchBound() = new(nothing, nothing, nothing, nothing, nothing,
                                      nothing, nothing, 
                                      nothing, nothing, nothing, nothing, nothing, nothing
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

# --------- attributes only for getting 
MOI.supports(::MultiObjectiveBranchBound, ::PrunedNodeCount) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::PrunedNodeCount)
    return something(alg.pruned_nodes, default(alg, attr))
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

"""
    Load coefficients matrices of the model. 
"""
function _loadMatrices(algorithm::MultiObjectiveBranchBound, model::Optimizer)
    varArray = MOI.get(model.inner, MOI.ListOfVariableIndices())
    N = length(varArray) ; algorithm.nb_vars = N 
    varIndex = Dict(varArray[i] => i for i=1:N)

    # read objective Qs # todo : only for multiple quadratic functions !! 
    if typeof(model.f) == MOI.VectorQuadraticFunction{Float64}
        algorithm.Qs = Vector{Matrix{Float64}}()
        for _ in 1:MOI.output_dimension(model.f) 
            push!(algorithm.Qs, zeros(N, N) )
        end

        for term in model.f.quadratic_terms
            i = varIndex[term.scalar_term.variable_1 ]; j = varIndex[term.scalar_term.variable_2 ]
            algorithm.Qs[term.output_index][i, j] = i==j ? term.scalar_term.coefficient/2 : term.scalar_term.coefficient
        end
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


function MOBB(algorithm::MultiObjectiveBranchBound, model::Optimizer, Bounds::Vector{Dict{MOI.VariableIndex, MOI.ConstraintIndex}},
            tree, node::Node, UBS::Vector{SupportedSolutionPoint}
)

    # println("\n\n -------------- node $(node.num) ")

    # get the actual node
    @assert node.activated == true "the actual node is not activated "
    node.activated = false

    # calculate the lower bound set 
    if computeLBS(node, model, algorithm, Bounds)
        prune!(node, INFEASIBILITY) ; algorithm.pruned_nodes += 1
        return
    end

    # update the upper bound set 
    if updateUBS(node, UBS)
        algorithm.pruned_nodes += 1 ; return 
    end 
    

    # test dominance 
    if fullyExplicitDominanceTest(node, UBS, model)
        prune!(node, DOMINANCE) ; algorithm.pruned_nodes += 1
        return
    end

    # otherwise this node is not fathomed, continue to branch on free variable
    assignment = getPartialAssign(node) ; var = pickUpAFreeVar(assignment, model)
    if var === nothing return end

    children =[ Node(model.total_nodes + 1, node.depth + 1, pred = node, var_idx = var, var_bound = 1.0, bound_type = 2),
                Node(model.total_nodes + 2, node.depth + 1, pred = node, var_idx = var, var_bound = 0.0, bound_type = 1)
    ]
    for child in children
        addTree(tree, algorithm, child) ; model.total_nodes += 1
        push!(node.succs, child)
    end

end

# -------------------------------------
# ----------- main program ------------
# -------------------------------------

function optimize_multiobjective!(
    algorithm::MultiObjectiveBranchBound,
    model::Optimizer,
    # verbose :: Bool = false,
)
    model.total_nodes = 0 ; algorithm.pruned_nodes = 0
    start_time = time()
    # step1 - set tolerance to inner model 
    if MOI.get(algorithm, Tolerance()) != default(algorithm, Tolerance())
        MOI.set(model, MOI.RawOptimizerAttribute("tol_inconsistent"), MOI.get(algorithm, Tolerance()))
    end

    # step2 - check lower bounds limit 
    if MOI.get(algorithm, LowerBoundsLimit()) < MOI.output_dimension(model.f)
        # at least p lower bounds optimized on each objective 
        MOI.set(algorithm, LowerBoundsLimit(), MOI.output_dimension(model.f) +1 )
    end

    # # -------------------------------
    # println(model)
    # # -------------------------------

    # step3 - LP relaxation 
    Bounds = relaxVariables(model)

    # -------------------------------
    println(model)
    # -------------------------------

    _loadMatrices(algorithm, model)

    # step4 - initialization
    UBS = Vector{SupportedSolutionPoint}()
    # global heuristic 
    if MOI.get(algorithm, Heuristic())
        algorithm.heuristic_time = heuristic(model, UBS, algorithm)
    end
    println("UBS = ", UBS)

     
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

        MOBB(algorithm, model, Bounds, tree, node_ref[], UBS)
        
        if node_ref[].deleted
            finalize(node_ref[])
        end
    end
    
    vars_idx = MOI.get(model, MOI.ListOfVariableIndices())
    # todo : tol rounding 
    return status, [SolutionPoint(
                                    Dict(vars_idx[i] => first(sol.x)[i] for i in 1:length(vars_idx) ) , sol.y
                    ) 
                    for sol in UBS
                ]
end
