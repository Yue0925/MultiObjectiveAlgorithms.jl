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

    # --------------- informations for getting attributes 
    pruned_nodes :: Union{Nothing, Int64}

    MultiObjectiveBranchBound() = new(nothing, nothing, nothing,
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

# --------- attributes only for getting 
MOI.supports(::MultiObjectiveBranchBound, ::PrunedNodeCount) = true

function MOI.get(alg::MultiObjectiveBranchBound, attr::PrunedNodeCount)
    return something(alg.pruned_nodes, default(alg, attr))
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
            # -------------------------------
            # println("t1 ", t1 , "  , t2 ", t2 )
            # println("ctr ", ctr_t)
            # -------------------------------

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

    println("\n\n -------------- node $(node.num) ")


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


    # step4 - initialization
    UBS = Vector{SupportedSolutionPoint}()
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
                                    Dict(vars_idx[i] => sol.x[1][i] for i in 1:length(vars_idx) ) , sol.y
                    ) 
                    for sol in UBS
                ]
end
