"""
   MultiObjectiveBranchBound()

`MultiObjectiveBranchBound` implements the multi-objective branch&bound framework.

## Supported optimizer attributes

* `MOA.LowerBoundsLimit()`: the maximum number of lower bounds calculated at each B&B node.

 ## Hypothesis :

 * only consider BINARY LINEAR programs for now (but not limited to)

 * no not deal with objective with type `FEASIBILITY_SENSE`

"""

mutable struct MultiObjectiveBranchBound <: AbstractAlgorithm
    lowerbounds_limit::Union{Nothing,Int}                   # the number of lower bounds solved at each node 
    
    MultiObjectiveBranchBound() = new(nothing)
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

# todo : add traverse attribute 

# -------------------------------------
# ----------- main program ------------
# -------------------------------------

function optimize_multiobjective!(
    algorithm::MultiObjectiveBranchBound,
    model::Optimizer,
    # verbose :: Bool = false,
)
    println("welcom to MOBB ...")
    start_time = time()

    # --------------------------------
    # -------- preamble --------------
    # --------------------------------

    # step1 - transfer Max to Min function 
    MAX_obj_reversed = false
    if MOI.get(model.inner, MOI.ObjectiveSense()) == MOI.MAX_SENSE
        MAX_obj_reversed = true ; MOI.set(model.inner, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        model.f = model.f * -1 
    end 

    # step2 - check lower bounds limit 
    if MOI.get(algorithm, LowerBoundsLimit()) < MOI.output_dimension(model.f)
        # at least p lower bounds optimized on each objective 
        MOI.set(algorithm, LowerBoundsLimit(), MOI.output_dimension(model.f) +1 )
    end


    return nothing, []
end
