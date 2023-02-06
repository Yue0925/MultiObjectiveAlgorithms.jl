#  Copyright 2019, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v.2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

"""
    EpsilonConstraint()

`EpsilonConstraint` implements the epsilon-constraint algorithm for
bi-objective programs.

## Supported optimizer attributes

 * `MOO.SolutionLimit()`: with a slight abuse of notation, `EpsilonConstraint`
   divides the width of the first-objective's domain in objective space by
   `SolutionLimit` to obtain the epsilon to use when iterating. Thus, there
   can be at most `SolutionLimit` solutions returned, but there may be fewer.
   If no value is set, the default is `100`, instead of the typical
   `default(::SolutionLimit)`.
"""
mutable struct EpsilonConstraint <: AbstractAlgorithm
    solution_limit::Union{Nothing,Int}

    EpsilonConstraint() = new(nothing)
end

default(::EpsilonConstraint, ::SolutionLimit) = 100

MOI.supports(::EpsilonConstraint, ::SolutionLimit) = true

function MOI.set(alg::EpsilonConstraint, ::SolutionLimit, value)
    alg.solution_limit = value
    return
end

function MOI.get(alg::EpsilonConstraint, attr::SolutionLimit)
    return something(alg.solution_limit, default(alg, attr))
end

function optimize_multiobjective!(
    algorithm::EpsilonConstraint,
    model::Optimizer,
)
    if MOI.output_dimension(model.f) != 2
        error("EpsilonConstraint requires exactly two objectives")
    end
    # Compute the bounding box ofthe objectives using Hierarchical().
    alg = Hierarchical()
    MOI.set.(Ref(alg), ObjectivePriority.(1:2), [1, 0])
    status, solution_1 = optimize_multiobjective!(alg, model)
    @assert status == MOI.OPTIMAL
    MOI.set(alg, ObjectivePriority(2), 2)
    status, solution_2 = optimize_multiobjective!(alg, model)
    @assert status == MOI.OPTIMAL
    a, b = solution_1[1].y[1], solution_2[1].y[1]
    left, right = min(a, b), max(a, b)
    # Compute the epsilon that we will be incrementing by each iteration
    n_points = MOI.get(algorithm, SolutionLimit())
    ε = abs(right - left) / (n_points - 1)
    solutions = ParetoSolution[]
    f1, f2 = MOI.Utilities.eachscalar(model.f)
    MOI.set(model.inner, MOI.ObjectiveFunction{typeof(f2)}(), f2)
    # Add epsilon constraint
    SetType = ifelse(
        MOI.get(model.inner, MOI.ObjectiveSense()) == MOI.MIN_SENSE,
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
    )
    ci = MOI.add_constraint(model, f1, SetType(left + ε))
    variables = MOI.get(model.inner, MOI.ListOfVariableIndices())
    for i in 1:n_points
        rhs = left + (i - 1) * ε
        MOI.set(model, MOI.ConstraintSet(), ci, SetType(rhs))
        MOI.optimize!(model.inner)
        if MOI.get(model.inner, MOI.TerminationStatus()) != MOI.OPTIMAL
            return MOI.OTHER_ERROR, nothing
        end
        X = Dict{MOI.VariableIndex,Float64}(
            x => MOI.get(model.inner, MOI.VariablePrimal(), x) for
            x in variables
        )
        Y = MOI.Utilities.eval_variables(x -> X[x], model.f)
        if isempty(solutions) || !(Y ≈ solutions[end].y)
            push!(solutions, ParetoSolution(X, Y))
        end
    end
    MOI.delete(model, ci)
    return MOI.OPTIMAL, unique(solutions)
end