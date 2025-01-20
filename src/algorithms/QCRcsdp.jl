using JuMP, LinearAlgebra, CSDP
include("MOBB.jl")


function QCR_csdp(Q, c, constant, model::Optimizer, varArray, varIndex, 
                                    algorithm::MultiObjectiveBranchBound, qcr_coeff::QCRcoefficients 
    )
    N = length(varArray)

    # SDP model ...
    model_sdp = Model(CSDP.Optimizer) ; JuMP.set_silent(model_sdp)
    set_attribute(model_sdp, "axtol", 1.0e-5) ; set_attribute(model_sdp, "atytol", 1.0e-5)

    @variable(model_sdp, x[1:N] )
    sdpIndex = Dict(x[i] => i for i=1:N)

    @variable(model_sdp, X[1:N, 1:N], Symmetric)
    @constraint(model_sdp, [1 x'; x X] in PSDCone())

    M = length(algorithm.b_eq)

    @objective(model_sdp, Min, tr(Q * X) + c'*x + constant )
    @constraint(model_sdp, algorithm.A_eq *x == algorithm.b_eq) 
    @constraint(model_sdp, algorithm.A_iq *x <= algorithm.b_iq) 

    con_μ = @constraint(model_sdp, [i in 1:N], X[i,i] - x[i] == 0)
    con_β = @constraint(model_sdp, sum( sum(algorithm.A_eq[k,i] * algorithm.A_eq[k,j]*X[i,j] for i in 1:N for j in 1:N) - 
                                    2 * sum(algorithm.b_eq[k] * algorithm.A_eq[k,j] * x[j] for j in 1:N) +
                                    algorithm.b_eq[k]^2  for k in 1:M
                                ) == 0 
            ) 

    optimize!(model_sdp)

    if termination_status(model_sdp)== OPTIMAL
        # obj = objective_value(model_sdp)

        μ = dual.(con_μ)
        β = dual(con_β)

        @expression(model_sdp, newf, x'*Q*x + c'*x + constant +
                                sum((-μ[i])* (x[i]^2 - x[i]) for i in 1:N) -
                                β * sum( (algorithm.A_eq[k, :]'*x - algorithm.b_eq[k] )^2 for k in 1:M)
        )
        qcr_coeff.Q = zeros(N, N) ; qcr_coeff.c = zeros(N) ; qcr_coeff.constant = newf.aff.constant

        for (p, v) in newf.terms
            i = sdpIndex[p.a] ; j = sdpIndex[p.b]
            if i == j 
                qcr_coeff.Q[i,j] = v
            else
                qcr_coeff.Q[i, j] = v/2 ; qcr_coeff.Q[j, i] = v/2
            end
        end

        for (k, v) in newf.aff.terms
            i = sdpIndex[k]
            qcr_coeff.c[i] = v
        end

        return true

    elseif termination_status(model_sdp) == MOI.ITERATION_LIMIT
            @warn("CSDP solver 100 iter limit reached ! ")
    end

    return false
end



"""
x'*Q*x + c'*x 
"""
function solve_weighted_sum(
    model::Optimizer,
    weights::Vector{Float64},
    QCR::Bool,
    algorithm::MultiObjectiveBranchBound, 
    qcr_coeff::QCRcoefficients 
)
    f = _scalarise(model.f, weights)
    MOI.set(model.inner, MOI.ObjectiveFunction{typeof(f)}(), f)

    # verify whether quadratic function is convex 
    if MOI.get(model.inner, MOI.ObjectiveFunctionType()) == MOI.ScalarQuadraticFunction{Float64}
        varArray = MOI.get(model.inner, MOI.ListOfVariableIndices())
        N = length(varArray)
        varIndex = Dict(varArray[i] => i for i=1:N)
        Q = zeros(N, N) ; c = zeros(N)
    
        for term in f.quadratic_terms
            if abs(term.coefficient) != 0.0
                i = varIndex[term.variable_1]; j = varIndex[term.variable_2]
                Q[i, j] = i==j ? term.coefficient/2 : term.coefficient
            end
        end

        for term in f.affine_terms
            i = varIndex[term.variable]
            c[i] = term.coefficient
        end
    
 
        if QCR 
            is_solved = QCR_csdp(Q, c, f.constant, model, varArray, varIndex, algorithm, qcr_coeff)

            if !is_solved
                return MOI.INFEASIBLE, nothing
            else
                quad_terms = MOI.ScalarQuadraticTerm{Float64}[
                    MOI.ScalarQuadraticTerm(
                        i==j ? qcr_coeff.Q[i, j]*2+0.0001 : qcr_coeff.Q[i, j],
                        varArray[i],
                        varArray[j],
                    ) for i in 1:N for j in 1:N
                ]

                affine_terms = MOI.ScalarAffineTerm{Float64}[
                    MOI.ScalarAffineTerm(
                        qcr_coeff.c[i],
                        varArray[i],
                    ) for i in 1:N 
                ]  
    
                QCR_f = MOI.ScalarQuadraticFunction(quad_terms, affine_terms, qcr_coeff.constant)
                MOI.set(model.inner, MOI.ObjectiveFunction{typeof(QCR_f)}(), QCR_f )
            end

        end
    end

    MOI.optimize!(model.inner)

    status = MOI.get(model.inner, MOI.TerminationStatus())
    if !_is_scalar_status_optimal(status)
        return status, nothing
    end
    variables = MOI.get(model.inner, MOI.ListOfVariableIndices())
    X, Y = _compute_point(model, variables, model.f)
    return status, SolutionPoint(X, Y)
end

