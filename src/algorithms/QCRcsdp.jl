using JuMP, LinearAlgebra, CSDP
include("MOBB.jl")


function QCR_csdp(Q, c, constant, model,
                                algorithm::MultiObjectiveBranchBound, qcr_coeff::QCRcoefficients 
    )

    N = algorithm.nb_vars
    varIndex = algorithm.variableIndex

    # SDP model ...
    model_sdp = Model(CSDP.Optimizer) ; JuMP.set_silent(model_sdp) 
    set_attribute(model_sdp, "axtol", 1.0e-4) ; set_attribute(model_sdp, "atytol", 1.0e-4)
    set_attribute(model_sdp, "maxiter", 200) 

    @variable(model_sdp, x[1:N] )
    sdpIndex = Dict(x[i] => i for i=1:N)

    @variable(model_sdp, X[1:N, 1:N], Symmetric)
    @constraint(model_sdp, [1 x'; x X] in PSDCone())

    # variable assignment
    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.GreaterThan{Float64}}())
    for ci in ctr_t
        i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
        @constraint(model_sdp, x[i] >= MOI.get(model.inner, MOI.ConstraintSet(), ci).lower )
    end

    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.LessThan{Float64}}())
    for ci in ctr_t
        i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
        @constraint(model_sdp, x[i] <= MOI.get(model.inner, MOI.ConstraintSet(), ci).upper )
    end

    M = length(algorithm.b_eq)

    @objective(model_sdp, Min, tr(Q * X)+ c'*x + constant  ) 
    length(algorithm.b_eq) > 0 ? @constraint(model_sdp, algorithm.A_eq *x == algorithm.b_eq) : nothing
    length(algorithm.b_iq) > 0 ? @constraint(model_sdp, algorithm.A_iq *x <= algorithm.b_iq) : nothing

    con_μ = @constraint(model_sdp, [i in 1:N], X[i,i] - x[i] == 0)

    con_β = nothing
    if length(algorithm.b_eq) > 0 
        con_β = @constraint(model_sdp, sum( sum(algorithm.A_eq[k,i] * algorithm.A_eq[k,j]*X[i,j] for i in 1:N for j in 1:N) - 
                                        2 * sum(algorithm.b_eq[k] * algorithm.A_eq[k,j] * x[j] for j in 1:N) +
                                        algorithm.b_eq[k]^2  for k in 1:M
                                    ) == 0 
                ) 
    end

    optimize!(model_sdp)

    if termination_status(model_sdp)== OPTIMAL
        μ = dual.(con_μ)

        if length(algorithm.b_eq) > 0 
            β = dual(con_β)

            @expression(model_sdp, newf, x'*Q*x + c'*x + constant +
                                    sum((-μ[i])* (x[i]^2 - x[i]) for i in 1:N) -
                                    β * sum( (algorithm.A_eq[k, :]'*x - algorithm.b_eq[k] )^2 for k in 1:M)
            )
        else
            @expression(model_sdp, newf, x'*Q*x + c'*x + constant +
                                    sum((-μ[i])* (x[i]^2 - x[i]) for i in 1:N)
            )
        end


        Q = zeros(N, N) ; c = zeros(N) ; constant = newf.aff.constant

        for (p, v) in newf.terms
            i = sdpIndex[p.a] ; j = sdpIndex[p.b]
            if i == j 
                Q[i,j] = v
            else
                Q[i, j] = v/2 ; Q[j, i] = v/2
            end
        end

        for (k, v) in newf.aff.terms
            i = sdpIndex[k]
            c[i] = v
        end

        push!(qcr_coeff.Q, Q) ; push!(qcr_coeff.c, c) ; push!(qcr_coeff.constant, constant)

        return true

    elseif termination_status(model_sdp) == MOI.ITERATION_LIMIT
            @warn("CSDP solver 200 iter limit reached ! ")
    end

    return false
end



"""
x'*Q*x + c'*x +constant
"""
function solve_weighted_sum(
    model::Optimizer,
    λ::Vector{Float64},
    QCR::Bool,
    qcr_coeff::QCRcoefficients,
    algorithm :: MultiObjectiveBranchBound
    ) 

    varArray = algorithm.variables
    N = length(varArray)
    varIndex = algorithm.variableIndex

    # --------------------------------------------
    # -- by defaut 
    if !QCR 
        f = _scalarise(model.f, λ)
        MOI.set(model.inner, MOI.ObjectiveFunction{typeof(f)}(), f)

        MOI.optimize!(model.inner)
    
        status = MOI.get(model.inner, MOI.TerminationStatus())
        if !_is_scalar_status_optimal(status)
            return status, nothing, nothing
        end
    
        X, Y = _compute_point(model, varArray, model.f)
        sol = zeros(N)
        for (x, val) in X
            sol[varIndex[x]] = val
        end
        return status, sol, Y
    end    

    # ----------------------------------------------
    # -- weighted sum of convexde coefficients
    Q = sum( λ[p].* qcr_coeff.Q[p] for p in 1:length(λ))
    c = sum( λ[p].* qcr_coeff.c[p] for p in 1:length(λ) )
    constant = λ'* qcr_coeff.constant

    quad_terms = MOI.ScalarQuadraticTerm{Float64}[
        MOI.ScalarQuadraticTerm(
            i==j ? Q[i, j]*2 +0.0001 : Q[i, j],
            varArray[i],
            varArray[j],
        ) for i in 1:N for j in 1:N
    ]

    affine_terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(
            c[i],
            varArray[i],
        ) for i in 1:N 
    ]  

    QCR_f = MOI.ScalarQuadraticFunction(quad_terms, affine_terms, constant)
    MOI.set(model.inner, MOI.ObjectiveFunction{typeof(QCR_f)}(), QCR_f )

    MOI.optimize!(model.inner)

    status = MOI.get(model.inner, MOI.TerminationStatus())
    if !_is_scalar_status_optimal(status)
        return status, nothing, nothing
    end
    
    sol = zeros(N)
    for x in varArray
        sol[varIndex[x]] = MOI.get(model.inner, MOI.VariablePrimal(), x)
    end

    Y = [sol'* qcr_coeff.Q[p] *sol + sol'* qcr_coeff.c[p] + qcr_coeff.constant[p] for p in 1:length(λ)]

    return status, sol, Y

        
end




function QCR_independent(Q, model, algorithm::MultiObjectiveBranchBound )
    varArray = MOI.get(model.inner, MOI.ListOfVariableIndices())
    N = length(varArray)
    varIndex = Dict(varArray[i] => i for i=1:N)

    # SDP model ...
    model_sdp = Model(CSDP.Optimizer) ; JuMP.set_silent(model_sdp) 
    set_attribute(model_sdp, "axtol", 1.0e-5) ; set_attribute(model_sdp, "atytol", 1.0e-5)
    set_attribute(model_sdp, "maxiter", 200) 

    @variable(model_sdp, x[1:N] )
    sdpIndex = Dict(x[i] => i for i=1:N)

    @variable(model_sdp, X[1:N, 1:N], Symmetric)
    @constraint(model_sdp, [1 x'; x X] in PSDCone())

    # variable assignment
    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.GreaterThan{Float64}}())
    for ci in ctr_t
    i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
    @constraint(model_sdp, x[i] >= MOI.get(model.inner, MOI.ConstraintSet(), ci).lower )
    end

    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.LessThan{Float64}}())
    for ci in ctr_t
    i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
    @constraint(model_sdp, x[i] <= MOI.get(model.inner, MOI.ConstraintSet(), ci).upper )
    end

    M = length(algorithm.b_eq)

    @objective(model_sdp, Min, tr(Q * X) )
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
        μ = dual.(con_μ)
        β = dual(con_β)

        @expression(model_sdp, newf, x'*Q*x +
            sum((-μ[i])* (x[i]^2 - x[i]) for i in 1:N) -
            β * sum( (algorithm.A_eq[k, :]'*x - algorithm.b_eq[k] )^2 for k in 1:M)
        )
        Q = zeros(N, N) ; c = zeros(N) ; constant = newf.aff.constant

        for (p, v) in newf.terms
            i = sdpIndex[p.a] ; j = sdpIndex[p.b]
            if i == j 
                Q[i,j] = v
            else
                Q[i, j] = v/2 ; Q[j, i] = v/2
            end
        end

        for (k, v) in newf.aff.terms
            i = sdpIndex[k]
            c[i] = v
        end

        return Q, c, constant

    elseif termination_status(model_sdp) == MOI.ITERATION_LIMIT
        @warn("CSDP solver 100 iter limit reached ! ")
    end

    return nothing, nothing, nothing
end
