using JuMP, LinearAlgebra, CSDP

#todo : improvments : coefficient matrices
function QCR_csdp(Q, c, constant, model::Optimizer, varArray, varIndex)
    N = length(varArray)

    # SDP model ...
    model_sdp = Model(CSDP.Optimizer) ; JuMP.set_silent(model_sdp)
    set_attribute(model_sdp, "axtol", 1.0e-5) ; set_attribute(model_sdp, "atytol", 1.0e-5)

    @variable(model_sdp, x[1:N] )
    sdpIndex = Dict(x[i] => i for i=1:N)

    @variable(model_sdp, X[1:N, 1:N], Symmetric)
    @constraint(model_sdp, [1 x'; x X] in PSDCone())

    # todo : coefficients matrices calculate once only 
    A = [] ; b = Vector{Float64}()          # ctr =
    A_inf = [] ; b_inf = Vector{Float64}()  # ctr <=

    for (t1, t2) in MOI.get(model.inner, MOI.ListOfConstraintTypesPresent())
        ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{t1,t2}())

        # ax = b 
        if t2 == MOI.EqualTo{Float64} && t1 == MOI.ScalarAffineFunction{Float64}
            for ci in ctr_t
                a = zeros(N)
                f_ = MOI.get(model.inner, MOI.ConstraintFunction(), ci) 

                for term in f_.terms
                    i = varIndex[term.variable]
                    a[i] = term.coefficient
                end
                push!(A, a)
                push!(b, MOI.get(model.inner, MOI.ConstraintSet(), ci).value )
            end
            #todo :  VectorAffineFunction
        end

        # ax <= b ScalarAffineFunction{Float64}-in-LessThan{Float64}
        if t2 == MOI.LessThan{Float64} && t1 == MOI.ScalarAffineFunction{Float64}
            for ci in ctr_t
                a = zeros(N)
                f_ = MOI.get(model.inner, MOI.ConstraintFunction(), ci) 

                for term in f_.terms
                    i = varIndex[term.variable]
                    a[i] = term.coefficient
                end
                push!(A_inf, a)
                push!(b_inf, MOI.get(model.inner, MOI.ConstraintSet(), ci).upper )
            end
            #todo :  VectorAffineFunction
        end

        # todo : ax >= b , add other superior constraint 

        # var assignment
        if t1 == MOI.VariableIndex
            if t2 == MOI.GreaterThan{Float64}
                for ci in ctr_t
                    i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
                    @constraint(model_sdp, x[i] >= MOI.get(model.inner, MOI.ConstraintSet(), ci).lower )
                end
            end

            if t2 == MOI.LessThan{Float64}
                for ci in ctr_t
                    i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
                    @constraint(model_sdp, x[i] <= MOI.get(model.inner, MOI.ConstraintSet(), ci).upper )
                end
            end
        end

    end

    A = reduce(vcat,transpose.(A))
    A_inf = reduce(vcat,transpose.(A_inf))
    M = length(b)

    @objective(model_sdp, Min, tr(Q * X) + c'*x + constant )
    @constraint(model_sdp, A*x == b) 
    @constraint(model_sdp, A_inf*x <= b_inf) 

    con_μ = @constraint(model_sdp, [i in 1:N], X[i,i] - x[i] == 0)
    con_β = @constraint(model_sdp, sum( sum(A[k,i]*A[k,j]*X[i,j] for i in 1:N for j in 1:N) - 
                                    2 * sum(b[k] * A[k,j] * x[j] for j in 1:N) +
                                    b[k]^2  for k in 1:M
                                ) == 0 
            ) 

    optimize!(model_sdp)

    if termination_status(model_sdp)== OPTIMAL
        obj = objective_value(model_sdp)

        # println("QCR obj = ", -obj )

        μ = dual.(con_μ)
        β = dual(con_β)

        # println("μ = ", μ)
        # println("β = ", β)

        @expression(model_sdp, newf, x'*Q*x + c'*x + constant +
                                sum((-μ[i])* (x[i]^2 - x[i]) for i in 1:N) -
                                β * sum( (A[k, :]'*x - b[k] )^2 for k in 1:M)
        )
        Q_ = zeros(N, N) ; c_ = zeros(N) ; const_ = newf.aff.constant

        # quad_terms = Array{MOI.ScalarQuadraticTerm{Float64}}(undef, 1)
        
        for (p, v) in newf.terms
            i = sdpIndex[p.a] ; j = sdpIndex[p.b]
            if i == j 
                # push!(quad_terms, MOI.ScalarQuadraticTerm(v*2, varArray[i], varArray[j] ) ) 
                Q_[i,j] = v
            else
                # push!(quad_terms, MOI.ScalarQuadraticTerm(v/2, varArray[i], varArray[j] ) ) 
                # push!(quad_terms, MOI.ScalarQuadraticTerm(v/2, varArray[j], varArray[i] ) ) 
                Q_[i, j] = v/2 ; Q_[j, i] = v/2
            end
        end

        # affine_terms = MOI.ScalarAffineTerm{Float64}[
        #     MOI.ScalarAffineTerm(
        #         v,
        #         varArray[sdpIndex[k]],
        #     ) for (k, v) in newf.aff.terms
        # ]   

        for (k, v) in newf.aff.terms
            i = sdpIndex[k]
            c_[i] = v
        end

        # QCRobj = MOI.ScalarQuadraticFunction(quad_terms, affine_terms, newf.aff.constant)
        # @info "QCRobj : $QCRobj"
        # println(model.inner)
        # MOI.set(model.inner, MOI.ObjectiveFunction{typeof(QCRobj)}(), QCRobj )

        return Q_, c_, const_

    elseif termination_status(model_sdp) == MOI.ITERATION_LIMIT
            @warn("CSDP solver 100 iter limit reached ! ")
    end

    return nothing, nothing, 0.0
end
