
using JuMP, CSDP
# include("MOBB.jl")
# include("MOBBTree.jl")

"""
For each level k ∈[0,n-1], convexify the ∑₍ᵢ₌ₖ₊₁₎ ∑₍ⱼ₌ₖ₊₁₎ xᵢ*qᵢⱼ*xⱼ term. 

`N` = n-k , the number of free variables 

return μ or nothing in case of failure. 
"""
function klevel_UQCR_csdp(N::Int64, Q::Matrix{Float64}, k::Int64, algorithm::MultiObjectiveBranchBound )

    # SDP model ...
    model_sdp = Model(CSDP.Optimizer) ; JuMP.set_silent(model_sdp) 
    set_attribute(model_sdp, "axtol", 1.0e-5) ; set_attribute(model_sdp, "atytol", 1.0e-5)
    set_attribute(model_sdp, "maxiter", 200) 

    @variable(model_sdp, x[1:N] )
    @variable(model_sdp, X[1:N, 1:N], Symmetric)
    @constraint(model_sdp, [1 x'; x X] in PSDCone())

    @objective(model_sdp, Min, tr(Q * X) )
    
    con_μ = @constraint(model_sdp, [i in 1:N], X[i,i] - x[i] == 0)

    # --------------------------
    # todo : relaxed ctr preproc=2
    if MOI.get(algorithm, Preproc()) == 2
        length(algorithm.b_iq) > 0 ? @constraint(model_sdp, algorithm.A_iq[:, k+1:end]* x ≤ algorithm.b_iq ) : nothing
        length(algorithm.b_eq) > 0 ? @constraint(model_sdp, algorithm.A_eq[:, k+1:end]* x ≤ algorithm.b_eq ) : nothing
    end

    optimize!(model_sdp)

    if termination_status(model_sdp)== OPTIMAL
        μ = dual.(con_μ)
    
        return μ
    elseif termination_status(model_sdp) == MOI.ITERATION_LIMIT
        @error("CSDP solver 200 iter limit reached for level $k ! ")
        return nothing
    else
        @error("CSDP solver failed for level $k ! ")
        return nothing
    end

    return nothing
end


function klevel_solve_weighted_sum( k :: Int64,
            node :: Node,
            model::Optimizer,
            λ::Vector{Float64},
            algorithm::MultiObjectiveBranchBound
    )
    varArray = MOI.get(model, MOI.ListOfVariableIndices())
    N = length(varArray)
    varArray_inner = MOI.get(model.inner, MOI.ListOfVariableIndices())

    varIndex_inner = Dict(varArray[i] => i for i=1:N)
    varIndex = Dict(i => varArray[i] for i=1:N)

    λQ = sum( λ[p].* algorithm.Qs[p] for p in 1:length(λ))

    Q = zeros(N, N) ; c = sum( λ[p].* algorithm.Ls[p] for p in 1:length(λ)) ; constant = sum( λ[p] * algorithm.Cs[p] for p in 1:length(λ)) 

    # i=1...k ; j=1...k
    for i in 1:k
        for j in 1:k
            constant += node.assignment[varIndex[i]] * λQ[i,j] * node.assignment[varIndex[j]] 
        end
    end

    # i=1...k, j=k+1...n 
    for i in 1:k
        for j in k+1:N
            c[j] += node.assignment[varIndex[i]] * λQ[i,j]
        end
    end

    # i=k+1...n, j=1...k 
    for i in k+1:N
        for j in 1:k
            c[i] += λQ[i,j] * node.assignment[varIndex[j]] 
        end
    end

    # i=k+1...n , j=k+1...n
    for i in k+1:N
        for j in k+1:N
            Q[i,j] += λQ[i,j]
        end
    end

    # ------------------------------
    # -- QCR coeff ∑₍ⱼ₌ₖ₊₁₎ μⱼ(xⱼ^2 - xⱼ)
    if k+1 <= algorithm.nb_vars
        μ = sum(λ[p] .* algorithm.preproc_μ[p][k+1] for p in 1:length(λ))
        for j in k+1:N
            Q[j, j] += -μ[j-k]
            c[j] += μ[j-k]
        end
    end

    quad_terms = MOI.ScalarQuadraticTerm{Float64}[
        MOI.ScalarQuadraticTerm(
            i==j ? Q[i, j]*2 + (i≥k+1 ? 0.0001 : 0.0) : Q[i, j],
            varArray_inner[i],
            varArray_inner[j],
        ) for i in 1:N for j in 1:N
    ]

    affine_terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(
            c[i],
            varArray_inner[i],
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
    for x in varArray_inner
        sol[varIndex_inner[x]] = MOI.get(model.inner, MOI.VariablePrimal(), x)
    end

    Y = [ sol'*algorithm.Qs[p]*sol + sol'*algorithm.Ls[p] + algorithm.Cs[p] + (k+1 <= algorithm.nb_vars ? sum(-μ[j-k]*(sol[j]^2 - sol[j]) for j in k+1:N) : 0 )
             for p in 1:length(λ)
        ]

    return status, sol, Y

end 