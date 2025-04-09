"""
    Write instance of kQKP. 
"""
function invokeHeuristic(λ, model::Optimizer, algorithm::MultiObjectiveBranchBound)
    Q = sum(λ[p] .* algorithm.Qs[p] for p in 1:MOI.output_dimension(model.f))
    fout = open("inst.data", "w")

    # n 
    println(fout, algorithm.nb_vars)

    # [i, j, val]
    mat = []
    for i in 1:algorithm.nb_vars
        for j in 1:algorithm.nb_vars
            val = -round(Q[i, j], digits = 5)
            if val != 0.0
                push!(mat, [i, j, val])
            end
        end
    end

    # println(λ)
    # println(Q)

    println(fout, length(mat))
    for l in mat
        println(fout, round(Int64, l[1]), " ", round(Int64, l[2]), " ", l[3])
    end

    if length(algorithm.b_iq) > 0 

        # i w 
        for i in 1:algorithm.nb_vars
            println(fout, i, " ", round(Int64, algorithm.A_iq[1, i]) ) 
        end

        # W 
        println(fout, round(Int64,algorithm.b_iq[1]) )

        # k-item 
        println(fout, round(Int64,algorithm.b_eq[1]) )

        close(fout)

        # run heuristic program 
        run(`./heurkQKP`)
    else
        c = sum(λ[p] .* algorithm.Ls[p] for p in 1:MOI.output_dimension(model.f))
        for i in 1:algorithm.nb_vars
            println(fout, -round(c[i], digits = 5), " " ) 
        end
        close(fout)

        run(`./heurUBQP`)
    end 

    # read heur sol
    x = []
    fout = open("heur_sol.data", "r")
    n = parse(Int64, split(readline(fout) , " ")[1] )

    for i in 1:n
        l = readline(fout) ; x_ = parse.(Int64, split(l, " ")[1:end-1] )
        # if sum(x_) == algorithm.b_eq[1] && x_'*algorithm.A_iq[1, :] <= algorithm.b_iq[1]
            push!(x, x_)
        # end
    end 

    close(fout)
    # @info "heur x = $x"
    run(`rm -f inst.data`) ; run(`rm -f heur_sol.data`)
    return x
end


function heuristic(model::Optimizer, UBS::Vector{SupportedSolutionPoint}, algorithm::MultiObjectiveBranchBound)
    start_time = time()
    p = MOI.output_dimension(model.f) ; Λ = Set{Vector{Float64}}()
    for i in 1:p
        λ = zeros(p) ; λ[i] = 1.0
        push!(Λ, λ) 
    end
    λ = [1/p for _ in 1:p] ; push!(Λ, λ)
    
    for i in 1:p
        λ_ = (λ .+ collect(Λ)[i]) ./2 ; push!(Λ, λ_) 
    end

    # todo : adjust
    λ_count = ( p * algorithm.nb_vars/5 )

    while λ_count > 0
        i = rand(1:length(Λ)) ; j = rand(1:length(Λ))
        if i != j
            λ_ = (collect(Λ)[j] + collect(Λ)[i]) ./2
            push!(Λ, λ_) ; λ_count -= 1
        end
    end

    # calculate mono 
    # sols = Set()
    for λ in Λ
        # println("heur λ = ", λ)

        X = invokeHeuristic(λ, model, algorithm)
        for x in X 
            push_filtering_dominance(UBS, SupportedSolutionPoint(Set([x]),
                                                        [x'*algorithm.Qs[k]'*x + x'*algorithm.Ls[k] + algorithm.Cs[k] for k in 1:p],
                                                        λ, true
                                        ) 
            )
        end

        # println("UBS = ", UBS)

    end
    
    heuristic_time = round(time() - start_time, digits = 2)
    println("heuristic iters ", length(Λ), " with time(s) ", heuristic_time)
    return heuristic_time
end




function invokeHeuristic_local(λ, model::Optimizer, algorithm::MultiObjectiveBranchBound, varidx, assignment)
    Q = sum(λ[p] .* algorithm.Qs[p] for p in 1:MOI.output_dimension(model.f))


    idx_free = [i for i in 1:algorithm.nb_vars] ; idx_assign =[varidx[x] for x in collect(keys(assignment)) ]
    deleteat!(idx_free, sort(idx_assign))



    fout = open("inst.data", "w")

    # n 
    println(fout, length(idx_free))

    # [i, j, val]
    mat = []
    for i in idx_free
        for j in idx_free
            val = -round(Q[i, j], digits = 5)
            if val != 0.0
                push!(mat, [i, j, val])
            end
        end
    end

    # println(λ)
    # println(Q)

    println(fout, length(mat))
    for l in mat
        println(fout, round(Int64, l[1]), " ", round(Int64, l[2]), " ", l[3])
    end

    if length(algorithm.b_iq) > 0 

        # i w 
        for i in idx_free
            println(fout, i, " ", round(Int64, algorithm.A_iq[1, i]) ) 
        end

        # W 
        W = algorithm.b_iq[1] ; K = algorithm.b_eq[1]
        for (x, val) in assignment
            i = varidx[x]
            W -= val * algorithm.A_iq[1, i]
            K -= val
        end

        println(fout, round(Int64, W) )

        # k-item 
        println(fout, round(Int64, K) )

        close(fout)

        # run heuristic program 
        run(`./heurKQKPprimal`)
    else
        c = sum(λ[p] .* algorithm.Ls[p] for p in 1:MOI.output_dimension(model.f))
        for i in idx_free
            println(fout, -round(c[i], digits = 5), " " ) 
        end
        close(fout)

        run(`./heurUBQP`)
    end 

    # read heur sol
    x = []
    fout = open("heur_sol.data", "r")
    n = parse(Int64, split(readline(fout) , " ")[1] )

    for _ in 1:n
        l = readline(fout) ; x_ = parse.(Int64, split(l, " ")[1:end-1] )
        x_star = zeros(Int64, algorithm.nb_vars) ; i=0
        for j in idx_free
            i +=1
            x_star[j] = x_[i]
        end
        for (var, val) in assignment
            j = varidx[var]
            x_star[j] = round(Int64, val)
        end
        # if sum(x_star) == algorithm.b_eq[1] && x_star'*algorithm.A_iq[1, :] <= algorithm.b_iq[1]
            push!(x, x_star)
        # end

    end 

    close(fout)
    # @info "heur x = $x"
    run(`rm -f inst.data`) ; run(`rm -f heur_sol.data`)
    return x
end





function heuristic_local(model::Optimizer, UBS::Vector{SupportedSolutionPoint}, algorithm::MultiObjectiveBranchBound,
                        assignment :: Dict{MOI.VariableIndex, Float64}, varidx
    )
    start_time = time()
    p = MOI.output_dimension(model.f) ; Λ = Set{Vector{Float64}}()
    for i in 1:p
        λ = zeros(p) ; λ[i] = 1.0
        push!(Λ, λ) 
    end
    λ = [1/p for _ in 1:p] ; push!(Λ, λ)
    
    for i in 1:p
        λ_ = (λ .+ collect(Λ)[i]) ./2 ; push!(Λ, λ_) 
    end

    # todo : adjust 
    λ_count = ( p * (algorithm.nb_vars - length(assignment) )/5 )

    while λ_count > 0
        i = rand(1:length(Λ)) ; j = rand(1:length(Λ))
        if i != j
            λ_ = (collect(Λ)[j] + collect(Λ)[i]) ./2
            push!(Λ, λ_) ; λ_count -= 1
        end
    end

    # calculate mono 
    
    for λ in Λ
        # println("heur λ = ", λ)

        X = invokeHeuristic_local(λ, model, algorithm, varidx, assignment)
        for x in X 
            push_filtering_dominance(UBS, SupportedSolutionPoint(Set([x]),
                                                        [x'*algorithm.Qs[k]'*x + x'*algorithm.Ls[k] + algorithm.Cs[k] for k in 1:p],
                                                        λ, true
                                        ) 
            )
        end

        # println("UBS = ", UBS)

    end
    
    heuristic_time = round(time() - start_time, digits = 4)
    # println("heuristic iters ", length(Λ), " with time(s) ", heuristic_time)
    return heuristic_time
end