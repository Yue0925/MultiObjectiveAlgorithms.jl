"""
    Write instance of kQKP. 
"""
function invokeSOheuristic(λ, model::Optimizer, algorithm::MultiObjectiveBranchBound)
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
    run(`./BestHybrid_objreel_yue`)

    # read heur sol
    x = zeros(Int64, algorithm.nb_vars)
    fout = open("heur_sol.data", "r")
    readline(fout) 
    lines = readlines(fout)
    for l in lines
        v = split(l, " ")
        i = parse(Int64, v[1])
        x[i] = parse(Int64, v[2])
    end

    close(fout)
    # @info "heur x = $x"
    run(`rm -f inst.data`) ; run(`rm -f heur_sol.data`)
    return x
end


function heuristic(model::Optimizer, UBS::Vector{SupportedSolutionPoint}, algorithm::MultiObjectiveBranchBound)
    start_time = time()
    p = MOI.output_dimension(model.f) ; Λ = []
    for i in 1:p
        λ = zeros(p) ; λ[i] = 1.0
        push!(Λ, λ) 
    end
    λ = [1/p for _ in 1:p] ; push!(Λ, λ)
    
    for i in 1:p
        λ_ = (λ .+ Λ[i]) ./2 ; push!(Λ, λ_) 
    end

    λ_count = 2*( p * algorithm.nb_vars/5 )

    while λ_count > 0
        i = rand(1:length(Λ)) ; j = rand(1:length(Λ))
        if i != j
            λ_ = (Λ[j] + Λ[i]) ./2
            push!(Λ, λ_) ; λ_count -= 1
        end
    end

    # calculate mono 
    # sols = Set()
    for λ in Λ
        x = invokeSOheuristic(λ, model, algorithm)

        push_filtering_dominance(UBS, SupportedSolutionPoint(Set([x]),
                                                             [x'*algorithm.Qs[k]'*x for k in 1:p],
                                                             λ, true
                                    ) )
    end
    
    heuristic_time = round(time() - start_time, digits = 4)
    println("heuristic iters ", length(Λ), " with time(s) ", heuristic_time)
    return heuristic_time
end