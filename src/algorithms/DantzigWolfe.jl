using JuMP, Gurobi

include("MOBB.jl")



function column_generation_algorithm(model, algorithm::MultiObjectiveBranchBound, Q, c, constant, Q_)
    # variable assignment
    N = algorithm.nb_vars
    varArray = MOI.get(model.inner, MOI.ListOfVariableIndices())
    varIndex = Dict(varArray[i] => i for i=1:N)
    assign = Dict{Int64, Int64}()

    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.GreaterThan{Float64}}())
    for ci in ctr_t
        if MOI.get(model.inner, MOI.ConstraintSet(), ci).lower == 1.0
            i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
            assign[i] = 1
        end
    end

    ctr_t =  MOI.get(model.inner, MOI.ListOfConstraintIndices{MOI.VariableIndex, MOI.LessThan{Float64}}())
    for ci in ctr_t
        if MOI.get(model.inner, MOI.ConstraintSet(), ci).upper == 0.0
            i = varIndex[ MOI.get(model.inner, MOI.ConstraintFunction(), ci) ]
            assign[i] = 0
        end
    end


    # init set knapsack polytope {x | ax <= b , x ∈{0,1}^n }
    P_kn = [ ]
    println("column_generation_algorithm  ... " )
    println("var assignment ", assign)

    while length(P_kn) == 0
        
        for _ in 1:algorithm.nb_vars
            x = zeros(Int64, algorithm.nb_vars)
            objects = [i for i in 1:algorithm.nb_vars]

            for (i, v) in assign
                if v ==1 x[i] = 1 end 
            end
            if length(assign)>0 deleteat!(objects, sort(collect(keys(assign))) ) end 
            
            if sum(x) == algorithm.b_eq[1]
                if x'* algorithm.A_iq[1, :] ≤ algorithm.b_iq[1] 
                    return OPTIMAL, x .* 1.0, [x'*algorithm.Qs[p]*x for p in 1:length(algorithm.Qs) ]
                else
                    return INFEASIBLE, nothing, nothing
                end
            end

            for _ in 1:(algorithm.b_eq[1] - sum(x))
                j = rand(1:length(objects))
                i = objects[j] ; deleteat!(objects, j)
                if algorithm.b_iq[1] - x'* algorithm.A_iq[1, :] - algorithm.A_iq[1, i] >= 0
                    x[i] = 1
                end
            end
            
            sum(x) == algorithm.b_eq[1] ? push!(P_kn, x) : nothing
        end

    end

    # @assert length(P_kn) > 0
    println("init set ", P_kn)

    iter = 0

    while iter <100 
        println("iter ", iter )
        iter += 1
        # solve master problem linear 
        master_model = Model(Gurobi.Optimizer) ; set_silent(master_model)
        n = length(P_kn)
        @variable(master_model, y[1:n] ≥ 0)
        
        @objective(master_model, Max, -sum( y[p] * (P_kn[p]'* Q_ *P_kn[p]) for p in 1:n) )

        con_π = @constraint(master_model, sum(y) == 1)
        con_μ = @constraint(master_model, sum( sum(P_kn[p] .* y[p]) for p in 1:n) == algorithm.b_eq[1])

        optimize!(master_model)

        if termination_status(master_model) != OPTIMAL
            @error("master problem is not solved ! (iter $iter) ")
            break 
        end

        master_obj = objective_value(master_model) ;  y_sol = value.(y)
        μ = dual(con_μ) ; π = dual(con_π)

        println("master obj = ", master_obj)
        println("y_sol = ", y_sol)
        println("μ  = $μ  π = $π ")

        # solve pricing problem quadratic
        pricing_model = Model(Gurobi.Optimizer) ; set_silent(pricing_model)
        @variable(pricing_model, x[1:algorithm.nb_vars], Bin)
        @constraint(pricing_model, algorithm.A_iq[1, :]'*x ≤ algorithm.b_iq[1])
        @objective(pricing_model, Max, π + sum(μ * x[i] for i in 1:algorithm.nb_vars) -
                                            (x'*Q*x + x'*c + constant ) 
                )

        for (i, v) in assign
            @constraint(pricing_model, x[i] == v)
        end

        optimize!(pricing_model)

        # todo : convex qcr objective de pricer
        

        if termination_status(pricing_model) != OPTIMAL
            @error("pricing problem is not solved ! (iter $iter) ")
            break 
        end

        pricing_obj = objective_value(pricing_model) ; x_sol = value.(x)

        println("pricing obj = ", pricing_obj)
        println("x = ", x_sol)

        # todo :  tolerance
        if pricing_obj ≥ 1e-4
            push!(P_kn, x_sol)

        else
            sol = [sum(P_kn[p][i] * y_sol[p] for p in 1:n) for i in 1:algorithm.nb_vars]
            println("sol = ", sol )
            return OPTIMAL, sol, [sol'*algorithm.Qs[p]*sol for p in 1:length(algorithm.Qs) ]

        end

    end

    return INFEASIBLE, nothing, nothing
end