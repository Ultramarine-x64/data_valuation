module UtilityFunction

    using Simulations
    using Solvers
    using DataFrames
    using Setfield

    function run(profile::Simulations.GridSettings, grid::Simulations.Grid)

        ## dataframe for solutions
        df = DataFrame(
            CaseStudy = String[],
            cost_DA = Float64[],
            cost_BM = Float64[],
            p_DA = Array{Float64,1}[],
            Δp_up = Array{Float64,1}[],
            Δp_down = Array{Float64,1}[],
        )

        ## No EV-aware
        case = "No EV-aware"
        DA_demand = vec(profile.demand.wo_EV)
        BM_demand = vec(profile.demand.wo_EV + profile.demand.true_EV)
        df = run_and_save(df, case, DA_demand, BM_demand, profile, grid)

        # Save the used demand
        grid_size = profile.grid_size
        case_id = profile.case_id
        open("results/demand-$grid_size-$case_id.txt", "w") do file
            print(file, profile.demand)
        end

        # Run for coalitions
        coalitions = Simulations.generate_coalitions(profile)
        df = run_for_coalitions(df, coalitions, profile, grid)
        v = get_utility(df)

        return v
    end

    function run_and_save(df::DataFrame, case, DA_demand::Vector{Float64}, BM_demand::Vector{Float64}, profile::Simulations.GridSettings, grid::Simulations.Grid)
        grid.NodesData.Pd = DA_demand
        grid = Simulations.update_data(grid)

        DA_profile = Solvers.Profile(solver=profile.solver, grid=grid, dispatch="DA")
        p_optimal_DA, total_cost_DA = Solvers.run(DA_profile)

        grid.NodesData.Pd = BM_demand
        grid = Simulations.update_data(grid)
        grid = @set grid.DA_optimal = p_optimal_DA;

        BM_profile = Solvers.Profile(solver=profile.solver, grid=grid, dispatch="BM")
        total_cost_BL, Δp⁺_optimal_BL_val, Δp⁻_optimal_BL_val = Solvers.run(BM_profile)

        push!(
            df,
            (
                case,
                total_cost_DA,
                total_cost_BL,
                round.(p_optimal_DA[[1, 4]], digits = 3),
                round.(Δp⁺_optimal_BL_val[vcat([1], profile.ind_EV_loc)], digits = 3),
                round.(Δp⁻_optimal_BL_val[vcat([1], profile.ind_EV_loc)], digits = 3),
            ),
        )  
        
        return df
    end

    function run_for_coalitions(df::DataFrame, coalitions::Matrix{Int64}, profile::Simulations.GridSettings, grid::Simulations.Grid)
        for c in 1:size(coalitions, 1)
            case = string("Coalition: ", coalitions[c, :])
            i_true_data = findall(coalitions[c, :] .== 1)
            i_no_data = findall(coalitions[c, :] .== 0)

            demand = zeros(size(grid.NodesData.Pd, 1), 1)
            demand[i_true_data] = vec(profile.demand.wo_EV[i_true_data] + profile.demand.true_EV[i_true_data])
            demand[i_no_data] = vec(profile.demand.wo_EV[i_no_data] + profile.demand.best_EV[i_no_data])
            # demand[i_true_data] = vec(profile.demand.wo_EV[i_true_data] + profile.demand.best_EV[i_true_data])


            DA_demand = vec(demand)
            BM_demand = vec(profile.demand.wo_EV + profile.demand.true_EV)
            df = run_and_save(df, case, DA_demand, BM_demand, profile, grid)   
        end

        return df
    end

    function get_utility(df::DataFrame)
        df[!, 2:3] = round.(df[:, 2:3], digits = 3)
        df[!, :TotalCost] = df[!, :cost_DA] + df[!, :cost_BM]
        df[!, :Δcost_no_data] = df[!, :TotalCost] .- df[end, :TotalCost]
        df[!, 7:8] = round.(df[:, 7:8], digits = 4)

        ini_index = 3
        vf = df[ini_index:end, :TotalCost]

        cNd = df[ini_index-1, :TotalCost]  # cost of no data
        cNd = cNd[1]
        v = cNd .- vf

        return v
    end

end