module Scenarios

    using Simulations

    function choice(a::Array)
        n = length(a)
        idx = rand(1:n)
        return a[idx]
    end

    function get_scenario_by_id(grid_id, case_id)
        
        # Common parameters
        c_ls = 100
        c_es = 20
        factor_line_capacity_scale = 0.3

        grid_id = Int(grid_id)
        filename = "nodes/indices_$grid_id.txt"
        ind_EV_loc = NaN

        open(filename) do file
            output = string(read(file, String))
            svec = split(strip(output, ['[', ']', ' ', '\n']), ',')
            ind_EV_loc = sort(map(x->tryparse(Int32, x), svec))
        end

        grid_size = grid_id
        alpha = 1/length(ind_EV_loc)

        # cost for balancing
        c_BM⁺ = zeros(grid_size,1)
        c_BM⁺[1] = 2.5  # substation
        c_BM⁺[2:end] .= 5.5  # distributed gencos

        c_BM⁻ = zeros(grid_size,1)
        c_BM⁻[1] = -0.5
        c_BM⁻[2:end] .= -0.5;

        demand_wo_EV = vec(zeros(grid_size,1))
        demand_wo_EV[2:end] .= 0.1

        true_EV_demand = zeros(grid_size, 1)
        true_EV_demand[ind_EV_loc] .= alpha * 1.5

        # rand_component = [1.16, 1.40, 1.56, 1.63, 1.22, 1.44, 1.53, 1.66, 1.86, 1.29]
        # true_EV_demand[ind_EV_loc] .= alpha * rand_component

        if grid_id == 5

            # best forecast demand
            best_EV_demand_forecast = zeros(grid_size, 1)
            if case_id == 1
                # original case: underestimation
                best_EV_demand_forecast[ind_EV_loc] .= alpha * 1.

            elseif case_id == 2
                # second case: overestimation
                best_EV_demand_forecast[ind_EV_loc] .= alpha * 2.

            elseif case_id == 3
                # third case: uneven estimation
                x = vec([1. 1.5 2.])
                best_EV_demand_forecast[ind_EV_loc] .= alpha * [choice(x) for i in 1:size(ind_EV_loc)[1]]

            elseif case_id == 4
                # four case: random
                rand_component = rand(size(ind_EV_loc)[1]) .+ 1
                #rand_component = [-0.0881, -0.3833, 0.399, -0.1]
                best_EV_demand_forecast[ind_EV_loc] .= alpha * rand_component

            elseif case_id == 5
                # original case: uneven understimation
                best_EV_demand_forecast[ind_EV_loc] .= alpha * range(1., 1.5, length = length(ind_EV_loc))
            end
        elseif grid_id == 33
            # best forecast demand
            best_EV_demand_forecast = zeros(grid_size, 1)
            if case_id == 1
                # original case: underestimation
                best_EV_demand_forecast[ind_EV_loc] .= alpha * 1.

            elseif case_id == 2
                # second case: overestimation
                best_EV_demand_forecast[ind_EV_loc] .= alpha * 2.

            elseif case_id == 3
                # third case: uneven estimation
                x = vec([1. 1.5 2.])
                best_EV_demand_forecast[ind_EV_loc] .= alpha * [choice(x) for i in 1:size(ind_EV_loc)[1]]

            elseif case_id == 4
                # four case: random
                rand_component = rand(size(ind_EV_loc)[1]) .+ 1
                #rand_component = [1.9311151512445586, 1.4389389593310216, 1.2468624804749107, 1.011819583479107, 1.0460428263964987]
                best_EV_demand_forecast[ind_EV_loc] .= alpha * rand_component

            elseif case_id == 5
                # original case: uneven underestimation
                best_EV_demand_forecast[ind_EV_loc] .= alpha * range(1.,1.5,length = length(ind_EV_loc))
            end 
        end

        demand = (wo_EV = demand_wo_EV,
            true_EV = true_EV_demand,
            best_EV = best_EV_demand_forecast
        );

        grid = Simulations.GridSettings(
            grid_size = grid_size,
            case_id = case_id,
            alpha = alpha,
            factor_line_capacity_scale = factor_line_capacity_scale,
            ind_EV_loc = ind_EV_loc,
            demand=demand,
            c_ls = c_ls,
            c_es = c_es,
            c_BM⁺ = c_BM⁺,
            c_BM⁻ = c_BM⁻
        );

    return grid
    end
end