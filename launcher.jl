#!/usr/bin/env julia

push!(LOAD_PATH, pwd());

# Internal Imports
using Simulations
using Solvers
using UtilityFunction
using Scenarios

# External Imports
using Setfield
using ArgParse
using JLD2

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--g"
            help = "grid id"
            required = true
        "--c"
            help = "case id"
            required = true
        "--a"
            help = "check dependence on alpha"
            default = false
            required = false        
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    grid_id = tryparse(Int64, parsed_args["g"])
    case_id = tryparse(Int64, parsed_args["c"])
    grid = Scenarios.get_scenario_by_id(grid_id, case_id)
    profile = Simulations.init(grid)
    grid_data = Simulations.load_data_from_file(profile)

    if parsed_args["a"] == "true"
        for alpha in 0.1:0.1:1

            print(grid.demand)

            demand = (wo_EV = grid.demand.wo_EV,
                true_EV = alpha .* grid.demand.true_EV .* grid.alpha^(-1),
                best_EV = alpha .* grid.demand.best_EV .* grid.alpha^(-1)
            )

            grid = @set grid.demand = demand
            grid = @set grid.alpha = alpha

            print(grid.demand)

            profile = Simulations.init(grid)
            grid_data = Simulations.load_data_from_file(profile)
            output = @timed UtilityFunction.run(profile.grid_settings, grid_data)
            v = output.value
            println("\n\nAlpha:$alpha")
            println("\n\nTime:$t")
            println("\n\nResults:$v")
        end
    else
        output = @timed UtilityFunction.run(profile.grid_settings, grid_data)
        v = output.value
        t = output.time
        println("\n\nTime:$t")
        println("\n\nResults:$v")
        
    end

    jldsave("results/workspace-$grid_id-$case_id.jld2"; profile, grid_data, output)
end

main()
