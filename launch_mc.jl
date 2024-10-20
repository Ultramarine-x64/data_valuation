#!/usr/bin/env julia

push!(LOAD_PATH, pwd());

# Internal Imports
using Simulations
using Solvers
using UtilityFunction
using Allocations
using Scenarios

# External Imports
using Setfield
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--i"
            help = "input file"
            required = true
        "--g"
            help = "grid id"
            required = true
        "--c"
            help = "case id"
            required = true
        "--m"
            help = "number of monte carlo simulations"
            default = 50
            arg_type = Int
            required = false
        "--n"
            help = "number of monte carlo runs for averaging"
            default = 10
            arg_type = Int
            required = false
    end

    return parse_args(s)
end

function read_file(filename)
    result = NaN
    open(filename) do file
        output = string(read(file, String))
        svec = split(strip(output, ['[', ']', ' ', '\n']), ',')
        result = map(x->tryparse(Float64, x), svec)
    end

    return result
end

function run_SV(grid, utility, exp)
    res = Allocations.run_monte_carlo(grid, utility, exp)
    res = round.(res; digits=4)
    print("MC_res:$res")
end

function main()
    parsed_args = parse_commandline()
    utility = read_file(parsed_args["i"])

    exp = Simulations.Experiment(
        m=parsed_args["m"],
        n=parsed_args["n"]
    )
    grid = Scenarios.get_scenario_by_id(tryparse(Float64, parsed_args["g"]), tryparse(Float64, parsed_args["c"]))
    run_SV(grid, utility, exp)
end

main()