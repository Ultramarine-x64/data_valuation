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
    end

    return parse_args(s)
end

using ArgParse

function read_file(filename)
    result = NaN
    open(filename) do file
        output = string(read(file, String))
        svec = split(strip(output, ['[', ']', ' ', '\n']), ',')
        result = map(x->tryparse(Float64, x), svec)
    end

    return result
end

function run_SV(grid, utility)
    profile = Simulations.init(grid)
    coalitions = Simulations.generate_coalitions(profile.grid_settings);

    #coalitions = coalitions[2:end, 2:end];
    Mat = Allocations.RunAllocations(coalitions, utility, profile.grid_settings.ind_EV_loc)

    SV_res = Mat[:,1]
    CIS_EANS = Mat[:,end]
    SV_res = round.(SV_res; digits=4)
    CIS_EANS = round.(CIS_EANS; digits=4)
    print("\n\nSV_res:$SV_res")
    print("\n\nCIS_EANS_res:$CIS_EANS")
end

function main()
    @show parsed_args = parse_commandline()
    utility = read_file(parsed_args["i"])
    grid = Scenarios.get_scenario_by_id(tryparse(Float64, parsed_args["g"]), tryparse(Float64, parsed_args["c"]))
    run_SV(grid, utility)
end

main()