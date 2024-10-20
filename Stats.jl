#!/usr/bin/env julia

push!(LOAD_PATH, pwd());

# External Imports
using ArgParse
using Statistics
using JLD2

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--n"
            help = "number of runs"
            required = true
        "--g"
            help = "grid id"
            required = true      
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    grid_id = tryparse(Int64, parsed_args["g"])
    n_runs = tryparse(Int64, parsed_args["n"])
    n_cases = 5

    root = "results"
    directory_root = "$root/run"
    filename_root = "workspace-$grid_id"

    times = []

    for run in 1:n_runs
        for case in 1:n_cases
            path = "$directory_root-$run/$filename_root-$case.jld2"
            vars = jldopen(path)

            times = [times; vars["output"].time]
        end
    end

    println("Time avg. over $(n_runs * n_cases) runs: $(mean(times)/60) mins")
end

main()

