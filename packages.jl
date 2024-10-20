using Pkg

dependencies = [
    "CSV",
    "Distributions",
    "StatsPlots",
    "JuMP",
    "Ipopt",
    "PyCall",
    "Combinatorics",
    "IJulia",
    "Plots",
    "DataFrames",
    "JLD2",
    "Setfield",
    "ArgParse",
]

Pkg.add(dependencies)

ENV["GUROBI_HOME"] = "/Library/gurobi950/mac64/"
Pkg.add("Gurobi")
Pkg.build("Gurobi")