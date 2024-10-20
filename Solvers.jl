module Solvers

    using JuMP, Gurobi, Ipopt   
    using Simulations: Grid

    Base.@kwdef struct Profile
        gurobi_env::Gurobi.Env = Gurobi.Env()
        solver::Int32 = 2
        OutputFlag_Gurobi::Int32 = 0 
        silence_solver::Bool = true
        grid::Grid
        dispatch::String #  DA, BM
        par_optimality_gap::Float64 = 2e-2
    end

    function run(p::Profile)
        if p.solver == 1
            return LinearDistFlow(p)
        elseif p.solver == 2
            return DistFlow(p)
        end
    end


    function LinearDistFlow(profile::Profile)

        @assert(profile.dispatch == "DA" || profile.dispatch == "BM")

        ##  Model LINDISTFLOW (linear)
        # m = Model(optimizer_with_attributes(
        #     Gurobi.Optimizer,
        #     "OutputFlag" => profile.OutputFlag_Gurobi,
        # ))

        m = Model(with_optimizer(Gurobi.Optimizer, profile.gurobi_env, OutputFlag=profile.OutputFlag_Gurobi, MIPGap=profile.par_optimality_gap))

        NI = profile.grid.NI
        NL = profile.grid.NL
        Nodes = profile.grid.Nodes
        R = profile.grid.R
        X = profile.grid.X
        Pdd = profile.grid.Pdd
        Qdd = profile.grid.Qdd
        Pmin = profile.grid.Pmin
        Pmax = profile.grid.Pmax
        Qmin = profile.grid.Qmin
        Qmax = profile.grid.Qmax
        Wmin = profile.grid.Wmin
        Wmax = profile.grid.Wmax
        SLmax = profile.grid.SLmax
        c_op = profile.grid.c_op
        c_ls = profile.grid.c_ls
        c_es = profile.grid.c_es
        c_BM⁺ = profile.grid.c_BM⁺
        c_BM⁻ = profile.grid.c_BM⁻
        DA_optimal = profile.grid.DA_optimal

        #define variables
        @variable(m, p[i = 1:NI]) # net active withdraw at node i
        @variable(m, q[i = 1:NI]) # net reactive withdraw at node i
        @variable(m, pl[l = 1:NL]) #active branch power flow from i to j
        @variable(m, ql[l = 1:NL]) #reactive branch power flow from i to j
        @variable(m, w[i = 1:NI]) # square of voltage at node i
        # @variable(m, IL[l = 1:NL]) # square of current from i to j
        if profile.dispatch == "DA"
            @variable(m, pG[i = 1:NI]) #  active geneartion at node i
        elseif profile.dispatch == "BM"
            @variable(m, pG[i=1:NI] == round.(DA_optimal[i],digits=6))
        end
        @variable(m, qG[i = 1:NI]) # reactive generation at node i
        @variable(m, ls[i=1:NI] >= 0) # load sheding at each node.
        @variable(m, es[i=1:NI] >= 0) # energy spillage at each node.

        if profile.dispatch == "BM"
            @variable(m, ΔpG⁺[i=1:NI] >= 0, start = 0) #  postive BALANCING active geneartion at node i
            @variable(m, ΔpG⁻[i=1:NI] >= 0, start = 0) #  negative BALANCING active geneartion at node i
        end

        @constraint( m, DistFlowEqP[l = 1:NL],
            pl[l] ==  p[Nodes[l, 2]] +  sum(pl[k] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
        @constraint(m, DistFlowEqQ[l = 1:NL],
            ql[l] == q[Nodes[l, 2]]   + sum(ql[k] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
        @constraint(m, VoltageLosses[l = 1:NL],
            w[Nodes[l, 2]] == w[Nodes[l, 1]] - 2 * (R[l] * pl[l] + X[l] * ql[l]))

        # nodal net withdraw
        if profile.dispatch == "DA"
            @constraint(m, NodalEqP[i = 1:NI], p[i] == -pG[i] + Pdd[i] - ls[i] + es[i]);
        elseif profile.dispatch == "BM"
            @constraint(m, NodalEqP[i=1:NI], p[i] == -(pG[i] + ΔpG⁺[i] - ΔpG⁻[i]) + Pdd[i] - ls[i] + es[i]);
        end
        @constraint(m, NodalEqQ[i = 1:NI], q[i] == -qG[i] + Qdd[i]);

        # limits on load sheding and energy spillager
        @constraint(m, [i=1:NI], ls[i] <= Pdd[i])
        @constraint(m, [i=1:NI], es[i] <= pG[i])

        #Technical generation limits
        if profile.dispatch == "DA"
            @constraint(m, GenLimitsP[i = 1:NI], Pmin[i] <= pG[i] <= Pmax[i])
        elseif profile.dispatch == "BM"
            @constraint(m, GenLimitsP[i=1:NI], Pmin[i] <= pG[i] + ΔpG⁺[i] - ΔpG⁻[i]<= Pmax[i])
        end
        @constraint(m, GenLimitsQ[i = 1:NI], Qmin[i] <= qG[i] <= Qmax[i])
        #Voltage limits
        @constraint(m, VoltageLimits[i = 1:NI], Wmin <= w[i] <= Wmax)
        #Line capacity limits
        @constraint(m, FlowLimitsMax[l = 1:NL], pl[l] <= SLmax[l])
        @constraint(m, FlowLimitsMin[l = 1:NL], pl[l] >= - SLmax[l])

        #Special case for 1st node
        if profile.dispatch == "DA"
            @constraint(m, LimitsP, pG[1]==pl[1])
        elseif profile.dispatch == "BM"
            @constraint(m, LimitsP, (pG[1] + ΔpG⁺[1] - ΔpG⁻[1]) == pl[1])
        end
        @constraint(m, LimitsQ, qG[1]==ql[1])
        @constraint(m, SELimitsW, w[1] == 1)

        # Minimization of the total cost
        if profile.dispatch == "DA"
            @expression(m, disp_cost, sum(c_op[i] * pG[i] for i = 1:NI) )
            @expression(m, emer_cost, sum(c_ls*ls[i] + c_es*es[i] for i = 1:NI) )
            @objective(m, Min, disp_cost + emer_cost )   
        elseif profile.dispatch == "BM"
            @objective(m, Min, sum(c_BM⁺[i]*ΔpG⁺[i] + c_BM⁻[i]*ΔpG⁻[i]  + c_ls*ls[i] + c_es*es[i] for i=1:NI))
        end
        
        set_optimizer_attribute(m, MOI.Silent(), profile.silence_solver)

        optimize!(m)

        if profile.dispatch == "DA"
            return value.(pG), value(disp_cost)
        elseif profile.dispatch == "BM"
            return objective_value(m), value.(ΔpG⁺), value.(ΔpG⁻)
        end
    end

    function DistFlow(profile::Profile)
        ##  Model DISTFLOW (non-linear)
        # m = Model(optimizer_with_attributes(
        #     Gurobi.Optimizer,
        #     "NonConvex" => 2,
        #     "OutputFlag" => profile.OutputFlag_Gurobi,
        # ))

        m = Model(with_optimizer(Gurobi.Optimizer, profile.gurobi_env, NonConvex = 2, OutputFlag=profile.OutputFlag_Gurobi, MIPGap=profile.par_optimality_gap))

        NI = profile.grid.NI
        NL = profile.grid.NL
        Nodes = profile.grid.Nodes
        R = profile.grid.R
        X = profile.grid.X
        Pdd = profile.grid.Pdd
        Qdd = profile.grid.Qdd
        Pmin = profile.grid.Pmin
        Pmax = profile.grid.Pmax
        Qmin = profile.grid.Qmin
        Qmax = profile.grid.Qmax
        Wmin = profile.grid.Wmin
        Wmax = profile.grid.Wmax
        SLmax = profile.grid.SLmax
        c_op = profile.grid.c_op
        c_ls = profile.grid.c_ls
        c_es = profile.grid.c_es
        c_BM⁺ = profile.grid.c_BM⁺
        c_BM⁻ = profile.grid.c_BM⁻
        DA_optimal = profile.grid.DA_optimal

        #define variables
        @variable(m, p[i = 1:NI]) # net active withdraw at node i
        @variable(m, q[i = 1:NI]) # net reactive withdraw at node i
        @variable(m, pl[l = 1:NL]) #active branch power flow from i to j
        @variable(m, ql[l = 1:NL]) #reactive branch power flow from i to j
        @variable(m, w[i = 1:NI]) # square of voltage at node i
        @variable(m, IL[l = 1:NL]) # square of current from i to j
        if profile.dispatch == "DA"
            @variable(m, pG[i = 1:NI]) #  active geneartion at node i
        elseif profile.dispatch == "BM"
            @variable(m, pG[i=1:NI] == round.(DA_optimal[i],digits=6))
        end
        @variable(m, qG[i = 1:NI]) # reactive generation at node i
        @variable(m, ls[i=1:NI] >= 0) # load sheding at each node.
        @variable(m, es[i=1:NI] >= 0) # energy spillage at each node.
        if profile.dispatch == "BM"
            @variable(m, ΔpG⁺[i=1:NI] >= 0, start = 0) #  postive BALANCING active geneartion at node i
            @variable(m, ΔpG⁻[i=1:NI] >= 0, start = 0) #  negative BALANCING active geneartion at node i
        end

        @constraint( m, DistFlowEqP[l = 1:NL],
            pl[l] ==  p[Nodes[l, 2]] + R[l] * IL[l] + sum(pl[k] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
        @constraint(m, DistFlowEqQ[l = 1:NL],
            ql[l] == q[Nodes[l, 2]] + X[l] * IL[l] + sum(ql[k] for k = 1:NL if Nodes[k, 1] == Nodes[l, 2]))
        @constraint(m, VoltageLosses[l = 1:NL],
            w[Nodes[l, 2]] == w[Nodes[l, 1]] + (R[l]^2 + X[l]^2) * IL[l] - 2 * (R[l] * pl[l] + X[l] * ql[l]))
        @constraint(m,PowerFlow_as_V_times_I[l = 1:NL], (pl[l])^2 + (ql[l])^2 == IL[l] * w[Nodes[l, 1]] )

        # nodal net withdraw
        if profile.dispatch == "DA"
            @constraint(m, NodalEqP[i = 1:NI], p[i] == -pG[i] + Pdd[i] - ls[i] + es[i]);
        elseif profile.dispatch == "BM"
            @constraint(m, NodalEqP[i=1:NI], p[i] == -(pG[i] + ΔpG⁺[i] - ΔpG⁻[i]) + Pdd[i] - ls[i] + es[i]);
        end
        @constraint(m, NodalEqQ[i = 1:NI], q[i] == -qG[i] + Qdd[i]);

        # limits on load sheding and energy spillager
        @constraint(m, [i=1:NI], ls[i] <= Pdd[i])
        @constraint(m, [i=1:NI], es[i] <= pG[i])

        #Technical generation limits
        if profile.dispatch == "DA"
            @constraint(m, GenLimitsP[i = 1:NI], Pmin[i] <= pG[i] <= Pmax[i])
        elseif profile.dispatch == "BM"
            @constraint(m, GenLimitsP[i=1:NI], Pmin[i] <= pG[i] + ΔpG⁺[i] - ΔpG⁻[i]<= Pmax[i])
        end
        @constraint(m, GenLimitsQ[i = 1:NI], Qmin[i] <= qG[i] <= Qmax[i])
        #Voltage limits
        @constraint(m, VoltageLimits[i = 1:NI], Wmin <= w[i] <= Wmax)
        #Line capacity limits
        @constraint(m, FlowLimits[l = 1:NL], (pl[l])^2 + (ql[l])^2 <= (SLmax[l])^2)

        #Special case for 1st node
        if profile.dispatch == "DA"
            @constraint(m, LimitsP, pG[1]==pl[1])
        elseif profile.dispatch == "BM"
            @constraint(m, LimitsP, (pG[1] + ΔpG⁺[1] - ΔpG⁻[1]) == pl[1])
        end
        @constraint(m, LimitsQ, qG[1]==ql[1])
        @constraint(m, SELimitsW, w[1] == 1)

        # Minimization of the total cost
        if profile.dispatch == "DA"
            @expression(m, disp_cost, sum(c_op[i] * pG[i] for i = 1:NI) )
            @expression(m, emer_cost, sum(c_ls*ls[i] + c_es*es[i] for i = 1:NI) )
            @objective(m, Min, disp_cost + emer_cost )
        elseif profile.dispatch == "BM"
            @objective(m, Min, sum(c_BM⁺[i]*ΔpG⁺[i] + c_BM⁻[i]*ΔpG⁻[i]  + c_ls*ls[i] + c_es*es[i] for i=1:NI))
        end

        set_optimizer_attribute(m, MOI.Silent(), profile.silence_solver)

        optimize!(m)

        print(value.(w))

        if profile.dispatch == "DA"
            return value.(pG), value(disp_cost)
        elseif profile.dispatch == "BM"
            return objective_value(m), value.(ΔpG⁺), value.(ΔpG⁻)
        end
    end
end