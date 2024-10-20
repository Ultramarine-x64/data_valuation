module Allocations

    using JuMP, Gurobi, Ipopt   
    using Combinatorics
    using Statistics
    using Simulations
    using Random


    function RunAllocations(coalitions,v, ind_EV_loc)
        ref_no_value = 0.
        coalitions = coalitions[2:end,ind_EV_loc]
        NC, NP = size(coalitions)

        if !(@isdefined printlevel_Ipopt)
            printlevel_Ipopt = 6
        end

        # indexes
        aux = sum(coalitions, dims=2)
        ind_grand_coal = findall(aux[:].== NP)
        ind_individual = findall(aux[:].== 1)

        # grand coalition value
        vN = v[ind_grand_coal[1]]
    ## Shapley Value
        SV = zeros(NP,1)
        for i=1:NP
            ind_S = findall(coalitions[:,i] .== 1)
            for s = 1:length(ind_S)
                c =  coalitions[ind_S[s],:]
                NS = sum(c)
                c[i] = 0
                ind_S_no_i = findall([sum([coalitions[k,j] .== c[j] for j=1:NP]) == NP for k=1:NC] .== true)

                val_S = v[ind_S[s]]
                if isempty(ind_S_no_i)
                    val_S_no_i = ref_no_value
                else
                    val_S_no_i = v[ind_S_no_i]
                end
                marginal_i_in_S = val_S[1] - val_S_no_i[1]
                weight_S =  (factorial(big(NP - NS))*factorial(big(NS-1)))/factorial(big(NP))

                SV[i] = SV[i] + weight_S*marginal_i_in_S
            end
        end


    ## CIS center of gravity of imputation set value

        CIS = zeros(NP,1)

        # remainder worth
        RW = vN - sum(v[ind_individual])

        for i=1:NP
            c = Int8[0 for i=1:NP]
            c[i] = 1
            ind_S_no_i = findall([sum([coalitions[k,j] .== c[j] for j=1:NP]) == NP for k=1:NC] .== true)

            CIS[i] = v[ind_S_no_i[1]] + (1/NP)*RW
        end


    ## EANS value, equal allocation of nonseparable cost value

        # separable contribution
        EANS = zeros(NP,1)
        SCont = zeros(NP,1)

        for i=1:NP
            c = Int8[1 for i=1:NP]
            c[i] = 0
            ind_S_no_i = findall([sum([coalitions[k,j] .== c[j] for j=1:NP]) == NP for k=1:NC] .== true)

            SCont[i] = vN- v[ind_S_no_i[1]]
        end

        total_SC = sum(SCont)
        for i=1:NP
            c = Int8[0 for i=1:NP]
            c[i] = 1
            ind_S_no_i = findall([sum([coalitions[k,j] .== c[j] for j=1:NP]) == NP for k=1:NC] .== true)

            EANS[i] = SCont[i] + (1/NP)*(vN -total_SC)
        end

        CIS_EANS_size = size(ind_individual)[1] + NP + 1


    
        ## (CIS + EANS)*0.5
        CIS_EANS = (CIS + EANS)*0.5

        ## Matrix with values
        Mat = hcat(SV, CIS_EANS)

        total_size = size(v)[1]
        println("Total size:$total_size")
        println("Used size:$CIS_EANS_size")

        return Mat
    end

    function RunMC_Shapley(ϕ, v, coalitions, grid, requests, exp)
        println(size(coalitions))
        ind_EV_loc = sort(grid.ind_EV_loc)
        N = size(ind_EV_loc)[1] 
        M_sim = exp.m
        ϕ_per_iteration = ϕ
        for i = 1:M_sim
            idx = rand(1:size(coalitions)[1])
            indp = coalitions[idx]
    
            # first coalition
            idx = findall(v -> v == indp[1], ind_EV_loc)[1]
            ϕ[idx] = ((i - 1) / i) * ϕ[idx] + 
                (1 / i) * compute_value_function([indp[1]], v, ind_EV_loc, requests)
    
            # third to last
            for j = 2:N
                idx = findall(v-> v == indp[j], ind_EV_loc)[1]
                v_j = compute_value_function(indp[1:j], v, ind_EV_loc, requests)
                v_j_1 = compute_value_function(indp[1:j-1], v, ind_EV_loc, requests) # can be replaced by the previous it. value j
                marginal_value = v_j - v_j_1
                ϕ[idx] = ((i - 1) / i) * ϕ[idx] + (1 / i) * marginal_value
            end
    
            ϕ_per_iteration = [ϕ_per_iteration; ϕ]
        end
    
        return ϕ_per_iteration
    end

    function compute_value_function(coalition, v, ind_EV_loc, requests)
        coalition = sort(coalition)
        combs = [p for p in combinations(ind_EV_loc)]
        idx = findall(x -> coalition == x, combs)
        push!(requests, coalition)
        return v[idx][1]
    end

    function run_monte_carlo(grid, v, exp::Simulations.Experiment)
        N = size(grid.ind_EV_loc)[1] # Number of players
        ind_EV_loc = grid.ind_EV_loc
        coalitions = [p for p in multiset_permutations(ind_EV_loc, size(ind_EV_loc, 1))];

        requests = []
        ϕ = zeros(Float32, (1, N));

        res = RunMC_Shapley(ϕ, v, coalitions, grid, requests, exp);
        ϕ_average = zeros(Float32, (1, N))

        Φ_m = []
        N_exp = exp.n

        for i = 1:N_exp
            ϕ = zeros(Float32, (1, N))
            
            result = RunMC_Shapley(ϕ, v, coalitions, grid, requests, exp)
            push!(Φ_m, result)
            
            if i == 1
                ϕ_average = result
            else
                ϕ_average = 0.5 .* (result + ϕ_average)
            end

        end

        s = size(unique(requests), 1)
        t_s = size(v)[1]

        println("Used size:$s")
        println("Total size:$t_s")

        return ϕ_average[end,1:end]
        
    end
end
