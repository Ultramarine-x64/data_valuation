module Simulations
    using CSV
    using DataFrames
    using Setfield
    using Combinatorics

    Base.@kwdef struct GridSettings
        grid_size::Int32
        case_id::Int32
        alpha::Float64
        factor_line_capacity_scale::Float64
        ind_EV_loc::Vector{Int64}
        demand::NamedTuple
        c_ls::Int32
        c_es::Int32
        c_BM⁺::Matrix{Float64}
        c_BM⁻::Matrix{Float64}
        solver::Int8 = 2 # 1 for LinDistFlow, 2 for DistFlow
    end

    struct Profile
        branch_file::String  # Branch data
        nodes_file::String # Nodes data
        grid_settings::GridSettings
    end

    Base.@kwdef struct Grid
        NI::Int8
        NL::Int8
        NodesData::DataFrame
        Nodes::Matrix{Int8}
        R::Vector{Float64}
        X::Vector{Float64}
        Pdd::Matrix{Float64}
        Qdd::Matrix{Float64}
        Pmin::Vector{Float64}
        Pmax::Vector{Float64}
        Qmin::Vector{Float64}
        Qmax::Vector{Float64}
        Wmin::Float64
        Wmax::Float64
        SLmax::Vector{Float64}
        c_op::Matrix{Float64}
        c_ls::Int32
        c_es::Int32
        c_BM⁺::Matrix{Float64}
        c_BM⁻::Matrix{Float64}
        DA_optimal::Vector{Float64} = vec([])
    end

    Base.@kwdef struct Experiment
        m::Int32 = 0 # number of monte carlo simulations
        n::Int32 = 0 # number of monte carlo runs for averaging
    end

    function init(grid::GridSettings)
        # -- Source data file
        case = grid.grid_size
        branch_file = "nodes/Branch_$case.csv"
        nodes_file = "nodes/Nodes_$case.csv"

        @assert(size(grid.demand.wo_EV, 1) == case)

        return Profile(branch_file, nodes_file, grid)
    end

    function load_data_from_file(s::Profile)
        ## Data from Lines
        BranchData = DataFrame(CSV.File(s.branch_file))

        NBuses = max(maximum(BranchData[:,:FromBus]), maximum(BranchData[:,:ToBus]));
        NLines = size(BranchData,1);

        BranchData.Z = BranchData.R + im*BranchData.X;
        aux = 1. / BranchData.Z
        BranchData.Y = aux';
        Ybus = zeros(Complex,NBuses,NBuses); #adjacenсy matrix

        #F or notational convenience, we introduce the concepts of parent nodes and the adjacency matrix.
        # The parent of j, ρ(j), is the adjacent node in the
        # direction toward node 0. The adjacency matrix, Ybus, encodes the network topology.

        Nodes = zeros(Int8,NLines,2); #node's matrix
        SLmax = zeros(NLines); #line limits

        for i=1:NLines
            Ybus[BranchData[i,:FromBus],BranchData[i,:ToBus]]=-BranchData[i,:Y]
            Ybus[BranchData[i,:ToBus],BranchData[i,:FromBus]]=-BranchData[i,:Y]
            Nodes[i,1] = BranchData[i,:FromBus]
            Nodes[i,2] = BranchData[i,:ToBus]
            SLmax[i] = BranchData[i,:Smax]
        end

        for i=1:NBuses
            for j=1:NLines
                if (i==BranchData[j,:FromBus])||(i==BranchData[j,:ToBus])
                    Ybus[i,i]+=BranchData[j,:Y]
                end
            end
        end

        Gbus = real(Ybus);
        Bbus = imag(Ybus);
        R = BranchData.R;
        X = BranchData.X;

        SLmax = SLmax.*s.grid_settings.factor_line_capacity_scale

        ## Data from Nodes
        NodesData = DataFrame(CSV.File(s.nodes_file))

        # Generator limits
        Pmin = zeros(NBuses);
        Pmax = zeros(NBuses);
        Qmin = zeros(NBuses);
        Qmax = zeros(NBuses);
        # Demand
        Pd = zeros(NBuses);
        Qd = zeros(NBuses);
        # Generation costs
        c = zeros(NBuses);

        #battery
        PmaxChar = zeros(NBuses);
        PmaxDisc = zeros(NBuses);
        Emin = zeros(NBuses);
        Emax = zeros(NBuses);
        efiChar = zeros(NBuses);
        efiDisc = zeros(NBuses);

        for i=1:NBuses
            Pmin[i] = NodesData[i,:Pmin]
            Pmax[i] = NodesData[i,:Pmax]
            Qmin[i] = NodesData[i,:Qmin]
            Qmax[i] = NodesData[i,:Qmax]
            Pd[i] = NodesData[i,:Pd]
            Qd[i] = NodesData[i,:Qd]
            c[i] = NodesData[i,:c]
            PmaxChar[i] = NodesData[i,:PmaxChar]
            PmaxDisc[i] = NodesData[i,:PmaxDisc]
            Emin[i] = NodesData[i,:Emin]
            Emax[i] = NodesData[i,:Emax]
            efiChar[i] = NodesData[i,:efiChar]
            efiDisc[i] = NodesData[i,:efiDisc]
        end


        #Angle limits
        θmax = pi/2;
        θmin = -pi/2;

        # Set size
        NI = length(Pmax); # nodes aka buses
        NL = length(SLmax); # edges aka lines

        # Voltage limits
        vmin= 0.95;
        vmax= 1.05;

        Wmin= vmin.^2
        Wmax= vmax.^2


        ## Load Profile
        a = [1.]
        T = length(a)

        Pdd=zeros(NI, T);
        Qdd=zeros(NI, T);

        for i = 1:NI
            for t=1:T
                Pdd[i,t]= NodesData[i,:Pd]*a[t]
            end
        end

        for i = 1:NI
            for t=1:T
            Qdd[i,t]= NodesData[i,:Qd]*a[t]
            end
        end

        ## Price Profile
        b = copy(a)
        b[:] = 0.5 .+ b[:].*0.5
        c_op=repeat(c, 1,T)

        for t=1:T
            c_op[1,t] = c[1]*b[t]
        end

        grid = Grid(
            NI = NI,
            NL = NL,
            NodesData = NodesData,
            Nodes = Nodes,
            R = R,
            X = X,
            Pdd = Pdd,
            Qdd = Qdd,
            Pmin = Pmin,
            Pmax = Pmax,
            Qmin = Qmin,
            Qmax = Qmax,
            Wmin = Wmin,
            Wmax = Wmax,
            SLmax = SLmax,
            c_op = c_op,
            c_ls = s.grid_settings.c_ls,
            c_es = s.grid_settings.c_es,
            c_BM⁺ = s.grid_settings.c_BM⁺,
            c_BM⁻ = s.grid_settings.c_BM⁻
        )

        return grid
    end

    function update_data(grid::Grid)
        NodesData = grid.NodesData
        NI = grid.NI

        ## Load Profile
        a = [1.]
        T = length(a)

        Pdd=zeros(NI, T);

        for i = 1:NI
            for t=1:T
                Pdd[i,t]= NodesData[i,:Pd]*a[t]
            end
        end

        grid = @set grid.Pdd = Pdd

        return grid
    end

    function generate_coalitions(grid::GridSettings)
        m_coalitions = zeros(Int, 1, grid.grid_size)
        for p in combinations(grid.ind_EV_loc)
            dummy = zeros(Int, 1, grid.grid_size)
            dummy[p] .= 1 
            m_coalitions = [m_coalitions; dummy]
        end

        return m_coalitions
    end

end # Module