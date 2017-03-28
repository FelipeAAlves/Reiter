
"""_
Holds information on collocation structure of the problem

### Fields
- `basis`   :

- `z_basis` :
- `p_basis` :

- `Φ_tensor`  :
- `Φ`         :
- `Φ_fac`     :

- `p_nodes` : in logs
- `n_coll_p`
- `z_nodes`
- `n_coll_z`
- `grid_nodes`
"""
immutable FirmColloc{BF<:BasisFamily,BP<:BasisParams}

    ##  Basis Information  ##
    p_basis::BasisMatrices.Basis{1,BF,BP}
    z_basis::BasisMatrices.Basis{1,BasisMatrices.Spline,BasisMatrices.SplineParams}
    basis::BasisMatrices.Basis{2}

    Φ::Union{SparseMatrixCSC{Float64,Int64},Array{Float64,2}}
    Φ_tensor::BasisMatrices.BasisMatrix{BasisMatrices.Tensor}
    Φ_fac::Union{Base.SparseArrays.UMFPACK.UmfpackLU{Float64,Int64},Base.LinAlg.LU{Float64,Array{Float64,2}}}

    ##  Nodes INFO  ##
    p_nodes::Vector{Float64}
    n_coll_p::Int64
    z_nodes::Vector{Int64}
    n_coll_z::Int64
    grid_nodes::Array{Float64,2}

end

function Base.show{BF1,BF2}(io::IO, fcoll::FirmColloc{BF1, BF2})

    m = """
    FirmColloc{$BF1,$BF2}

    """
    print(io, m)
    show(io, fcoll.basis)
end

function FirmColloc{BF<:BasisFamily}(::Type{BF})

    ##  WARN:  CARAFUL with construction     ##
    ##   gridmake(p_nodes, z_nodes) will generate
    ##      [ p_1 z_1]
    ##      [ p_2 z_1]
    ##      [ p_3 z_1]
    ##      [ p_1 z_2]
    ##      [ p_2 z_2]
    ##      [ p_3 z_2]
    @getPar __pars
    println(β)

    #== Collocation Grids ==#
    n_p = 30 # _default_ 30
    w_guess_high = 0.95; # _default_ 0.90
    w_guess_low  = 0.50; # _default_ -.70 _OK_
    p_low   = w_guess_low / z_vals[end]
    p_high  = 1.05 * ϵ/(ϵ-1) * (w_guess_high/ z_vals[1])

    #== More points in the middle ==#
    # p_grid_begin  = collect( linspace( p_low , 0.8   , 5) )
    # p_grid_end    = collect( linspace( 1.20  , p_high, 5) )
    # p_grid_middle = collect( linspace( 0.80  , 1.20  , n_p-8) )
    # p_grid = vcat(p_grid_begin[1:end-1], p_grid_middle, p_grid_end[2:end])

    #== Usual grid ==#
    p_grid = collect( linspace( p_low , p_high, n_p) )

    ##  WARN WARN:  nodes in p̃ are in logs  ##
    p_basis = Basis(BF(), log(p_grid), 0, 3) # using cubic spline
    z_basis = Basis(Spline(), z_vals, 0, 1)
    basis   = Basis(p_basis, z_basis)
    # .....................................................................................

    #== Useul objects ==#
    Φ_tensor = BasisMatrix(basis, Tensor())
    Φ_direct = BasisMatrix(basis, Direct())
    Φ        = convert(Expanded, Φ_direct).vals[1]
    Φ_fac    = factorize(Φ)

    #== GRIDS that I will use in the solution ==#
    p_nodes = nodes(p_basis)[1]
    z_nodes = collect(1:length(z_vals))
    grid_nodes = gridmake(p_nodes, z_nodes)

    n_coll_p    = length(p_nodes)
    n_coll_z    = length(z_nodes)

    FirmColloc(p_basis, z_basis, basis,
    Φ, Φ_tensor, Φ_fac,
    p_nodes, n_coll_p, z_nodes, n_coll_z, grid_nodes)
end

"""
Intends to hold the solution of the PROBLEM
### Fields
- `coeff` : coeff for the value function
- `pstar` : optimal price in case of adjustment
- `ξstar` :
"""
immutable FirmSolution
    ##  STEADY-STATE solutions  ##
    coeff::Array{Float64,2}
    pstar::Vector{Float64}
    ξstar::Vector{Float64}

    ##  Solutions at other points  ##
    pol_at::Array{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}},1} # DEPRECATED

end

function FirmSolution(fcoll::FirmColloc)

    FirmSolution( zeros(Float64, fcoll.n_coll_p*fcoll.n_coll_z, 3),
                  zeros(fcoll.n_coll_z), zeros(length(fcoll.z_basis)*length(fcoll.p_basis)),
                  Array(Tuple{Vector{Float64},Vector{Float64},Vector{Float64}},0) )

end


"""
    Holds information on the steady-state
"""
# Observation: Let it be type and not immutable to facilitate UPDATE in steady-state fnc
type StstHistogram

    ## Grid ##
    p_hist_nodes::Vector{Float64}
    z_hist_nodes::Vector{Int64}
    hist_nodes::Matrix{Float64}

    ##  Distribution ##
    vHistogram::Vector{Float64}
    mHistogram::Matrix{Float64}

    ##  All Ss values and Dictionary  ##
    xstst::Vector{Float64}
    ystst::Vector{Float64}
    ixvar::Dict{Symbol,UnitRange{Int64}}
    iyvar::Dict{Symbol,UnitRange{Int64}}

    Zstst::Vector{Float64}
    iZvar::Dict{Symbol,UnitRange{Int64}}

    ##  Auxiliary  ##
    χ::Float64
end


function StstHistogram(fcoll::FirmColloc)

    @getPar __pars

    p_nodes  = fcoll.p_nodes
    n_coll_p = fcoll.n_coll_p

    #== Grids ==#
    p_hist_nodes = linspace( p_nodes[1],p_nodes[end], 4*n_coll_p) # _default_ 3* n_coll, increase if having trouble finding the stst
    z_hist_nodes = 1:length(z_vals)
    hist_nodes = gridmake(p_hist_nodes, z_hist_nodes)

    nHistogram = length(p_hist_nodes) * n_z
    vHistogram = zeros(nHistogram)
    mHistogram = zeros(length(p_hist_nodes), n_z)

    ###  Create the Dicitionary for x,y ###
    # *************************************************************************************
    ixvar = Dict{Symbol,UnitRange{Int64}}()
    iyvar = Dict{Symbol,UnitRange{Int64}}()

    ##  CASE:  State variables  ##  [vHistogram[2:end], Z, ϵ]
    nx = 0
    # ixvar[:endo_aggr]  = nn:5;                                 nn += 5;
    ixvar[:histogram] = nx+1:nx+(nHistogram-1);   nx += (nHistogram-1)  ###  WARN : histogram - 1
    # ixvar[:histogram] = nx+1:nx+(nHistogram);   nx += (nHistogram)    ###  WARN : last element treatment histogram
    ixvar[:exog_aggr] = nx+1:nx+2;              nx += 2;                ###  [Z_t, ϵ_i]

    ##  CASE:  controls  ## [Y_t, N_t, i_t, Π_t, w_t, Ve, pᵃ]
    ny = 0
    iyvar[:aggr]      = ny+1:5;                            ny += 5;   # [Y_t, N_t, i_t, Π_t, w_t]
    iyvar[:value_fnc] = ny+1:ny+(n_coll_p*n_z);            ny += (n_coll_p*n_z)
    iyvar[:foc]       = ny+1:ny+(n_z);                     ny += (n_z)
    # -------------------------------------------------------------------------------------

    ###  Create the Dicitionary for Z  ###
    # **************************************************************************************
    # Zstst = [ystst; xstst; ystst; xstst; 0; 0]
    #
    iZvar = Dict{Symbol, UnitRange{Int64}}()

    nExog = 2
    # nEta  = length(fcoll.grid_nodes[:,1]) + length(ss_histogram.vHistogram) + n_z

    #== iZvar ==#
    iZvar[:x′]    = 1:nx
    iZvar[:y′]    = nx+1:ny+nx

    iZvar[:y] = (nx+ny) + (iZvar[:y′])
    iZvar[:x] = (nx+ny) + (iZvar[:x′])
    nn = 2*(nx+ny)

    iZvar[:eps] = nn+1:nn+nExog
    # -------------------------------------------------------------------------------------

    StstHistogram(p_hist_nodes, z_hist_nodes, hist_nodes, vHistogram, mHistogram,
        Array(Float64,nx), Array(Float64,ny), ixvar, iyvar,
        Array(Float64,nn+nExog),iZvar, 0.0)
end

function BasisMatrices.nodes(ss_histogram::StstHistogram)
    return ss_histogram.hist_nodes, (ss_histogram.p_hist_nodes, ss_histogram.z_hist_nodes)
end
