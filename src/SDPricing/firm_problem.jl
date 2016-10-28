
"""
### Fields
- `basis`   :
- `z_basis` :
- `p̃_basis`
- `Φ_tensor`  :
- `Φ`         :
- `Φ_fac`     :
"""
type CollocStruct{BF1<:BasisFamily,BF2<:BasisFamily}

    basis::CompEcon.Basis{2}
    z_basis::CompEcon.Basis{1,BF1}
    p̃_basis::CompEcon.Basis{1,BF2}

    Φ_tensor::CompEcon.BasisStructure{CompEcon.Tensor}

    Φ::Array{Float64,2}
    Φ_fac::Base.LinAlg.LU{Float64,Array{Float64,2}}

    # Nodes
    z_nodes::Vector{Float64}
    n_z::Int64
    p̃_nodes::Vector{Float64}
    n_p̃::Int64
    grid_nodes::Array{Float64,2}
    # coeffs
    coeff::Array{Float64,2}
    ξstar::Vector{Float64}
    # Indices
    ind_z_x_z::Array{Int64,2}
    ind_z_x_p̃::Array{Int64,2}
end

immutable FirmProblem{BF1<:BasisFamily,BF2<:BasisFamily}
    # Parameters
    β::Float64
    ϵ::Float64

    #Stochastic Process
    n_z::Int64
    Π_z::Matrix{Float64}
    z_vals::Vector{Float64}
    # Menu Cost
    ξbar::Float64
    H::Function
    cond_mean::Function

    # Basis Object
    mbasis::CollocStruct{BF1,BF2}
end

function FirmProblem{BF1<:BasisFamily,BF2<:BasisFamily}(::Type{BF1}, ::Type{BF2}, β = 0.2, ϵ = 5.0)

    #== Stochastic productivity ==#
    n_z = 10
    mc  = rouwenhorst(n_z, 0.9, 0.1)
    z_vals = exp(collect(mc.state_values))                                     ##  IMPORTANT:  nodes are the exponencial     ##
    Π_z    = mc.p

    #== Menu Cost ==#
    ξbar = 2.5
    H(ξ) = 1.0/( ξbar - 0.0 )*(ξ - 0.0)
    cond_mean(ξ) = 1/(2*ξbar) * (ξ.^2)

    #== Create Basis ==#
    n_p̃ = 30
    z_basis = Basis(BF1(), z_vals, 0, 1)

    p_low   = 0.5 / z_vals[end]
    p_high  = 1.1 * ϵ/(ϵ-1) * 0.5/ z_vals[1]
    p̃_basis = Basis(BF2(), n_p̃, log( p_low ), log( p_high ) )
    mbasis::CollocStruct{BF1,BF2} = CollocStruct(z_basis, p̃_basis)

    FirmProblem{BF1,BF2}(β, ϵ, n_z, Π_z, z_vals,
                         ξbar, H, cond_mean, mbasis)
end

function Base.show(io::IO, fp::FirmProblem)

    m = """
    FirmProblem

    """
    print(io, m)
    show(io, fp.mbasis.basis)
end

function CollocStruct{BF1<:BasisFamily,BF2<:BasisFamily,BP1,BP2}(z_basis::CompEcon.Basis{1,BF1,BP1}, p̃_basis::CompEcon.Basis{1,BF2,BP2})

    ##  IMPORTANT:  careful with construction     ##
    basis = z_basis × p̃_basis

    Φ_tensor = BasisStructure(basis,Tensor())

    Φ_direct = BasisStructure(basis,Direct())
    Φ        = convert(Expanded, Φ_direct).vals[1]
    Φ_fac    = factorize(Φ)

    #== GRIDS that I will use in the solution ==#
    z_nodes = nodes(z_basis)[1]
    p̃_nodes = nodes(p̃_basis)[1]
    grid_nodes = gridmake(z_nodes, p̃_nodes)

    n_z    = length(z_nodes)
    n_p̃    = length(p̃_nodes)

    raw_ind = [collect(1:n) for n in (n_z, n_z)]
    ind_z_x_z = gridmake(raw_ind...)

    raw_ind = [collect(1:n) for n in (n_z, n_p̃)]
    ind_z_x_p̃ = gridmake(raw_ind...)

    CollocStruct{BF1,BF2}(basis, z_basis, p̃_basis,
                         Φ_tensor, Φ, Φ_fac,
                         z_nodes, n_z, p̃_nodes, n_p̃, grid_nodes,
                         zeros(n_z*n_p̃,3), zeros(n_z*n_p̃,1), ind_z_x_z, ind_z_x_p̃)
end
