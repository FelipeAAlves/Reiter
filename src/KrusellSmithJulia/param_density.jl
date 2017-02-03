


"""
Evals g_ϵ ( a; 1, ρ ), where

    g_ϵ(a; 1, ρ) = exp \{ ρ¹(a-m¹) + ∑ ρⁱ [ (a-m¹)ⁱ - mⁱ] \}

#### Return

- `dens::Vector{T}` : density value at {a_j} points
"""
function density{T<:Real}(ρ::Vector{T}, moments::Vector{T}, asset_points::Vector{Float64}, cp::ConsumerProblem)

    nAssets  = length(asset_points)
    nMoments = length(moments)

    aGridMoments = Float64[
                    (asset_points[j] - moments[1])^i - (i>1)*moments[i]
                    for j=1:nAssets, i=1:nMoments]

    dens = exp(aGridMoments * ρ)

    return dens
end

"""
Returns the value of
        ∫ g_ϵ ( a; 1, ρ ) * (a - m)ⁱ da
where
    g_ϵ(a; 1, ρ) = exp \{ ρ¹(a-m¹) + ∑ ρⁱ [ (a-m¹)ⁱ - mⁱ] \}

### INPUT
    - `ρ`                   : parameters of the distribution
    - `moments`
    - `cp::ConsumerProblem`
    - `iMom`                :

NOTE : for `iMom`=0 it just evaluates the density.
"""
function density_int{T<:Real}(ρ::Vector{T}, moments::Vector{T}, cp::ConsumerProblem, iMom::Int64 = 0, g0::Float64 = 1.0)

    nAssets    = cp.nQuad
    nMoments = length(moments)

    asset_quad_points, quad_weights  =  cp.asset_quad_points, cp.quad_weights
    AssetQuadGrid = Float64[
                    (asset_quad_points[j] - (iMom>1)*moments[1])^iMom
                    for j=1:nAssets]

    asset_quad_points, quad_weights  =  cp.asset_quad_points, cp.quad_weights

    return dot( quad_weights, g0 * density(ρ, moments, asset_quad_points, cp) .* AssetQuadGrid )
end
"""

### INPUT
- `ρ`       : point to eval
- `∇dens`   : gradient
"""
function density_int_g!{T<:Real}(ρ::Vector{T}, ∇dens::Vector{T}, moments::Vector{T}, cp::ConsumerProblem)

    nAssets    = cp.nQuad
    nMoments = length(moments)

    asset_quad_points, quad_weights  =  cp.asset_quad_points, cp.quad_weights
    aGridMoments = Float64[
                    (asset_quad_points[j] - moments[1])^i - (i>1)*moments[i]
                    for j=1:nAssets, i=1:nMoments]

    dens = density(ρ, moments, asset_quad_points, cp)

    for iMom = 1 : nMoments
        ∇dens[iMom] = dot( quad_weights , dens .* aGridMoments[:,iMom] )
    end

    return Void
end
