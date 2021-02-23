abstract type AbstractFullScheme <: AbstractScheme end

mutable struct FullScheme <: AbstractFullScheme
    X::AbstractSpatial
    T::AbstractTemporal
    function FullScheme(X::Tx, T::Tt) where {Tx<:AbstractSpatial,Tt<:AbstractTemporal}
        new(X, T)
    end
end

CFL_limit(F::T; kwargs...) where {T<:AbstractFullScheme} = CFL_limit(F.X, F.T; kwargs...)

function CFL_limit(X::Tx, T::Tt; kwargs...) where {Tx<:AbstractFRSpatial,Tt<:AbstractERKTemporal}
    function max_spectral_radius(t, K, X, T)
        rho_max = 0.
        for k in K
            rho_max = max(rho_max, maximum(abs.(eigvals(amplification_factor(t*Qmatrix(X, k), T)))))
        end
        return rho_max
    end

    t0 = haskey(kwargs, :t_start) ? getindex(kwargs, :t_start) : 1e-5
    t1 = haskey(kwargs, :t_end) ? getindex(kwargs, :t_end) : 20

    k_min = haskey(kwargs, :k_min) ? getindex(kwargs, :k_min) : 1e-5
    k_max = haskey(kwargs, :k_max) ? getindex(kwargs, :k_max) : nyquist_wavenumber(X)
    n_k = haskey(kwargs, :n_k) ? getindex(kwargs, :n_k) : 100
    K = LinRange(k_min, k_max, n_k)

    t_tol = haskey(kwargs, :t_tol) ? getindex(kwargs, :t_tol) : 1e-5
    r_tol = haskey(kwargs, :r_tol) ? getindex(kwargs, :r_tol) : 1e-5

    while t1-t0 > t_tol
        t2 = (t1 + t0)/2
        rho_2 = max_spectral_radius(t2, K, X, T)
        t1 = rho_2 >= 1+r_tol ? t2 : t1
        t0 = rho_2 >= 1+r_tol ? t0 : t2
    end

    return (t1 + t0)/2
end
