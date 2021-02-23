abstract type AbstractRKTemporal <: AbstractTemporal end

function get_RK_scheme(scheme::String)
    append = occursin(r"\(([^()]|(?R))*\)", scheme) ? "" : "()"
    A_scheme = Meta.parse(string("A_", scheme, append))
    b_scheme = Meta.parse(string("b_", scheme, append))
    A = @eval $A_scheme
    b = @eval $b_scheme
    return A, vec(b), vec(sum(A, dims=2))
end


function amplification_factor(λ::N, X::T) where {T<:AbstractRKTemporal, N<:Number}
    1 + λ*transpose(X.b)*inv(UniformScaling(1) - λ*X.A)*ones(eltype(X.b), size(X.b))
end

@memoize function amplification_factor_poly(X::T) where {T<:AbstractRKTemporal}
    Λ, ~ = gausslegendre(length(X.b) + 1)
    P = amplification_factor.(Λ, Ref(X))
    fit(Λ, P).coeffs
end
function amplification_factor(λ::Matrix{N}, X::T) where {T<:AbstractRKTemporal, N<:Number}
    P_modes = amplification_factor_poly(X)
    matrixpolyval(P_modes, λ)
end

function pseudo_amplification_factor(λ::N, ω, dt, dτ, X::T) where {T<:AbstractRKTemporal, N<:Number}
    1 + dτ*(λ - 1/(dt*ω))*transpose(X.b)*inv(UniformScaling(1) - dτ*λ*X.A)*ones(eltype(X.b), size(X.b))
end

@memoize function pseudo_amplification_factor_poly(ω, dt, dτ, X::T) where {T<:AbstractRKTemporal}
    Λ, ~ = gausslegendre(length(X.b) + 1)
    P = pseudo_amplification_factor.(Λ, Ref(ω), Ref(dt), Ref(dτ), Ref(X))
    fit(Λ, P).coeffs
end

function pseudo_amplification_factor(Λ::Matrix{N}, ω, dt, dτ, X::T) where {T<:AbstractRKTemporal, N<:Number}
    P_modes = pseudo_amplification_factor_poly(ω, dt, dτ, X)
    matrixpolyval(P_modes, Λ)
end

function pseudo_source_factor(λ::N, ω, dτ, X::T) where {T<:AbstractRKTemporal, N<:Number}
    dτ*transpose(X.b)*inv(UniformScaling(1) - λ*dτ*X.A)*ones(eltype(X.b), size(X.b))
end

function pseudo_source_factor(Λ::Matrix{N}, ω, dτ, X::T) where {T<:AbstractRKTemporal, N<:Number}
    P_modes = pseudo_source_factor_poly(ω, dτ, X)
    matrixpolyval(P_modes, Λ)
end

@memoize function pseudo_source_factor_poly(ω, dτ, X::T) where {T<:AbstractRKTemporal}
    Λ, ~ = gausslegendre(length(X.b) + 1)
    P = pseudo_source_factor.(Λ, Ref(ω), Ref(dτ), Ref(X))
    fit(Λ, P).coeffs
end
