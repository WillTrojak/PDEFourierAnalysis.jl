abstract type RKTemporal <: AbstractTemporal end

function get_RK_scheme(scheme::String)
    append = occursin(r"\(([^()]|(?R))*\)", scheme) ? "" : "()"
    A_scheme = Meta.parse(string("A_", scheme, append))
    b_scheme = Meta.parse(string("b_", scheme, append))
    A = @eval $A_scheme
    b = @eval $b_scheme
    return A, vec(b), vec(sum(A, dims=2))
end

function amplification_factor(lambda::N, X::T) where {T<:RKTemporal, N<:Number}
    e = ones(eltype(X.b), size(X.b))
    I = UniformScaling(1)
    1 + lambda*transpose(X.b)*inv(I - lambda*X.A)*e
end

function amplification_factor(lambda::Matrix{N}, X::T) where {T<:RKTemporal, N<:Number}
    P_modes = amplification_factor_poly(X)
    matrixpolyval(P_modes, lambda)
end

@memoize function amplification_factor_poly(X::T) where {T<:RKTemporal}
    Lambda, ~ = gausslegendre(length(X.b) + 1)
    P = amplification_factor.(Lambda, Ref(X))
    Poly.fit(Lambda, P).coeffs
end

function pseudo_amplification_factor(lambda::N, omega, dt, dtau, X::T) where {T<:RKTemporal, N<:Number}
    e = ones(eltype(X.b), size(X.b))
    I = UniformScaling(1)
    1 + dtau*(lambda - 1/(dt*omega))*transpose(X.b)*inv(I - lambda*X.A)*e
end

#function pseudo_amplification_factor(lambda::Matrix{N}, omega, dt, dtau, X::T) where {T<:RKTemporalScheme, N<:Number}
#
#end

function pseudo_source_factor(lambda::N, omega, dtau, X::T) where {T<:RKTemporal, N<:Number}
    e = ones(eltype(X.b), size(X.b))
    I = UniformScaling(1)
    dtau*transpose(X.b)*inv(I - lambda*dtau*X.A)*e
end
