abstract type RKTemporalScheme <: TemporalScheme end

function get_RK_scheme(scheme::String)
    append = occursin(r"\(([^()]|(?R))*\)", scheme) ? "" : "()"
    A_scheme = Meta.parse(string("A_", scheme, append))
    b_scheme = Meta.parse(string("b_", scheme, append))
    A = @eval $A_scheme
    b = @eval $b_scheme
    return A, vec(b), vec(sum(A, dims=2))
end

function amplification_factor(lambda, X::T) where {T<:RKTemporalScheme}
    e = ones(eltype(X.b), size(X.b))
    I = LA.UniformScaling(1)
    1 + lambda*transpose(X.b)*inv(I - lambda*X.A)*e
end

function pseudo_amplification_factor(lambda::N, omega, dt, dtau, X::T) where {T<:RKTemporalScheme, N<:Number}
    e = ones(eltype(X.b), size(X.b))
    I = LA.UniformScaling(1)
    1 + dtau*(lambda - 1/(dt*omega))*transpose(X.b)*inv(I - lambda*X.A)*e
end

function pseudo_amplification_factor(lambda::Matrix{N}, omega, dt, dtau, X::T) where {T<:RKTemporalScheme, N<:Number}

end

function pseudo_source_factor(lambda, omega, dtau, X::T) where {T<:RKTemporalScheme}
    e = ones(eltype(X.b), size(X.b))
    I = LA.UniformScaling(1)
    dtau*transpose(X.b)*inv(I - lambda*dtau*X.A)*e
end
