abstract type FRSpatialScheme <: SpatialScheme end

Dmatrix(X::T) where {T<:FRSpatialScheme} = [legendre_d(j-1, i-1) for i=1:X.p+1, j=1:X.p+1]

P2Xmatrix(X::T) where {T<:FRSpatialScheme} = [SPoly.basis.(SPoly.Legendre, i)(X.xs[j]) for i=0:X.p, j=1:X.p+1]
X2Pmatrix(X::T) where {T<:FRSpatialScheme} = inv(P2Xmatrix(X))
