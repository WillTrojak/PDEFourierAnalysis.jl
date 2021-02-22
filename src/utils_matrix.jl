function matrixpolyval(p, X::Matrix{T}) where {T<:Number}
    l, W = eigen(X)
    L = Diagonal(l)
    Y = ones(eltype(X), size(X))*p[1]
    for i=2:length(p)
        Y += p[i]*L^(i-1)
    end
    Y = W*Y*inv(W)
end
