function matrixpolyval(p, X::Matrix{T}) where {T<:Number}
    l, W = eigen(X)
    L = Diagonal(l)
    Y = zeros(eltype(X), size(X)) + UniformScaling(p[1])
    for i=2:length(p)
        Y += p[i]*L^(i-1)
    end
    Y = W*Y*inv(W)
end

perms(a) = reverse(collect(permutations(a)))

function organise_cols!(A::Matrix{T}) where {T<:Number}
    A[2,:] = closest(A[1,:], A[2,:])
    for i=3:size(A, 1)
        y = 2A[i-1,:] - A[i-2,:] #linear extrapolation of last 2 cols
        A[i,:] = closest(y, A[i,:])
    end
end

function organise_rows!(A::Matrix{T}) where {T<:Number}
    A[:,2] = closest(A[:,1], A[:,2])
    for i=3:size(A, 2)
        y = 2A[:,i-1] - A[:,i-2] #linear extrapolation of last 2 cols
        A[:,i] = closest(y, A[:,i])
    end
end

function closest(a, b)
    c = b
    min_dif = sum(abs.(b - a))
    for p in permutations(b)
        dif = sum(abs.(p - a))
        if min_dif > dif
            min_dif = dif
            c = p
        end
    end
    return c
end
