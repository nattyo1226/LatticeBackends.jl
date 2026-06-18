function build_matrix(::Type{N}, pr::AbstractOperatorPrimitive) where {N<:Number}
    throw(ArgumentError("Unsupported operator type: $(typeof(pr))"))
end

function build_matrix(pr::AbstractOperatorPrimitive)
    return build_matrix(ComplexF64, pr)
end

function build_matrix(::Type{N}, ::Identity{T}) where {N<:Number,T<:AbstractSystemTag}
    return reshape(N[1, 0, 0, 1], (2, 2))
end

function build_matrix(::Type{N}, ::PauliX) where {N<:Number}
    return reshape(N[0, 1, 1, 0], (2, 2))
end

function build_matrix(::Type{N}, ::PauliY) where {N<:Number}
    return reshape(N[0, -im, im, 0], (2, 2))
end

function build_matrix(::Type{N}, ::PauliZ) where {N<:Number}
    return reshape(N[1, 0, 0, -1], (2, 2))
end

function build_matrix(::Type{N}, pr::SummedOperatorPrimitive) where {N<:Number}
    return sum(build_matrix(N, pr_sub) for pr_sub in pr.prs)
end

function build_matrix(::Type{N}, pr::ProductedOperatorPrimitive) where {N<:Number}
    return foldl(*, [build_matrix(N, pr_sub) for pr_sub in pr.prs]; init=I)
end
