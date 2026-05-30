function build_local_matrix(::Type{T}, num_sites::Int, ids_site::NTuple{K,Int}, matrices::NTuple{K,<:AbstractMatrix}) where K where T<:Number
    dim = 2 ^ num_sites
    dim_matrix = 2 ^ K

    rows = Int[]
    cols = Int[]
    vals = T[]

    for col in 0:(dim-1)
        col_local = 0
        for (n, id) in enumerate(ids_site)
            bit = (col >> (id - 1)) & 1
            col_local |= bit << (n - 1)
        end

        for row_local in 0:(dim_matrix-1)
            val = one(T)
            for (n, matrix) in enumerate(matrices)
                row_bit = (row_local >> (n - 1)) & 1
                col_bit = (col_local >> (n - 1)) & 1
                val *= matrix[row_bit+1, col_bit+1]
            end

            iszero(val) && continue

            row = col
            for (n, id) in enumerate(ids_site)
                row_bit = (row_local >> (n - 1)) & 1
                mask = 1 << (id - 1)
                row = (row & ~mask) | (row_bit << (id - 1))
            end

            push!(rows, row + 1)
            push!(cols, col + 1)
            push!(vals, val)
        end
    end

    return sparse(rows, cols, vals)
end

function build_local_matrix(::Type{T}, num_sites::Int, id_site::Int, matrix::AbstractMatrix) where T<:Number
    return build_local_matrix(T, num_sites, (id_site,), (matrix,))
end

function build_zero_matrix(::Type{T}, num_sites::Int) where T<:Number
    dim = 2 ^ num_sites
    return spzeros(T, (dim, dim))
end

function build_term_type(term::OnsiteTerm)
    return eltype(op_to_matrix(term.op))
end

function build_term_type(term::PairTerm)
    T1 = eltype(op_to_matrix(term.op1))
    T2 = eltype(op_to_matrix(term.op2))
    return promote_type(T1, T2)
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::OnsiteTerm) where T<:Number
    matrix = op_to_matrix(term.op)

    num_sites = nsites(lattice)
    matrix_result = build_zero_matrix(T, num_sites)

    for site in 1:nsites(lattice)
        matrix_result += term.coeff * build_local_matrix(T, num_sites, site, matrix)
    end

    return matrix_result
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::PairTerm) where T<:Number
    matrix1 = op_to_matrix(term.op1)
    matrix2 = op_to_matrix(term.op2)

    num_sites = nsites(lattice)

    matrix_result = build_zero_matrix(T, num_sites)
    for (site1, site2) in neighbor_pairs(lattice, term.shell)
        matrix_result += term.coeff * build_local_matrix(T, num_sites, (site1, site2), (matrix1, matrix2))
    end

    return matrix_result
end

function build_term_matrix(lattice::Lattice, term::AbstractTerm)
    T = build_term_type(term)
    return build_term_matrix(T, lattice, term)
end

function build_hamiltonian_matrix(lattice::Lattice, model::GenericModel)
    num_sites = nsites(lattice)
    T = promote_type(build_term_type.(model.terms)...)
    matrix_result = build_zero_matrix(T, num_sites)

    for term in model.terms
        matrix_result += build_term_matrix(T, lattice, term)
    end

    return matrix_result
end
