module LatticeBackends

include("matrices/Matrices.jl")
using .Matrices
export State, apply, apply!, dense_matrix, sparse_matrix, reverse_bits

end
