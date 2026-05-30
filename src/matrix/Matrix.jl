module Matrix

using LatticeGeometry
using LatticeModel
using SparseArrays

include("../utils.jl")

include("builder.jl")
export build_local_matrix, build_onsite_matrix, build_pair_matrix, build_hamiltonian_matrix

end
