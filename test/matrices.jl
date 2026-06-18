function test_matrix_1()
    space = Space(
        SpinHalfSpace(),
        Hypercubic(
            (2, 2),
            OpenBoundary(2),
        ),
    )
    model = TFIHamiltonian(space, -1.0, -2.0)

    hamiltonian_matrix = build_matrix(space, model)
    hamiltonian_matrix_true = Float64[
        -4 -2 -2 0 -2 0 0 0 -2 0 0 0 0 0 0 0
        -2 0 0 -2 0 -2 0 0 0 -2 0 0 0 0 0 0
        -2 0 0 -2 0 0 -2 0 0 0 -2 0 0 0 0 0
        0 -2 -2 0 0 0 0 -2 0 0 0 -2 0 0 0 0
        -2 0 0 0 0 -2 -2 0 0 0 0 0 -2 0 0 0
        0 -2 0 0 -2 0 0 -2 0 0 0 0 0 -2 0 0
        0 0 -2 0 -2 0 4 -2 0 0 0 0 0 0 -2 0
        0 0 0 -2 0 -2 -2 0 0 0 0 0 0 0 0 -2
        -2 0 0 0 0 0 0 0 0 -2 -2 0 -2 0 0 0
        0 -2 0 0 0 0 0 0 -2 4 0 -2 0 -2 0 0
        0 0 -2 0 0 0 0 0 -2 0 0 -2 0 0 -2 0
        0 0 0 -2 0 0 0 0 0 -2 -2 0 0 0 0 -2
        0 0 0 0 -2 0 0 0 -2 0 0 0 0 -2 -2 0
        0 0 0 0 0 -2 0 0 0 -2 0 0 -2 0 0 -2
        0 0 0 0 0 0 -2 0 0 0 -2 0 -2 0 0 -2
        0 0 0 0 0 0 0 -2 0 0 0 -2 0 -2 -2 -4
    ]
    @test size(hamiltonian_matrix) == (16, 16)
    @test hamiltonian_matrix == hamiltonian_matrix_true
end

function test_matrix_2()
    t = 1.0
    u = 2.0

    space = Space(
        SpinfulFermionSpace(),
        Hypercubic(
            (2,),
            OpenBoundary(1),
        ),
    )
    model = HubbardHamiltonian(space, t, u)

    hamiltonian_matrix = build_matrix(space, model)


    hamiltonian_matrix_true = zeros(ComplexF64, 16, 16)
    hamiltonian_matrix_true[4, 4] += u
    hamiltonian_matrix_true[8, 8] += u
    hamiltonian_matrix_true[12, 12] += u
    hamiltonian_matrix_true[13, 13] += u
    hamiltonian_matrix_true[14, 14] += u
    hamiltonian_matrix_true[15, 15] += u
    hamiltonian_matrix_true[16, 16] += 2u

    hamiltonian_matrix_true[2, 5] += -t
    hamiltonian_matrix_true[5, 2] += -t

    hamiltonian_matrix_true[3, 9] += -t
    hamiltonian_matrix_true[9, 3] += -t

    hamiltonian_matrix_true[4, 7] += t
    hamiltonian_matrix_true[7, 4] += t

    hamiltonian_matrix_true[4, 10] += -t
    hamiltonian_matrix_true[10, 4] += -t

    hamiltonian_matrix_true[7, 13] += t
    hamiltonian_matrix_true[13, 7] += t

    hamiltonian_matrix_true[8, 14] += t
    hamiltonian_matrix_true[14, 8] += t

    hamiltonian_matrix_true[10, 13] += -t
    hamiltonian_matrix_true[13, 10] += -t

    hamiltonian_matrix_true[12, 15] += t
    hamiltonian_matrix_true[15, 12] += t

    hamiltonian_matrix_true += -I(16) * u * 2.0 / 4

    @test size(hamiltonian_matrix) == (16, 16)
    @test hamiltonian_matrix == hamiltonian_matrix_true
end

@testset "Matrix builder functions" begin
    test_matrix_1()
    test_matrix_2()
end
