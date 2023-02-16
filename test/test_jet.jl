@testitem "JET stability" begin
    using JET

    function test(mode)
        latvecs = lattice_vectors(1,1,2,90,90,90)
        crystal = Crystal(latvecs, [[0,0,0]])
        L = 2
        sys = System(crystal, (L,L,1), [SpinInfo(1, S=1)], mode)
        set_exchange!(sys, -1.0, Bond(1,1,(1,0,0)))
        polarize_spins!(sys, (0,0,1))

        # Test stability of LocalSampler
        sampler = LocalSampler(kT=0.2, propose=propose_flip)
        @test_opt step!(sys, sampler)

        # Test stability with mixing
        propose = @mix_proposals 0.5 propose_flip 0.5 propose_delta(0.2)
        sampler = LocalSampler(kT=0.2; propose)
        @test_opt step!(sys, sampler)
    end

    test(:dipole)
    test(:SUN)
end
