# Classical motion regression test: trap frequency of atom in a Gaussian tweezer.
#
# Validates that a Ytterbium-171 atom displaced from the beam centre oscillates
# at the analytically expected trap frequency.  Fails if the force calculation,
# polarizability model, or classical integrator is broken.
#
# Pattern follows atom_sorting.jl: classical-only dynamics via play(sys, seq)
# without an initial quantum state.

@testset "Classical motion: trap oscillation at correct frequency" begin
    c_v  = 2.997_924_58e8        # m/s
    ε0_v = 8.854_187_812_8e-12   # F/m
    amu  = 1.660_539_066_60e-27  # kg
    m    = 171 * amu             # Yb-171 mass

    # Tweezer parameters
    λ_nm = 759.0
    λ    = λ_nm * 1e-9
    w0   = 2e-6    # 2 µm waist
    P    = 0.1     # 100 mW

    # Analytical trap frequency
    α_SI   = polarizability_si(AtomTwin.YB171_POLARIZABILITY_1S0, λ_nm)
    I0     = 2 * P / (π * w0^2)
    ω_trap = sqrt(4 * α_SI * I0 / (m * c_v * ε0_v * w0^2))
    T_trap = 2π / ω_trap

    # Initial displacement: 10% of waist in x
    x0 = 0.10 * w0

    # Single-level atom (no quantum dynamics) displaced from the beam centre
    atom = Ytterbium171Atom(;
        levels = [Level("1S0")],
        x_init = [x0, 0.0, 0.0],
        v_init = [0.0, 0.0, 0.0],
    )
    tweezer = GaussianBeam(λ, w0, P)
    sys = System(atom, tweezer)
    add_detector!(sys, MotionDetectorSpec(atom; dims = [1], name = "x"))

    # Simulate for 3 trap periods using Wait (no quantum instruction needed)
    T_sim = 3 * T_trap
    dt    = min(T_trap / 200, 100e-9)
    seq   = Sequence(dt)
    @sequence seq begin
        Wait(T_sim)
    end

    # No initial_state → classical dynamics only, matching atom_sorting pattern.
    # Explicitly expect (and capture) the resulting warning so it doesn't appear as noise.
    out = @test_logs (:warn, r"Initial state not specified") play(sys, seq)
    x_traj = out.detectors["x"][:, 1]
    t_traj = out.times

    # 1. Atom must stay trapped (not escape)
    @test maximum(abs.(x_traj)) < 2 * x0

    # 2. First zero crossing must exist: x(t)=x0·cos(ω_trap·t) → zero at T_trap/4
    sign0 = sign(x_traj[1])
    first_zero_idx = findfirst(i -> sign(x_traj[i]) != sign0, 2:length(x_traj))
    @test !isnothing(first_zero_idx)

    if !isnothing(first_zero_idx)
        t_first_zero = t_traj[first_zero_idx]
        ω_measured   = π / (2 * t_first_zero)   # quarter-period → full frequency

        # 3. Measured frequency within 10% of analytical value
        @test isapprox(ω_measured, ω_trap, rtol = 0.10)
    end
end
