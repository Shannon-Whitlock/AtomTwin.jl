@testset "getqstate raises before simulation" begin
    g, e = Level("g"), Level("e")
    atom = Atom(; levels = [g, e])
    sys = System(atom)
    @test_throws ErrorException getqstate(sys)
end

@testset "play returns valid detector output" begin
    g, e = Level("g"), Level("e")
    atom = Atom(; levels = [g, e])
    sys = System(atom)

    coupling = add_coupling!(sys, atom, g => e, 2π * 1e6; active = false)
    add_detector!(sys, PopulationDetectorSpec(atom, e; name = "P_e"))

    seq = Sequence(1e-9)
    @sequence seq begin
        Pulse(coupling, 100e-9)
    end

    out = play(sys, seq; initial_state = g)
    @test haskey(out.detectors, "P_e")
    @test length(out.times) > 0
    @test all(x -> 0.0 - 1e-6 ≤ x ≤ 1.0 + 1e-6, out.detectors["P_e"])
end
