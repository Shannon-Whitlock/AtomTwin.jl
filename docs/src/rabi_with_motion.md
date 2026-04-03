# Rabi oscillations with atomic motion

In this example we simulate Rabi oscillations of a single atom trapped in a
tweezer, including its thermal motion in the trap.

The atom is initialized with a Maxwell–Boltzmann velocity distribution, and
we monitor both the internal excited–state population and the motional
degrees of freedom during a resonant Rabi pulse.

````julia
using AtomTwin
using AtomTwin.Units
using StatsBase
using Plots
````

## Parameters

````julia
Ω              = 2π * 0.05MHz   # Rabi frequency (rad/s)
temperature    = 5µK           # Initial temperature (K)

pulse_duration = 200µs          # Pulse duration (s)
dt             = 10ns            # Time step (s)
````

## System construction

We construct a simple two–level atom (|g⟩, |e⟩) and load it into a single
tweezer site, with motional degrees of freedom determined by the trap
parameters and the thermal velocity distribution.

Two-level atom with thermal velocity distribution

````julia
g, e = Level("1S0"), Level("3P0")
atom = Ytterbium171Atom(;
    levels = [g, e],
    x_init = [0.0, 0.0, 0.0],
    v_init = maxwellboltzmann(T = temperature),
)

display(atom)
````

Single-site tweezer array with specified geometry and powers

````julia
tweezer = GaussianBeam(
    λ    = 759nm,
    w0   = 1.0µm,
    P   = 50mW
)
````

Resonant planar beam driving the |g⟩ ↔ |e⟩ transition

````julia
beam = PlanarBeam(578e-9, 1.0, [1.0, 0.0, 0.0], [0, 1, 0])
````

Build the full system and add a coherent coupling between |g⟩ and |e⟩

````julia
system   = System(atom, tweezer)
coupling = add_coupling!(system, atom, g => e, Ω; beam = beam, active = false)
````

## Build Sequence

We measure the excited–state population as well as the atomic motion in the
transverse directions, while applying a single resonant pulse.

````julia
add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))
add_detector!(system, MotionDetectorSpec(atom; dims = [1, 2], name = "atom"))

seq = Sequence(dt)
@sequence seq begin
    Pulse(coupling, pulse_duration)
end
````

## Run simulations

````julia
out = play(system, seq; initial_state = g, shots = 100)
````

## Plot results

The plot shows the individual quantum trajectories (red, faint) and their
mean excited–state population as a function of time (black line).

````julia
tlist = out.times
plt = Plots.plot(
    tlist .* 1e6,
    out.detectors["P_e"],
    label     = "",
    xlabel    = "Time (μs)",
    ylabel    = "Population",
    title     = "Rabi oscillations with atomic motion",
    linewidth = 0.2,
    alpha     = 0.1,
    color     = :red,
    legend    = :none,
)

Plots.plot!(
    plt,
    tlist .* 1e6,
    mean(out.detectors["P_e"], dims = 2),
    color     = :black,
    linewidth = 3,
)

plt
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

