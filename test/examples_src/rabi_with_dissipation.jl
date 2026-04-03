# # Rabi oscillations with dissipation
#
# In this example we simulate driven Rabi oscillations in a two-level atom,
# including spontaneous decay and dephasing, and monitor the excited-state population.
#
# We compare two simulation methods:
#
# - Master equation (Lindblad) evolution of the density matrix.
# - Quantum trajectories (Monte Carlo wavefunction method) with many shots.

# Internal notes for test scripts (not included in docs/examples):         #src
# - Files in `test/examples_src/` are run as tests and also used to        #src
#   generate docs and runnable examples.                                   #src
# - Lines containing `#src` are removed by `make.jl` when generating       #src
#   docs/examples, but are present when running tests.                     #src
# - To omit a line from tests but keep it in docs/examples, wrap it in     #src
#   an `if false ... end` block that is itself tagged with `#src`:         #src
#                                                                          #src
#       if false            #src                                           #src
#           using Plots     # omitted in tests, kept in docs/examples      #src
#       end                 #src                                           #src

using AtomTwin
using StatsBase
if false            #src
using Plots      
end                 #src

# ## Parameters

Ω     = 2π * 0.5e6          # Rabi frequency (rad/s)
Gamma = 2π * 50e3           # Spontaneous emission rate |e⟩ → |g⟩ (rad/s)
gamma = 2π * 50e3           # Pure dephasing rate of the excited state (rad/s)

pulse_duration = 10e-6      # Total pulse duration (s)
dt             = 5e-9       # Time step (s)
shots          = 400        # Number of quantum trajectories

descriptor = "Rabi with dissipation: Ω/2π = $(Ω/2π/1e6) MHz, Gamma = $(Gamma/2π/1e3) kHz, gamma = $(gamma/2π/1e3) kHz" #src

# ## System definition

# Define a simple two-level atom
g, e = Level("g"), Level("e")
atom = Atom(; levels = [g, e])
system = System(atom) 
 
# Add a coherent resonant drive between |g⟩ and |e⟩ with Rabi frequency Ω,
# as well as decay and dephasing of the excited state (rates Gamma and gamma respectively)
coupling = add_coupling!(system, atom, g => e, Ω; active = false) 
decay = add_decay!(system, atom, e => g, Gamma; active = true) 
deph = add_dephasing!(system, atom, e, gamma; active = true) 


# ## Build sequence
# Register population detector on |e⟩
add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e")) 

seq = Sequence(dt)
@sequence seq begin
    Pulse(coupling, pulse_duration)
end 

# ## Run simulations
runtime = @elapsed begin                                    #src
out_me = play(system, seq; initial_state = g, density_matrix = true) # master equation
out_qt = play(system, seq; initial_state = g, shots = 400) # quantum trajectories
end     
checksum_data  = out_me.detectors["P_e"]                    #src
                                                            #src
## Validate physical correctness                           #src
# 1. Population is bounded in [0, 1] at all times         #src
@assert all(x -> 0.0 - 1e-6 ≤ x ≤ 1.0 + 1e-6, out_me.detectors["P_e"]) "Master equation population out of [0,1]" #src
# 2. Dissipation damps the amplitude over time:           #src
#    excited-state population at t=end should be less     #src
#    than at the first Rabi peak (around t=0.5/Ω ≈ 1 μs) #src
pe_me = out_me.detectors["P_e"]                           #src
peak_idx = argmax(pe_me)                                  #src
@assert pe_me[end] < pe_me[peak_idx] "Dissipation should reduce late-time peak population" #src
# 3. QT mean should track master equation result (mean absolute deviation < 5%) #src
pe_qt_mean = vec(mean(out_qt.detectors["P_e"], dims = 2))                        #src
@assert mean(abs.(pe_qt_mean .- pe_me)) < 0.05 "QT mean deviates from master equation by more than 5% on average" #src

# ## Plot results

# Collect data
if false  #src
tlist = out_qt.times                 # Time grid (s)
t_us  = tlist .* 1e6                 # Convert to microseconds

P_traj = out_qt.detectors["P_e"]     # excited state populations for each trajectory
P_traj_mean = mean(P_traj, dims = 2) # average excited state population
P_me = out_me.detectors["P_e"]       # average from master equation simulation

# Create a plot showing:
#   - Individual trajectories (red, faint)
#   - Trajectory average (black)
#   - Master-equation result (blue)

plt = Plots.plot(
    t_us, P_traj;
    label      = "",
    xlabel     = "Time (μs)",
    ylabel     = "Excited-state population",
    title      = "Rabi oscillations with dissipation",
    linewidth  = 0.2,
    alpha      = 0.2,
    color      = :red,
    legend     = :bottomleft,
)

Plots.plot!(plt, t_us, P_traj_mean;
      label = "Quantum trajectories (mean)",
      color = :black,
      linewidth = 3)

Plots.plot!(plt, t_us, P_me;
      label = "Master equation",
      color = :blue,
      linewidth = 3)

plt
end  #src
