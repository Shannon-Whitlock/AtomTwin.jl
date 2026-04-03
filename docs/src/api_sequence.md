# Sequences and Simulations

AtomTwin allows to specify simulations in a very similar way to how a hardware control system operates. The main object here is the `Sequence` which contains a set of instructions that define a simulation

## Building Sequences
```@autodocs
Modules = [AtomTwin]
Pages = ["sequence.jl"]
Public = true
```

## Instructions
AtomTwin provides a small set of instruction types that describe low-level control operations on couplings and tweezer arrays (pulses, on/off switching, motion, ramps, etc.).

### Pulse and switching instructions
```@docs
AbstractInstruction
Pulse
On
Off
Wait
```

### Tweezer instructions
```@docs
MoveCol
MoveRow
RampCol
AmplCol
AmplRow
FreqCol
FreqRow
```

### Example: simple pulse sequence
```@example

using AtomTwin

Ω = 2*pi*1e6

g, e = Level("g"), Level("e")
atom = Atom(; levels = [g, e])
system = System(atom)

coupling = add_coupling!(
    system, atom, g => e, Ω;
    active = false)

seq = Sequence(1e-9)           # 1 ns timestep
@sequence seq begin
    Pulse(coupling, 1e-6)      # 1 μs pulse
    Wait(0.5e-6)               # 0.5 μs wait
end
```

## Running Simulations
```@docs
play
compile
recompile!
AtomTwin.Dynamiq.evolve!
```
