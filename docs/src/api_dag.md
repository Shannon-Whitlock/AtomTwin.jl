# DAG system

AtomTwin builds system models from a list of *nodes* — each node encodes
the recipe for one physical quantity (a coupling, a detuning, a decay channel, a
beam, …) together with its parametric dependencies. This list is called the
computation graph, or DAG (*directed acyclic graph*), of the system.

## How it works

When you call `add_coupling!`, `add_detuning!`, etc., AtomTwin creates one or
more nodes and appends them to `sys.nodes`. Each node holds the physical
parameters (which may be `Parameter` or `ParametricExpression` objects) and a
reference to the compiled Dynamiq object it produces.

### Node lifecycle

Every node goes through three stages:

| Stage | Function | When called |
|-------|----------|-------------|
| Build | `build_node!(node, basis)` | At `add_*!` time, using **default** parameter values. |
| Compile | `compile_node!(node, basis, rng, param_values)` | Once per `compile` call (per `play`). Resolves overrides and samples noise. |
| Recompile | `recompile_node!(node, obj, rng, param_values)` | Once per Monte Carlo shot on the thread-local job copy. |

Nodes are compiled **in insertion order**. Dependencies (e.g. `BeamNode`) must
be added before the nodes that depend on them (e.g. the `CouplingNode` that
reads the resolved beam).

### Value resolution protocol

Any value stored inside a node - including values stored inside the node's
sub-objects - follows the same two-method protocol:

- `_resolve_node_default(x)`: evaluates `x` using parameter defaults (called at build time)
- `_resolve_node_value(x, param_values, rng)`: evaluates `x` with override lookup and optional noise sampling (called at compile/recompile time)

This protocol is implemented for `Number`, `Parameter`, `ParametricExpression`,
and extended by types such as `GaussianPosition`, `MaxwellBoltzmann`,
`MixedPolarization`, and `BeamRabiFrequency`.

---

## Node types

### BeamNode

Resolves a `ParametricBeam` (or a concrete beam) to a concrete `AbstractBeam`
at compile time. Must be added to `sys.nodes` before any `CouplingNode` that
depends on it. Created automatically by `add_coupling!` when a beam argument is
provided.

---

### CouplingNode

Created by `add_coupling!` for uniform (spatially flat) couplings. Holds a
parametric Rabi frequency `Ω` and maps to a `GlobalCoupling` Hamiltonian term.

---

### NoisyCouplingNode

Created by `add_coupling!` when a `LaserPhaseNoiseModel` is supplied. Wraps
the `GlobalCoupling` in a `NoisyField` so that per-shot phase noise is applied
during the modifier loop in `recompile!`.

---

### PlanarCouplingNode

Created by `add_coupling!` for position-dependent couplings driven by a
`PlanarBeam`. Produces a `PlanarCoupling` Hamiltonian term.

---

### DetuningNode

Created by `add_detuning!`. Holds a parametric detuning `delta` (rad/s) and
maps to a `Detuning` diagonal Hamiltonian term.

---

### DecayNode

Created by `add_decay!` and `add_dephasing!`. Holds a parametric rate `Gamma`
and maps to a `Jump` Lindblad operator.

---

### InteractionNode

Created by `add_interaction!`. Holds a parametric interaction strength `V` and
maps to an `Interaction` two-atom Hamiltonian term.

---

## Value types

Value types implement the resolution protocol and can be used anywhere a
parametric quantity is accepted (e.g. as a field inside a node or as a
constructor argument).

### GaussianPosition

```@docs
GaussianPosition
gaussian
```

### MaxwellBoltzmann

```@docs
MaxwellBoltzmann
maxwellboltzmann
```

### MixedPolarization

```@docs
MixedPolarization
estimate_PER
estimate_PER_dB
```
