"""
Define the abstract interfaces for atoms, levels, fields, and instructions.
"""

"""
    AbstractLevel

Abstract supertype for atomic energy levels (ground, excited, hyperfine, etc.).
All level types must implement this interface.
"""
abstract type AbstractLevel end

"""
    AbstractManifold

Abstract supertype for manifolds of atomic levels (e.g. `FineManifold`,
`HyperfineManifold`). Provides a common interface for iteration and
access to the underlying levels.
"""
abstract type AbstractManifold end

"""
    AbstractField

Abstract supertype for electromagnetic fields (e.g. Rabi couplings, etc.)
that couple atomic levels.
"""
abstract type AbstractField end

"""
    AbstractNoiseModel

Abstract supertype for noise models that can be applied to fields
(e.g., laser phase noise).
"""
abstract type AbstractNoiseModel end

"""
    AbstractInstruction

Abstract supertype for all low-level control instructions (pulses,
switching, waits, tweezer moves, ramps, etc.).
"""
abstract type AbstractInstruction end

