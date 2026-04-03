"""
Describe low-level hardware actions as instruction types pulses, waits, on/off switching, and movements.
"""

"""
    Switchable

Union of field-like or coupling objects that can be switched, pulsed,
or otherwise controlled by the instruction layer.
"""
const Switchable = Union{Dynamiq.AbstractField, GlobalCoupling, PlanarCoupling, Detuning, NoisyField}

"""
    Pulse(couplings::Vector{<:Switchable}, duration, ampl, amplitudes)

Apply laser coupling for a fixed duration or with a shaped amplitude profile.

# Fields
- `couplings::Vector{Switchable}` – list of `Switchable` objects (e.g. `GlobalCoupling{NLevelAtom}`) to be modulated.
- `duration::Float64` – pulse duration in seconds. The compiled pulse will span exactly this duration.
- `ampl::ComplexF64` – overall coupling strength (relative to the default). Can be complex. It scales the entire amplitude shape.
- `amplitudes::Vector{ComplexF64}` – vector of complex amplitudes defining the pulse shape. If empty, the pulse has constant amplitude `ampl` over `duration`. If non‑empty, the shape defined by `amplitudes` is resampled and interpolated over the full duration `duration` (using third‑order Lagrange interpolation).
"""
struct Pulse <: AbstractInstruction
    couplings::Vector{Switchable}
    duration::Float64
    ampl::ComplexF64
    amplitudes::Vector{ComplexF64} 
end

"""
    Pulse(couplings::Vector{<:Switchable}, duration; ampl=1.0, amplitudes=ComplexF64[])

Construct a `Pulse` acting on multiple couplings simultaneously, all with the same duration and relative amplitude `ampl`.

# Arguments
- `couplings`: list of `Switchable` objects to be modulated.
- `duration`: pulse duration (seconds). The compiled pulse will span exactly this duration.
- `ampl`: overall coupling strength (relative to default). Scales the entire amplitude shape.
- `amplitudes`: optional vector of amplitudes defining the pulse shape. If empty, the pulse has constant amplitude `ampl` over `duration`. If non‑empty, the shape is resampled and interpolated over the full duration `duration` (using third‑order Lagrange interpolation).

See `Pulse(couplings, duration, ampl, amplitudes)` for details on how the shape is interpolated.
"""
function Pulse(couplings::Vector, duration::Float64; ampl=1.0, amplitudes=ComplexF64[])
    Pulse(couplings, duration, ampl, amplitudes)
end

"""
    Pulse(coupling::Switchable, duration; ampl=1.0, amplitudes=ComplexF64[])

Convenience constructor for a single‑coupling pulse. Equivalent to `Pulse([coupling], duration; ampl, amplitudes)`.

See `Pulse(couplings, duration, ampl, amplitudes)` for details on how the shape is interpolated.
"""
Pulse(coupling::Switchable, duration::Float64; ampl=1.0, amplitudes=ComplexF64[]) = Pulse([coupling], duration, ampl, amplitudes)

"""
    On(couplings)

Turn on one or more couplings.

- `couplings` is a vector of `Switchable` objects to be switched on.
"""
struct On <: AbstractInstruction
    couplings::Vector{Switchable}
end

"""
    On(coupling::Switchable)

Convenience constructor for a single coupling. Equivalent to `On([coupling])`.
"""
On(coupling::Switchable) = On([coupling])

"""
    Off(couplings)

Turn off one or more couplings.

- `couplings` is a vector of `Switchable` objects to be switched off.
"""
struct Off <: AbstractInstruction
    couplings::Vector{Switchable}
end

"""
    Off(coupling::Switchable)

Convenience constructor for a single coupling. Equivalent to `Off([coupling])`.
"""
Off(coupling::Switchable) = Off([coupling])

"""
    Wait(duration)

Idle period with no explicit control actions.

- `duration` is the wait time (seconds).
"""
struct Wait <: AbstractInstruction
    duration::Float64
end

"""
    MoveCol(tweezers, cols, delta, duration; sweep = :min_jerk)

Move one or more columns of tweezers simultaneously.

- `cols` can be a single integer or a vector of integers specifying the columns.
- `delta` is the frequency shift (Hz).
- `duration` is the total time for the move (seconds).
- `sweep` is the sweep profile: built-in symbol (`:linear`, `:cosine`, `:min_jerk`)
  or custom function `s -> f(s)` mapping [0,1] → [0,1].
"""
mutable struct MoveCol <: AbstractInstruction
    tweezers::TweezerArray               
    cols::AbstractVector{Int}
    delta::Union{Float64, Parameter, ParametricExpression} 
    duration::Float64
    sweep::Union{Symbol, Function}
end

"""
    MoveCol(tweezers, col, delta, duration; sweep = :min_jerk)

Convenience constructor for moving a single column. `col` can be an `Int`
or a parametric expression; it is stored in the general `cols` field.
"""
MoveCol(ta::TweezerArray, col, delta, duration; sweep=:min_jerk) = MoveCol(
    ta, col isa Int ? [col] : col, delta, duration, sweep)

"""
    MoveRow(tweezers, rows, delta, duration; sweep = :min_jerk)

Move one or more rows of tweezers simultaneously.

- `rows` can be a single integer or a vector of integers specifying the rows.
- `delta` is the frequency shift (Hz).
- `duration` is the total time for the move (seconds).
- `sweep` is the sweep profile: built-in symbol (`:linear`, `:cosine`, `:min_jerk`)
  or custom function `s -> f(s)` mapping [0,1] → [0,1].
"""
struct MoveRow <: AbstractInstruction
    tweezers::TweezerArray
    rows::AbstractVector{Int}
    delta::Union{Float64, Parameter, ParametricExpression} 
    duration::Float64
    sweep::Union{Symbol, Function}   # Symbol or custom function
end

"""
    MoveRow(tweezers, row, delta, duration; sweep = :min_jerk)

Convenience constructor for moving a single row. Wraps `row` into a vector.
"""
MoveRow(ta::TweezerArray, row::Int, delta::Float64, duration::Float64; sweep=:min_jerk) = MoveRow(ta, [row], delta, duration, sweep)

"""
    RampRow(tweezers, rows, final_amplitude, ramp_time)

Ramp the amplitude of one or more rows of tweezers simultaneously.

- `rows` can be a single integer or a vector of integers specifying the rows.
- `final_amplitude` is the value to ramp each row to (scalar or vector).
- `ramp_time` is the total duration of the ramp (seconds).
"""
struct RampRow <: AbstractInstruction
    tweezers::TweezerArray
    rows::AbstractVector{Int}
    final_amplitude::Union{Float64,AbstractVector{Float64}}
    ramp_time::Float64
end

"""
    RampRow(tweezers, row, final_amplitude, ramp_time)

Convenience constructor for ramping a single row. Wraps `row` into a vector.
"""
RampRow(ta::TweezerArray, row::Int, final_amplitude::Float64, ramp_time::Float64) =
    RampRow(ta, [row], final_amplitude, ramp_time)

"""
    RampCol(tweezers, cols, final_amplitude, ramp_time)

Ramp the amplitude of one or more columns of tweezers simultaneously.

- `cols` can be a single integer or a vector of integers specifying the columns.
- `final_amplitude` is the value to ramp each column to (scalar or vector).
- `ramp_time` is the total duration of the ramp (seconds).
"""
struct RampCol <: AbstractInstruction
    tweezers::TweezerArray
    cols::AbstractVector{Int}
    final_amplitude::Union{Float64,AbstractVector{Float64}}
    ramp_time::Float64
end

"""
    RampCol(tweezers, col, final_amplitude, ramp_time)

Convenience constructor for ramping a single column. Wraps `col` into a vector.
"""
RampCol(ta::TweezerArray, col::Int, final_amplitude::Float64, ramp_time::Float64) =
    RampCol(ta, [col], final_amplitude, ramp_time)

"""
    AmplCol(tweezers, col, ampl)

Set the relative amplitude of a single tweezer column.

- `col` is the column index.
- `ampl` is the new relative amplitude.
"""
struct AmplCol <: AbstractInstruction
    tweezers::TweezerArray
    col::Int
    ampl::Float64
end

"""
    AmplRow(tweezers, row, ampl)

Set the relative amplitude of a single tweezer row.

- `row` is the row index.
- `ampl` is the new relative amplitude.
"""
struct AmplRow <: AbstractInstruction
    tweezers::TweezerArray
    row::Int
    ampl::Float64
end

"""
    FreqCol(tweezers, col, freq)

Set the frequency of a single tweezer column.

- `col` is the column index.
- `freq` is the absolute frequency (Hz).
"""
struct FreqCol <: AbstractInstruction
    tweezers::TweezerArray
    col::Int
    freq::Float64
end

"""
    FreqRow(tweezers, row, freq)

Set the frequency of a single tweezer row.

- `row` is the row index.
- `freq` is the absolute frequency (Hz).
"""
struct FreqRow <: AbstractInstruction
    tweezers::TweezerArray
    row::Int
    freq::Float64
end

"""
    Parallel(parts::Vector{<:AbstractInstruction})

Group several instructions into a single instruction that is executed
in parallel.

Each sub-instruction is compiled with the same time interval. The
resulting `Parallel` segment:

- Uses a common time span from the earliest start to the latest
  end time among all sub-instructions.
- Concatenates all emitted modifiers from the sub-instructions.
- Relies on each modifier's `update!` method being a no-op once the
  modifier has exhausted its internal data, so shorter instructions
  naturally hold their final value while longer ones continue.

This effectively "pads" shorter instructions in time up to the longest
instruction without modifying the modifiers themselves.
"""
struct Parallel <: AbstractInstruction
    parts::Vector{AbstractInstruction}
end
