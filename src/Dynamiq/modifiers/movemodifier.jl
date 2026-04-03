"""
    MoveModifier <: AbstractModifier

Modifier that updates a beam position incrementally over time to realize a
smooth displacement according to a user-defined schedule.

Unlike `PositionModifier`, which overwrites positions, `MoveModifier` stores
per-step increments that are added to the beam position.

# Fields

- `beam::AbstractBeam`: Beam whose position `r0` will be moved.
- `vals::Vector{Vector{Float64}}`: Per-step displacement increments.
- `tspan::Vector{Float64}`: Internal time grid (left edges of increments).
- `dims::Vector{Int}`: Components of `r0` to move at each step.

# Constructor

    MoveModifier(beam, displacement, tspan;
                 dims = [1:3...],
                 schedule = s -> s)

- `displacement::Vector{Float64}`: Total displacement to be applied over `tspan`.
- `tspan::Vector{Float64}`: Time grid from start to end (length ≥ 2).
- `schedule::Function`: Mapping `s ∈ [0, 1]` ↦ scalar factor that shapes the
  trajectory; `s` is normalized time `(t - t0) / (t_end - t0)`.
"""
struct MoveModifier <: AbstractModifier
    beam::AbstractBeam
    vals::Vector{Vector{Float64}}
    tspan::Vector{Float64}
    dims::Vector{Int}

    function MoveModifier(beam::AbstractBeam,
                          displacement::Vector{Float64},
                          tspan::Vector{Float64};
                          dims = [1:3...],
                          schedule = s -> s)

        # Guard against empty or single-point tspan
        length(tspan) ≥ 2 ||
            error("MoveModifier: tspan must contain at least two time points")

        # Total duration (might be small)
        T = tspan[end] - tspan[1]
        @assert T>0 "Move interval cannot be zero"

        # Build route with normalized time s ∈ [0,1]
        # s_k = (t_k - t_0) / T, so s_1 = 0, s_end = 1
        route = Vector{Vector{Float64}}(undef, length(tspan))
        t0    = tspan[1]
        for (i, t) in enumerate(tspan)
            s = (t - t0) / T
            route[i] = schedule(s) * displacement
        end

        # Increments between successive route points
        vals = diff(route)  # length = length(tspan) - 1

        # Use the left edges of each increment as internal tspan
        internal_tspan = tspan[1:end-1]

        return new(beam, vals, internal_tspan, dims)
    end
end

"""
    update!(m::MoveModifier, i)

Increment the beam position `r0` at time step `i` by the stored displacement
`m.vals[i]` on the components listed in `m.dims`, if `i` is within bounds.
"""
function update!(m::MoveModifier, i::Int)
    if i <= length(m.vals)
        for d in m.dims
            m.beam.r0[d] += m.vals[i][d]
        end
    end
end
