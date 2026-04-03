"""
    PositionModifier <: AbstractModifier

Modifier that overwrites the position of a beam at each time step using
precomputed trajectories.

`PositionModifier` is typically used to impose an externally defined
trajectory on a beam, rather than integrating motion from forces.

# Fields

- `beam::AbstractBeam`: Beam whose position `r0` will be updated.
- `vals::Vector{Vector{Float64}}`: Per-step target positions (full 3D or more).
- `tspan::Vector{Float64}`: Time grid associated with `vals`.
- `dims::Vector{Int}`: Components of `r0` to overwrite at each step.

# Constructor

- `PositionModifier(beam, vals, tspan; dims = [1:3...])`
"""
struct PositionModifier <: AbstractModifier
    beam::AbstractBeam
    vals::Vector{Vector{Float64}}
    tspan::Vector{Float64}
    dims::Vector{Int}

    function PositionModifier(beam::AbstractBeam,
                              vals::Vector{Vector{Float64}},
                              tspan::Vector{Float64};
                              dims = [1:3...])
        new{}(beam, vals, tspan, dims)
    end
end

"""
    update!(m::PositionModifier, i)

Overwrite the beam position components `r0[d]` at time step `i` using
the stored trajectory `m.vals[i]`.

Only the indices listed in `m.dims` are updated.
"""
function update!(m::PositionModifier, i::Int)
    beam = m.beam
    for d in m.dims
        beam.r0[d] = m.vals[i][d]
    end
end
