"""
    PhotoDetector{A} <: AbstractDetector{A}

Simple photon-counting detector that registers integer counts per time step.

This detector is typically used together with stochastic jump processes;
each time a relevant event occurs, `write!(d, i)` is called to increment
the count at time index `i`.

# Fields

- `atom::A`: Associated atom or source object.
- `vals::Vector{Int}`: Per-timestep integer count of detected events.
- `tspan::Vector{Float64}`: Time vector corresponding to the bins in `vals`.

# Constructor

- `PhotoDetector(atom, tspan)`
"""
struct PhotoDetector{A} <: AbstractDetector
    atom::A
    vals::Vector{Int}
    tspan::Vector{Float64}

    function PhotoDetector(atom::AbstractAtom,
                           tspan::Vector{Float64})
        new{typeof(atom)}(atom, zeros(Int, length(tspan)), tspan)
    end
end

"""
    write!(d::PhotoDetector, i)

Register a single detection event at time step `i` by incrementing the
count in `d.vals[i]`.
"""
function write!(d::PhotoDetector{A}, i::Int) where {A}
    d.vals[i] += 1
end

"""
    reset!(d::PhotoDetector)

Reset all recorded counts to zero while preserving the length and `tspan`.
"""
function reset!(d::PhotoDetector)
    fill!(d.vals, zero(eltype(d.vals)))
end
