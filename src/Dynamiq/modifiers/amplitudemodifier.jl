"""
    AmplitudeModifier{F} <: AbstractModifier

Generic modifier that updates the complex amplitude `_coeff[]` of a field
or beam as a function of time.

# Fields

- `field::F`: Target field or beam (must have a `_coeff::Ref{ComplexF64}`).
- `vals::Vector{ComplexF64}`: Time series of complex coefficients aligned with simulation time grid.

# Constructors

- `AmplitudeModifier(field::AbstractField, vals)`
- `AmplitudeModifier(beam::AbstractBeam, vals)`
"""
struct AmplitudeModifier{F} <: AbstractModifier
    field::F
    vals::Vector{ComplexF64}

    function AmplitudeModifier(field::Union{AbstractField, AbstractBeam},
                               vals::AbstractVector{<:Number})
        # Convert once - no allocation if already ComplexF64
        new{typeof(field)}(field, convert(Vector{ComplexF64}, vals))
    end
end


"""
    update!(m::AmplitudeModifier, i)

Set the complex amplitude `_coeff[]` of the underlying field or beam
to `m.vals[i]` at time step `i`, if `i` is within bounds.
"""
@inline function update!(m::AmplitudeModifier{F}, i::Int) where {F}
    @inbounds if i <= length(m.vals)
        m.field._coeff[] = m.vals[i] 
    end
end
