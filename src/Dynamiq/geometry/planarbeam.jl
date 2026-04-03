"""
    PlanarBeam <: AbstractBeam

Plane-wave beam with fixed intensity, propagation direction, and polarization.

Fields:
- `λ::Float64`: Wavelength in meters.
- `I::Float64`: Intensity in \\(\\mathrm{W/m^2}\\).
- `E_field::Float64`: Electric-field amplitude in natural units (scaled by \\(\\hbar\\)).
- `k::NTuple{3,Float64}`: Wavevector components in \\(\\mathrm{m^{-1}}\\).
- `unit_k::NTuple{3,Float64}`: Normalized propagation direction.
- `polarization::NTuple{3,ComplexF64}`: Polarization in the spherical q-basis
  (q = -1, 0, 1).
- `_coeff::Base.RefValue{ComplexF64}`: Complex amplitude envelope used by
  time-dependent modifiers.

The constructor takes a Cartesian propagation `direction` and a Jones
polarization vector and converts them into the q-basis via `BasisPolarization`.
"""
struct PlanarBeam <: AbstractBeam
    λ::Float64                # Wavelength [m]
    I::Float64                # Intensity [W/m^2]
    E_field::Float64          # Amplitude of Electric Field [N/c] (scaled)
    k::NTuple{3,Float64}      # k vector [1/m]
    unit_k::NTuple{3,Float64} # Normalized k vector
    polarization::NTuple{3,ComplexF64} # Polarization in q-basis (q = -1, 0, 1)
    _coeff::Base.RefValue{ComplexF64}

    function PlanarBeam(λ, I, direction, polarization)
        e_field       = sqrt(2 * I / (c * ε0)) / hbar
        unit_k        = normalize(direction)
        q_polarization = BasisPolarization(unit_k, polarization)
        k             = (2 * pi / λ) * unit_k
        new(λ, I, e_field, Tuple(k), Tuple(unit_k), Tuple(q_polarization),
            Ref(Complex(1.0)))
    end
end

"""
    getbeams(pb::PlanarBeam) -> Vector{AbstractBeam}

Return a vector containing the beam pb
"""
getbeams(pb::PlanarBeam) = AbstractBeam[pb]