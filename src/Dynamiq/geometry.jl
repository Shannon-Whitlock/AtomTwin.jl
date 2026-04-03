#####################################################################################
# Atoms
#####################################################################################

"""
    Atom(x, v, m, λs, γs, ls)

Simple atomic model with position, velocity, and scalar polarizability.

Fields:
- `x::Vector{Float64}`: atomic position in real space.
- `v::Vector{Float64}`: atomic velocity.
- `m::Float64`: atomic mass.
- `λ::NTuple{N,Float64}`: wavelengths of relevant transitions (in meters).
- `γ::NTuple{N,Float64}`: spontaneous emission rates (in Hz).
- `ls::NTuple{N,Int}`: line-strength factors for each transition.
- `alpha::F`: polarizability as a function of probe wavelength.

The constructor builds a dispersion-like polarizability

\\[
\\alpha(\\omega) \\propto \\sum_i \\frac{\\text{ls}_i \\gamma_i}{\\omega - \\omega_i}
\\]

from the transition data, where `ω(λ) = 2πc/λ`.
"""
struct Atom{N,F}
    x::Vector{Float64}
    v::Vector{Float64}
    m::Float64
    λ::NTuple{N,Float64}  # wavelength of relevant transitions in m
    γ::NTuple{N,Float64}  # spontaneous emission rates Hz
    ls::NTuple{N,Int}     # linestrength factors
    alpha::F              # polarizability as a function of wavelength

    function Atom(x::Vector{Float64}, v::Vector{Float64}, m::Float64, λs, γs, ls)
        N = length(λs)

        w(x) = 2 * pi * c / x
        d(λ, x) = w(x) - w(λ)

        alpha(x) = pi / 2 * c^2 / w(x)^3 * sum(ls .* γs ./ d.(λs, x))

        new{N,typeof(alpha)}(x, v, m, Tuple(λs), Tuple(γs), Tuple(ls), alpha)
    end
end

#####################################################################################
# Lasers
#####################################################################################

"""
    BasisPolarization(unit_k, polarization)

Convert a Jones vector `polarization` into the spherical (q-basis) polarization
components for a beam with propagation direction `unit_k`.

- `unit_k`: 3D unit vector along the beam propagation direction.
- `polarization`: transverse Jones vector in the lab frame.

Returns a 3-component complex vector corresponding to \\(\\sigma^+\\), \\(\\pi\\),
and \\(\\sigma^-\\) components in the quantization basis aligned with `unit_k`.
"""
function BasisPolarization(unit_k, polarization)
    new_pol = polarization / norm(polarization)
    theta = acos(unit_k[3])

    if theta != 0
        phi = acos(unit_k[1] / sin(theta))

        ex = Vector{ComplexF64}([
            cos(theta) * exp(im * phi) / sqrt(2),
            -sin(theta),
            -cos(theta) * exp(-im * phi) / sqrt(2),
        ])
        ey = Vector{ComplexF64}([
            im * exp(im * phi) / sqrt(2),
            0,
            im * exp(-im * phi) / sqrt(2),
        ])
    else
        phi = 0

        ex = Vector{ComplexF64}([
            1 / sqrt(2),
            0,
            -1 / sqrt(2),
        ])
        ey = Vector{ComplexF64}([
            im / sqrt(2),
            0,
            im / sqrt(2),
        ])
    end

    polarization_q = new_pol[1] * ex + new_pol[2] * ey
    return polarization_q
end
