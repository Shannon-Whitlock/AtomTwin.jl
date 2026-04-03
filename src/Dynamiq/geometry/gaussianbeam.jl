using .Units

#------------------------------------------------------------------------------
# Axis-aligned Gaussian beam
#------------------------------------------------------------------------------

"""
    GaussianBeam <: AbstractBeam

Axis-aligned (z-propagating) Gaussian beam with a 3D harmonic envelope.

Fields:
- `λ::Float64`: Wavelength in meters.
- `w0::Float64`: Transverse beam waist (radius) in meters.
- `P::Float64`: Optical power in watts.
- `I0::Float64`: Peak intensity at the waist.
- `w0z::Float64`: Effective axial waist (harmonic approximation) in meters.
- `r0::Vector{Float64}`: Beam center position \\((x_0, y_0, z_0)\\).
- `_coeff::Base.RefValue{ComplexF64}`: Complex amplitude envelope used by
  time-dependent modifiers.

The axial waist `w0z` is chosen such that the simple 3D Gaussian intensity
\\(I \\propto \\exp[-2(x^2/w_0^2 + y^2/w_0^2 + z^2/w_{0z}^2)]\\) matches the
quadratic expansion of a paraxial Gaussian near the focus.
"""
mutable struct GaussianBeam <: AbstractBeam
    λ::Float64      # wavelength in m
    w0::Float64     # beam waist in m
    P::Float64      # power in watts
    I0::Float64     # peak intensity
    w0z::Float64    # axial waist (harmonic approx.)
    r0::Vector{Float64}
    _coeff::Base.RefValue{ComplexF64}
end

"""
    GaussianBeam(λ, w0, P; r0 = [0, 0, 0])

Construct a `GaussianBeam` from wavelength, transverse waist, and power.

- `λ::Float64`: Wavelength (m).
- `w0::Float64`: Waist radius (m).
- `P::Float64`: Optical power (W).
- `r0::Vector{Float64}`: Beam center position (default `[0, 0, 0]`).

`I0` and `w0z` are derived automatically from these parameters.
"""
function GaussianBeam(λ::Float64, w0::Float64, P::Float64; r0 = [0.0, 0.0, 0.0])
    I0  = 2 * P / (π * w0^2)
    w0z = sqrt(2) * π * w0^2 / λ  # effective axial "waist" from quadratic expansion
    return GaussianBeam(λ, w0, P, I0, w0z, r0, Ref(ComplexF64(1.0)))
end

"""
    GaussianBeam(; λ, w0, P, r0 = [0, 0, 0])

Keyword-only constructor for `GaussianBeam`. Parameters are the same
as the positional constructor.
"""
function GaussianBeam(; λ::Float64,
                      w0::Float64,
                      P::Float64,
                      r0::Vector{Float64} = [0.0, 0.0, 0.0])
    I0  = 2 * P / (π * w0^2)
    w0z = sqrt(2) * π * w0^2 / λ
    return GaussianBeam(λ, w0, P, I0, w0z, r0, Ref(ComplexF64(1.0)))
end

"""
    copy(b::GaussianBeam)

Create a deep copy of a `GaussianBeam`, including a new `r0` vector and
a new `_coeff` reference with the same complex value.
"""
function copy(b::GaussianBeam)
    GaussianBeam(
        b.λ,
        b.w0,
        b.P,
        b.I0,
        b.w0z,
        copy(b.r0),
        Ref(b._coeff[]),
    )
end

"""
    intensity(b::GaussianBeam, r)

Intensity of the axis-aligned Gaussian beam at position `r`.

Implements a 3D harmonic Gaussian profile

\\[
I(r) = I_0 \\exp\\bigl[-2(x^2/w_0^2 + y^2/w_0^2 + z^2/w_{0z}^2)\\bigr]
\\, |c|^2,
\\]

with an 8-waist cutoff in the approximate ellipsoid
\\(x^2/w_0^2 + y^2/w_0^2 + z^2/w_{0z}^2 \\le 16\\).
"""
function intensity(b::GaussianBeam, r)
    dx = r[1] - b.r0[1]
    dy = r[2] - b.r0[2]
    dz = r[3] - b.r0[3]

    w0  = b.w0
    wz0 = b.w0z

    # 8-waist cutoff in 3D harmonic approximation
    if dx^2 / w0^2 + dy^2 / w0^2 + dz^2 / wz0^2 > 16
        return 0.0
    else
        @fastmath m = b.I0 * exp(-2 * (dx^2 / w0^2 + dy^2 / w0^2 + dz^2 / wz0^2)) *
                      abs2(b._coeff[])
        return m
    end
end

"""
    dIdx(b::GaussianBeam, r) -> (dI_dx, dI_dy, dI_dz)

Gradient of the intensity of the axis-aligned Gaussian beam at position `r`.

Returns \\((\\partial I/\\partial x, \\partial I/\\partial y, \\partial I/\\partial z)\\)
for the same 3D harmonic profile used in `intensity`.
"""
function dIdx(b::GaussianBeam, r)
    dx = r[1] - b.r0[1]
    dy = r[2] - b.r0[2]
    dz = r[3] - b.r0[3]

    w0  = b.w0
    wz0 = b.w0z

    if dx^2 / w0^2 + dy^2 / w0^2 + dz^2 / wz0^2 > 16
        return 0.0, 0.0, 0.0
    else
        @fastmath I = b.I0 * exp(-2 * (dx^2 / w0^2 + dy^2 / w0^2 + dz^2 / wz0^2)) *
                      abs2(b._coeff[])
        @fastmath dIdx = -4 * dx / w0^2 * I
        @fastmath dIdy = -4 * dy / w0^2 * I
        @fastmath dIdz = -4 * dz / wz0^2 * I
        return dIdx, dIdy, dIdz
    end
end

"""
    Efield(b::GaussianBeam, r)

Complex electric-field amplitude of an axis-aligned (z-propagating)
Gaussian beam at position `r`.

The envelope is a 3D Gaussian,

\\[
E(r) = \\sqrt{I_0} \\exp\\bigl[-(x^2 + y^2)/w_0^2 - z^2/w_{0z}^2\\bigr]
      \\exp(i k z)\\, c,
\\]

with \\(k = 2\\pi/\\lambda\\) and complex coefficient `c = b._coeff[]`.
An 8-waist cutoff is applied in the same ellipsoid used in `intensity`.
"""
@inline function Efield(b::GaussianBeam, r::Vector{Float64})
    dx = r[1] - b.r0[1]
    dy = r[2] - b.r0[2]
    dz = r[3] - b.r0[3]

    w0  = b.w0
    wz0 = b.w0z

    # 8-waist cutoff in the approximate ellipsoid
    if dx^2 / w0^2 + dy^2 / w0^2 + dz^2 / wz0^2 > 16
        return 0.0 + 0.0im
    else
        k   = 2π / b.λ
        env = sqrt(2/(c*ε0)) * sqrt(b.I0) * exp(-(dx^2 + dy^2) / w0^2 - dz^2 / wz0^2)
        return env * cis(k * dz) * b._coeff[]
    end
end

"""
    efield_scalar(b::GaussianBeam, r) -> ComplexF64

Scalar complex field amplitude at position `r`. Identical to `Efield` for
`GaussianBeam` (which already returns a scalar).
"""
@inline efield_scalar(b::GaussianBeam, r::Vector{Float64}) = Efield(b, r)

"""
    dEdr(b::GaussianBeam, r) -> (dE_dx, dE_dy, dE_dz)

Gradient \\(\\nabla E\\) of the complex electric field amplitude of an
axis-aligned Gaussian beam at position `r`.

Returns a tuple of three complex numbers
\\((\\partial E/\\partial x, \\partial E/\\partial y, \\partial E/\\partial z)\\)
for the same 3D Gaussian model used in `Efield`.
"""
@inline function dEdr(b::GaussianBeam, r::Vector{Float64})
    dx = r[1] - b.r0[1]
    dy = r[2] - b.r0[2]
    dz = r[3] - b.r0[3]

    w0  = b.w0
    wz0 = b.w0z

    if dx^2 / w0^2 + dy^2 / w0^2 + dz^2 / wz0^2 > 16
        return (0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)
    else
        k      = 2π / b.λ
        env    = sqrt(2/(c*ε0)) * sqrt(b.I0) * exp(-(dx^2 + dy^2) / w0^2 - dz^2 / wz0^2)
        ephase = cis(k * dz)

        coeff = b._coeff[]
        E     = env * ephase * coeff

        # ∂env/∂x, ∂env/∂y, ∂env/∂z
        dEnv_dx = -2 * dx / w0^2  * env
        dEnv_dy = -2 * dy / w0^2  * env
        dEnv_dz = -2 * dz / wz0^2 * env

        # Full derivatives: envelope + phase along z
        dEdx = dEnv_dx * ephase * coeff
        dEdy = dEnv_dy * ephase * coeff
        dEdz = (dEnv_dz * ephase + env * (im * k * ephase)) * coeff

        return (dEdx, dEdy, dEdz)
    end
end

"""
    update!(gb::GaussianBeam, i)

Update the internal coefficient `_coeff` of a `GaussianBeam` via any attached
modifiers at time step `i`.

This assumes that external code has associated modifiers with the beam and
defined `update!(beam, modifier, i)` methods.
"""
function update!(gb::GaussianBeam, i)
    for m in gb.modifiers
        update!(gb, m, i)
    end
end

"""
    getwavelength(gb::GaussianBeam) -> Float64

Return the beam wavelength in meters.
"""
getwavelength(gb::GaussianBeam) = gb.λ

"""
    getposition(gb::GaussianBeam) -> Vector{Float64}

Return the beam center position vector `gb.r0`.
"""
getposition(gb::GaussianBeam) = gb.r0

"""
    getbeams(gb::GaussianBeam) -> Vector{AbstractBeam}

Return a vector containing the beam gb
"""
getbeams(gb::GaussianBeam) = AbstractBeam[gb]

#------------------------------------------------------------------------------
# General oriented elliptical Gaussian beam
#------------------------------------------------------------------------------

"""
    GeneralGaussianBeam <: AbstractBeam

Elliptical paraxial Gaussian beam with arbitrary propagation direction.

Fields:
- `λ::Float64`: Wavelength in meters.
- `w0x::Float64`: Waist along local x (m).
- `w0y::Float64`: Waist along local y (m).
- `P::Float64`: Optical power (W).
- `I0::Float64`: Peak intensity.
- `r0::Vector{Float64}`: Waist center in global coordinates.
- `k::Vector{Float64}`: Propagation direction (normalized).
- `u::Vector{Float64}`: Local x axis (perpendicular to `k`).
- `v::Vector{Float64}`: Local y axis (perpendicular to `k` and `u`).
- `pol::Vector{complexF64}`: complex polarization vector in global coordinates
- `_coeff::Base.RefValue{ComplexF64}`: Complex amplitude envelope.

The beam profile is elliptical in the transverse plane spanned by `u` and `v`.
"""
struct GeneralGaussianBeam <: AbstractBeam
    λ::Float64
    w0x::Float64
    w0y::Float64
    P::Float64
    I0::Float64
    r0::Vector{Float64}
    k::Vector{Float64}                 # unit propagation direction
    u::Vector{Float64}
    v::Vector{Float64}
    pol::Vector{ComplexF64}            # complex polarization (unit, ⟂ k)
    _coeff::Base.RefValue{ComplexF64}
end


"""
    GeneralGaussianBeam(λ, w0x, w0y, P, k, pol; r0 = [0,0,0])

Positional constructor for `GeneralGaussianBeam`.

- `λ::Float64`: Wavelength (m).
- `w0x::Float64`, `w0y::Float64`: Transverse waists along local x and y.
- `P::Float64`: Optical power (W).
- `k::Vector{Float64}`: Initial propagation direction (normalized internally).
- `pol::Vector{ComplexF64}`: Polarization vector in global coordinates
- `r0::Vector{Float64}`: Waist center in global coordinates (default [0,0,0])

The local axes `u` and `v` are constructed orthonormal to `k`.
"""
function GeneralGaussianBeam(λ::Float64,
                             w0x::Float64,
                             w0y::Float64,
                             P::Float64,
                             k::Vector{<:Real},
                             pol::Vector{<:Number};
                             r0::Vector{<:Real} = [0.0, 0.0, 0.0]
                             )

    # normalize k
    k_norm = sqrt(sum(x -> x^2, k))
    k_norm == 0 && error("k vector must be nonzero")
    k̂ = [x / k_norm for x in k]

    # polarization: ensure nonzero and transverse
    p = ComplexF64.(pol)
    p_norm = sqrt(sum(abs2, p))
    p_norm == 0 && error("polarization vector must be nonzero")
    p ./= p_norm

    # project out any small component along k̂
    pk = sum(p[i] * k̂[i] for i in 1:3)
    p .-= pk .* ComplexF64.(k̂)

    p_norm2 = sqrt(sum(abs2, p))
    p_norm2 == 0 && error("polarization must be orthogonal to k vector")
    p ./= p_norm2

    # existing u,v construction
    ref = abs(k̂[1]) < 0.99 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    u = cross(k̂, ref)
    u_norm = sqrt(sum(x -> x^2, u))
    u = [x / u_norm for x in u]
    v = cross(k̂, u)

    I0 = 2 * P / (π * w0x * w0y)

    return GeneralGaussianBeam(λ, w0x, w0y, P, I0,
                               copy(r0), k̂, u, v, p, Ref(ComplexF64(1.0)))
end


"""
    GeneralGaussianBeam(; λ, w0x, w0y, P, r0 = [0,0,0], k = [0,0,1])

Keyword-only constructor for `GeneralGaussianBeam`. Parameters are the same
as the positional constructor.
"""
function GeneralGaussianBeam(; λ::Float64,
                             w0x::Float64,
                             w0y::Float64,
                             P::Float64,
                             r0::Vector{Float64} = [0.0, 0.0, 0.0],
                             k::Vector{Float64} = [0.0, 0.0, 1.0],
                             pol::Vector{ComplexF64} = ComplexF64[1.0+0im, 0.0+0im, 0.0+0im])
    GeneralGaussianBeam(λ, w0x, w0y, P; r0 = r0, k = k, pol = pol)
end


"""
    copy(b::GeneralGaussianBeam)

Create a deep copy of a `GeneralGaussianBeam`, including all direction
and position vectors and the `_coeff` reference.
"""
function copy(b::GeneralGaussianBeam)
    GeneralGaussianBeam(
        b.λ,
        b.w0x,
        b.w0y,
        b.P,
        b.I0,
        copy(b.r0),
        copy(b.k),
        copy(b.u),
        copy(b.v),
        copy(b.pol),           
        Ref(b._coeff[]),
    )
end


"""
    local_coords(r, b) -> (x′, y′)

Convert a global position `r` into local transverse coordinates \\((x', y')\\)
in the frame defined by a `GeneralGaussianBeam`.

The local z coordinate \\(z'\\) is the projection along `b.k`.
"""
@inline function local_coords(r::Vector{Float64}, b::GeneralGaussianBeam)
    Δ = [r[i] - b.r0[i] for i in 1:3]
    x′ = sum(Δ[i] * b.u[i] for i in 1:3)
    y′ = sum(Δ[i] * b.v[i] for i in 1:3)
    return x′, y′
end

"""
    intensity(b::GeneralGaussianBeam, r)

Intensity of the elliptical Gaussian beam at position `r`.

In local coordinates \\((x', y')\\),

\\[
I(r) = I_0 \\exp\\bigl[-2(x'^2 / w_{0x}^2 + y'^2 / w_{0y}^2)\\bigr] \\,
|c|^2,
\\]

with an 8-waist cutoff ellipse
\\(x'^2 / w_{0x}^2 + y'^2 / w_{0y}^2 \\le 16\\).
"""
@inline function intensity(b::GeneralGaussianBeam, r::Vector{Float64})
    x′, y′ = local_coords(r, b)

    # 8w0 cutoff ellipse
    if x′^2 / b.w0x^2 + y′^2 / b.w0y^2 > 16
        return 0.0
    else
        return b.I0 * exp(-2 * (x′^2 / b.w0x^2 + y′^2 / b.w0y^2)) *
               abs2(b._coeff[])
    end
end

"""
    dIdx(b::GeneralGaussianBeam, r) -> (dI_dx, dI_dy, dI_dz)

Gradient of the intensity of an elliptical Gaussian beam at global position `r`.

The computation is performed in local coordinates \\((x', y')\\) and then
mapped back to global coordinates via the local axes `u` and `v`.
Returns a tuple of three real numbers.
"""
@inline function dIdx(b::GeneralGaussianBeam, r::Vector{Float64})
    x′, y′ = local_coords(r, b)

    if x′^2 / b.w0x^2 + y′^2 / b.w0y^2 > 16
        return (0.0, 0.0, 0.0)
    else
        m = -4 * b.I0 * exp(-2 * (x′^2 / b.w0x^2 + y′^2 / b.w0y^2)) *
            abs2(b._coeff[])

        # Gradient in global coordinates
        grad = [m * (x′ / b.w0x^2) * b.u[i] + m * (y′ / b.w0y^2) * b.v[i]
                for i in 1:3]
        return (grad[1], grad[2], grad[3])
    end
end

"""
    update!(gb::GeneralGaussianBeam, i)

Update the internal coefficient `_coeff` of a `GeneralGaussianBeam` via any
attached modifiers at time step `i`.
"""
function update!(gb::GeneralGaussianBeam, i)
    for m in gb.modifiers
        update!(gb, m, i)
    end
end

"""
    getwavelength(gb::GeneralGaussianBeam) -> Float64

Return the beam wavelength in meters.
"""
getwavelength(gb::GeneralGaussianBeam) = gb.λ

"""
    getposition(gb::GeneralGaussianBeam) -> Vector{Float64}

Return the beam waist center `gb.r0` in global coordinates.
"""
getposition(gb::GeneralGaussianBeam) = gb.r0

"""
    getbeams(ggb::GeneralGaussianBeam) -> Vector{AbstractBeam}

Return a vector containing the beam ggb
"""
getbeams(ggb::GeneralGaussianBeam) = AbstractBeam[ggb]

"""
    Efield(b::GeneralGaussianBeam, r)

Complex electric-field amplitude of an elliptical paraxial Gaussian beam
at position `r`.

This function accounts for:

- z'-dependent waists `wx(z')`, `wy(z')`.
- Curvatures `Rx(z')`, `Ry(z')`.
- Gouy phases `ζx(z')`, `ζy(z')`.
- Arbitrary beam orientation via local coordinates.
- Elliptical transverse profile.

Returns zero if the point lies outside the 8-waist cutoff ellipse.

# Arguments

- `b::GeneralGaussianBeam`: Beam object.
- `r::Vector{Float64}`: 3D global position.

# Returns

- `E::ComplexF64`: Complex electric field at `r`.
"""
@inline function Efield(b::GeneralGaussianBeam, r::Vector{Float64})
    x′, y′ = local_coords(r, b)
    z′ = sum((r[i] - b.r0[i]) * b.k[i] for i in 1:3)

    λ = b.λ
    w0x, w0y = b.w0x, b.w0y
    k = 2π / λ
    zRx = π * w0x^2 / λ
    zRy = π * w0y^2 / λ

    wx2 = w0x^2 * (1 + (z′ / zRx)^2)
    wy2 = w0y^2 * (1 + (z′ / zRy)^2)
    wx = sqrt(wx2)
    wy = sqrt(wy2)

    cutoff = x′^2 / wx2 + y′^2 / wy2
    if cutoff > 16
        return ComplexF64[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]
    end

    Rx = z′ == 0.0 ? Inf : z′ * (1 + (zRx / z′)^2)
    Ry = z′ == 0.0 ? Inf : z′ * (1 + (zRy / z′)^2)
    ζx = atan(z′ / zRx)
    ζy = atan(z′ / zRy)

    amp = sqrt(2/(c*ε0)) * sqrt(b.I0) * (w0x / wx) * (w0y / wy)
    G = exp(- (x′^2 / wx2 + y′^2 / wy2))
    φ = k * z′ + k * (x′^2 / (2 * Rx) + y′^2 / (2 * Ry)) - (ζx + ζy)

    scalar = amp * G * cis(-φ) * b._coeff[]

    # NEW: vector field = scalar envelope × polarization vector
    return ComplexF64[scalar * b.pol[1],
                      scalar * b.pol[2],
                      scalar * b.pol[3]]
end

"""
    efield_scalar(b::GeneralGaussianBeam, r) -> ComplexF64

Scalar complex field amplitude at position `r`, without constructing the
polarization vector. Used by `force` to avoid a temporary allocation.
"""
@inline function efield_scalar(b::GeneralGaussianBeam, r::Vector{Float64})
    x′, y′ = local_coords(r, b)
    z′ = sum((r[i] - b.r0[i]) * b.k[i] for i in 1:3)

    λ = b.λ
    w0x, w0y = b.w0x, b.w0y
    k = 2π / λ
    zRx = π * w0x^2 / λ
    zRy = π * w0y^2 / λ

    wx2 = w0x^2 * (1 + (z′ / zRx)^2)
    wy2 = w0y^2 * (1 + (z′ / zRy)^2)

    if x′^2 / wx2 + y′^2 / wy2 > 16
        return 0.0 + 0.0im
    end

    wx = sqrt(wx2)
    wy = sqrt(wy2)
    Rx = z′ == 0.0 ? Inf : z′ * (1 + (zRx / z′)^2)
    Ry = z′ == 0.0 ? Inf : z′ * (1 + (zRy / z′)^2)
    ζx = atan(z′ / zRx)
    ζy = atan(z′ / zRy)

    amp = sqrt(2/(c*ε0)) * sqrt(b.I0) * (w0x / wx) * (w0y / wy)
    G   = exp(-(x′^2 / wx2 + y′^2 / wy2))
    φ   = k * z′ + k * (x′^2 / (2 * Rx) + y′^2 / (2 * Ry)) - (ζx + ζy)

    return amp * G * cis(-φ) * b._coeff[]
end

"""
    dEdr(b::GeneralGaussianBeam, r) -> (dE_dx, dE_dy, dE_dz)

Gradient \\(\\nabla E\\) of the complex electric field amplitude of an elliptical
paraxial Gaussian beam at global position `r`.

This function accounts for:

- z'-dependent waists, curvatures, and Gouy phases.
- All cross derivatives with respect to local coordinates.
- Arbitrary beam orientation via the chain rule to global coordinates.
- Elliptical transverse profile and 8-waist cutoff.

Returns a tuple of three complex numbers.
"""
@inline function dEdr(b::GeneralGaussianBeam, r::Vector{Float64})
    x′, y′ = local_coords(r, b)
    z′ = sum((r[i] - b.r0[i]) * b.k[i] for i in 1:3)

    λ        = b.λ
    w0x, w0y = b.w0x, b.w0y
    k        = 2π / λ
    zRx      = π * w0x^2 / λ
    zRy      = π * w0y^2 / λ
    wx2      = w0x^2 * (1 + (z′ / zRx)^2)
    wy2      = w0y^2 * (1 + (z′ / zRy)^2)
    wx       = sqrt(wx2)
    wy       = sqrt(wy2)
    wx3      = wx2 * wx
    wy3      = wy2 * wy
    cutoff   = x′^2 / wx2 + y′^2 / wy2

    if cutoff > 16
        return (0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)
    end

    Rx = z′ == 0.0 ? Inf : z′ * (1 + (zRx / z′)^2)
    Ry = z′ == 0.0 ? Inf : z′ * (1 + (zRy / z′)^2)
    ζx = atan(z′ / zRx)
    ζy = atan(z′ / zRy)

    amp = sqrt(2/(c*ε0)) * sqrt(b.I0) * (w0x / wx) * (w0y / wy)
    G   = exp(- (x′^2 / wx2 + y′^2 / wy2))
    φ   = k * z′ + k * (x′^2 / (2 * Rx) + y′^2 / (2 * Ry)) - (ζx + ζy)
    cphase = cis(-φ) * b._coeff[]
    E      = amp * G * cphase

    # Derivatives with respect to x′ and y′
    dφ_dx′ = k * x′ / Rx
    dE_dx′ = E * (-2 * x′ / wx2 - im * dφ_dx′)
    dφ_dy′ = k * y′ / Ry
    dE_dy′ = E * (-2 * y′ / wy2 - im * dφ_dy′)

    # Derivative with respect to z′
    d_wx_dz′ = w0x * (z′ / zRx^2) / sqrt(1 + (z′ / zRx)^2)
    d_wy_dz′ = w0y * (z′ / zRy^2) / sqrt(1 + (z′ / zRy)^2)
    d_amp_dz′ = amp * (-d_wx_dz′ / wx - d_wy_dz′ / wy)
    dG_dz′ = G * (2 * x′^2 / wx3 * d_wx_dz′ + 2 * y′^2 / wy3 * d_wy_dz′)

    dRx_dz′ = z′ != 0.0 ? (1 + (zRx / z′)^2) - 2 * zRx^2 / z′^3 : 0.0
    dRy_dz′ = z′ != 0.0 ? (1 + (zRy / z′)^2) - 2 * zRy^2 / z′^3 : 0.0
    dζx_dz′ = zRx / (zRx^2 + z′^2)
    dζy_dz′ = zRy / (zRy^2 + z′^2)
    dφ_dz′  = k + k * (x′^2 / (2 * Rx^2) * dRx_dz′ +
                       y′^2 / (2 * Ry^2) * dRy_dz′) -
              (dζx_dz′ + dζy_dz′)

    dE_dz′ = (d_amp_dz′ * G + amp * dG_dz′) * cphase - im * dφ_dz′ * E

    # Chain rule to global coordinates
    dEdr_vec = (
        dE_dx′ * b.u[1] + dE_dy′ * b.v[1] + dE_dz′ * b.k[1],
        dE_dx′ * b.u[2] + dE_dy′ * b.v[2] + dE_dz′ * b.k[2],
        dE_dx′ * b.u[3] + dE_dy′ * b.v[3] + dE_dz′ * b.k[3],
    )

    return dEdr_vec
end