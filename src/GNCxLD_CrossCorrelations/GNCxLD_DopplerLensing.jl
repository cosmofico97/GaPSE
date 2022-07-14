# -*- encoding: utf-8 -*-
#
# This file is part of GaPSE
# Copyright (C) 2022 Matteo Foglieni
#
# GaPSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GaPSE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GaPSE. If not, see <http://www.gnu.org/licenses/>.
#


@doc raw"""
     integrand_ξ_GNCxLD_Doppler_Lensing(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the Lensing-Doppler cross-correlation function 
``\\xi^{\\kappa v_{\\parallel}} (s_1, s_2, \\cos{\\theta})``, i.e. the function 
``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = 
     \\mathcal{H}_0^2 \\Omega_{M0} D(s_2) f(s_2) \\mathcal{H}(s_2) \\mathcal{R}(s_2) 
     \\frac{ D(\\chi_1) (\\chi_1 - s_1) }{a(\\chi_1) s_1} 
     \\left(
          J_{00} I^0_0(\\Delta\\chi_1) + J_{02} I^0_2(\\Delta\\chi_1) 
          + J_{04} I^0_4(\\Delta\\chi_1) + J_{20} I^2_0(\\Delta\\chi_1)
     \\right)
```

where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_1 = \\sqrt{\\chi_1^2 + s_2^2 - 2\\chi_1s_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
     J_{00} & = \\frac{1}{15}(\\chi_1^2 y + \\chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \\\\
     J_{02} & = \\frac{1}{42 \\Delta\\chi_1^2} 
          (4 \\chi_1^4 y + 4 \\chi_1^3 (2 y^2 - 3) s_2 + \\chi_1^2 y (11 - 23 y^2) s_2^2 + 
          \\chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \\\\
     J_{04} & = \\frac{1}{70 \\Delta\\chi_1^2}
          (2 \\chi_1^4 y + 2 \\chi_1^3 (2y^2 - 3) s_2 - \\chi_1^2 y (y^2 + 5) s_2^2 + 
          \\chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \\\\
     J_{20} & = y \\Delta\\chi_1^2
\\end{align*}
```

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ2`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_GNCxLD_Doppler_Lensing`](@ref), [`int_on_mu_Lensing_Doppler`](@ref)
[`integral_on_mu`](@ref), [`ξ_LD_multipole`](@ref)
"""
function integrand_ξ_GNCxLD_Doppler_Lensing(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1, D_s1, f_s1, ℋ_s1, ℛ_s1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ_GNC
     s2 = P2.comdist
     χ2, D2, a2 = IP.comdist, IP.D, IP.a
     Ω_M0 = cosmo.params.Ω_M0


     Δχ2_square = χ2^2 + s1^2 - 2 * χ2 * s1 * y
     Δχ2 = Δχ2_square > 0.0 ? √(Δχ2_square) : 0.0

     common = - ℋ0^2 * Ω_M0 * D2 * (χ2 - s2) / (s2 * a2)
     factor = D_s1 * f_s1 * ℋ_s1 * ℛ_s1

     new_J00 = 1.0 / 15.0 * (χ2^2 * y + χ2 * (4 * y^2 - 3) * s1 - 2 * y * s1^2)
     new_J02 = 1.0 / (42 * Δχ2^2) * (
          4 * χ2^4 * y + 4 * χ2^3 * (2 * y^2 - 3) * s1
          + χ2^2 * y * (11 - 23 * y^2) * s1^2
          + χ2 * (23 * y^2 - 3) * s1^3 - 8 * y * s1^4)
     new_J04 = 1.0 / (70 * Δχ2^2) * (
          2 * χ2^4 * y + 2 * χ2^3 * (2 * y^2 - 3) * s1
          -
          χ2^2 * y * (y^2 + 5) * s1^2
          +
          χ2 * (y^2 + 9) * s1^3 - 4 * y * s1^4)
     new_J20 = y * Δχ2^2

     I00 = cosmo.tools.I00(Δχ2)
     I20 = cosmo.tools.I20(Δχ2)
     I40 = cosmo.tools.I40(Δχ2)
     I02 = cosmo.tools.I02(Δχ2)

     #println("J00 = $new_J00, \t I00(Δχ2) = $(I00)")
     #println("J02 = $new_J02, \t I20(Δχ2) = $(I20)")
     #println("J31 = $new_J31, \t I13(Δχ2) = $(I13)")
     #println("J22 = $new_J22, \t I22(Δχ2) = $(I22)")

     parenth = (
          new_J00 * I00 + new_J02 * I20 +
          new_J04 * I40 + new_J20 * I02
     )

     first = common * factor * parenth

     #new_J31 = -3 * χ2^2 * y * f0 * ℋ0
     #I13 = cosmo.tools.I13(χ2)
     #second = new_J31 * I13 * common

     return first
end


function integrand_ξ_GNCxLD_Doppler_Lensing(
     χ2::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ2, cosmo)
     return integrand_ξ_GNCxLD_Doppler_Lensing(IP, P1, P2, y, cosmo)
end


"""
     ξ_GNCxLD_Doppler_Lensing(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Int = 100):: Float64

Return the Lensing-Doppler cross-correlation function 
``\\xi^{\\kappa v_{\\parallel}} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{\\kappa v_{\\parallel}} (s_1, s_2, \\cos{\\theta}) = 
     \\mathcal{H}_0^2 \\Omega_{M0} D(s_2) f(s_2) \\mathcal{H}(s_2) \\mathcal{R}(s_2) 
     \\int_0^{s_1} \\mathrm{d} \\chi_1 
     \\frac{ D(\\chi_1) (\\chi_1 - s_1) }{a(\\chi_1) s_1} 
     \\left(
          J_{00} I^0_0(\\Delta\\chi_1) + J_{02} I^0_2(\\Delta\\chi_1) 
          + J_{04} I^0_4(\\Delta\\chi_1) + J_{20} I^2_0(\\Delta\\chi_1)
     \\right)
```

where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_1= \\sqrt{\\chi_1^2 + s_2^2 - 2 \\chi_1 s_2 \\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by:

```math
\\begin{align*}
     J_{00} & = \\frac{1}{15}(\\chi_1^2 y + \\chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \\\\
     J_{02} & = \\frac{1}{42 \\Delta\\chi_1^2} 
          (4 \\chi_1^4 y + 4 \\chi_1^3 (2 y^2 - 3) s_2 + \\chi_1^2 y (11 - 23 y^2) s_2^2 + 
          \\chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \\\\
     J_{04} & = \\frac{1}{70 \\Delta\\chi_1^2}
          (2 \\chi_1^4 y + 2 \\chi_1^3 (2y^2 - 3) s_2 - \\chi_1^2 y (y^2 + 5) s_2^2 + 
          \\chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \\\\
     J_{20} & = y \\Delta\\chi_1^2
\\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_LD_GNC-LD_Doppler_Lensing`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `Δχ_min::Float64 = 1e-6` : when ``\\Delta\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2 \\, \\chi_1 \\chi_2 y} \\to 0^{+}``,
  some ``I_\\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
     \\lim_{\\chi\\to0^{+}} (J_{00} \\, I^0_0(\\chi) + J_{02} \\, I^0_2(\\chi) + 
          J_{31} \\, I^3_1(\\chi) + J_{22} \\, I^2_2(\\chi)) = 
          \\frac{4}{15} \\, (5 \\, \\sigma_2 + \\frac{2}{3} \\, σ_0 \\,s_1^2 \\, \\chi_2^2)
  ```

  So, when it happens that ``\\chi < \\Delta\\chi_\\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

- `N_χs::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ2`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_GNCxLD_Doppler_Lensing`](@ref), [`int_on_mu_Lensing_Doppler`](@ref)
[`integral_on_mu`](@ref), [`ξ_LD_multipole`](@ref)
"""
function ξ_GNCxLD_Doppler_Lensing(s1, s2, y, cosmo::Cosmology;
     en::Float64=1e6, N_χs::Int=100)

     χ2s = s2 .* range(1e-6, 1.0, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNCxLD_Doppler_Lensing(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ2s, int_ξs)
     #println("res = $res")
     return res / en
end



##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_Lensing_Doppler(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNCxLD_Doppler_Lensing(s2, s1, y, cosmo; kwargs...)
end


