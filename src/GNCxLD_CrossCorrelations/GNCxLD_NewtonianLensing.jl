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
     integrand_ξ_GNCxLD_Newtonian_Lensing(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the Doppler-LocalGP cross-correlation function 
``\\xi^{v_{\\parallel}\\int\\phi} (s_1, s_2, \\cos{\\theta})``, i.e. the function 
``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = 
     3 \\mathcal{H}(s_1) f(s_1) D(s_1) \\mathcal{H_0}^2 \\Omega_{M0} 
     \\mathcal{R}(s_1) J_{31} I^3_1(\\chi)
```
where ``\\mathcal{H} = a H``, 
``\\chi = \\sqrt{s_1^2 + \\chi_2^2 - 2 s_1 \\chi_2 \\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and:

```math
J_{31} = 
     \\frac{D(\\chi_2) (s_1 - \\chi_2 \\cos{\\theta})}{a(\\chi_2)} \\chi^2 
     \\left(
          \\frac{1}{s_2} - \\mathcal{R}(s_2) \\mathcal{H}(\\chi_2) (f(\\chi_2) - 1)
     \\right)
```

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_GNCxLD_Newtonian_Lensing`](@ref), [`int_on_mu_Newtonian_Lensing`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function integrand_ξ_GNCxLD_Newtonian_Lensing(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1, D_s1, f_s1 = P1.comdist, P1.D, P1.f
     s2 = P2.comdist
     χ2, D2, a2 = IP.comdist, IP.D, IP.a
     b_s1 = cosmo.params.b
     Ω_M0 = cosmo.params.Ω_M0

     Δχ2_square = s1^2 + χ2^2 - 2 * s1 * χ2 * y
     Δχ2 = Δχ2_square > 0 ? √(Δχ2_square) : 0

     common = - D_s1 * ℋ0^2 * Ω_M0 * D2 * (χ2 - s2) / (a2 * s2)

     new_J00 = 1 / 5 * (f_s1 * χ2 * (3 * y^2 - 1) - 3 * y * s1 * f_s1 - 5 * y * s1 * b_s1)
     new_J02 = 1 / (14 * Δχ2^2) * (
          7 * s1 * b_s1 * (-2 * χ2^2 * y + χ2 * s1 * (y^2 + 3) - 2 * y * s1^2) +
          f_s1 * (
               4 * χ2^3 * (3 * y^2 - 1) - 2 * χ2^2 * y * s1 * (3 * y^2 + 8)
               +
               χ2 * s1^2 * (9 * y^2 + 11) - 6 * y * s1^3
          )
     )
     new_J04 = 3 / (70 * Δχ2^4) * f_s1 * (
                    χ2^5 * (6 * y^2 - 2) + 6 * χ2^4 * y * s1 * (y^2 - 3)
                    -
                    χ2^3 * s1^2 * (y^4 + 12 * y^2 - 21)
                    +
                    2 * χ2^2 * y * s1^3 * (y^2 + 3) - 12 * χ2 * s1^4
                    +
                    4 * y * s1^5
               )

     I00 = cosmo.tools.I00(Δχ2)
     I20 = cosmo.tools.I20(Δχ2)
     I40 = cosmo.tools.I40(Δχ2)

     return common * (new_J00 * I00 + new_J02 * I20 + new_J04 * I40)
end


function integrand_ξ_GNCxLD_Newtonian_Lensing(
     χ2::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ2, cosmo)
     return integrand_ξ_GNCxLD_Newtonian_Lensing(IP, P1, P2, y, cosmo; kwargs...)
end


@doc raw"""
     ξ_GNCxLD_Newtonian_Lensing(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Int = 100):: Float64

Return the Doppler-LocalGP cross-correlation function 
``\\xi^{v_{\\parallel}\\int\\phi} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{v_{\\parallel}\\int\\phi} (s_1, s_2, \\cos{\\theta}) = 
     3 \\mathcal{H}(s_1) f(s_1) D(s_1) \\mathcal{H_0}^2 \\Omega_{M0} \\mathcal{R}(s_1) 
     \\int_0^{s_2} \\mathrm{d}\\chi_2 \\,  J_{31} \\,  I^3_1(\\chi)
```

where ``\\mathcal{H} = a H``, 
``\\chi = \\sqrt{s_1^2 + \\chi_2^2 - 2 s_1 \\chi_2 \\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and:

```math
J_{31} = 
     \\frac{D(\\chi_2) (s_1 - \\chi_2 \\cos{\\theta})}{a(\\chi_2)} \\chi^2 
     \\left(
          - \\frac{1}{s_2} + \\mathcal{R}(s_2) \\mathcal{H}(\\chi_2) (f(\\chi_2) - 1)
     \\right)
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNCxLD_Newtonian_Lensing`.


## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `N_χs::Int = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_GNCxLD_Newtonian_Lensing`](@ref), [`int_on_mu_Newtonian_Lensing`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function ξ_GNCxLD_Newtonian_Lensing(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100)

     χ2s = s2 .* range(1e-6, 1.0, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNCxLD_Newtonian_Lensing(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ2s, int_ξs)
     #println("res = $res")
     return res / en
end



##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_Lensing_Newtonian(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNCxLD_Newtonian_Lensing(s2, s1, y, cosmo; kwargs...)
end

