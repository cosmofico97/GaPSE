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
     integrand_ξ_LD_Doppler_IntegratedGP(
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


See also: [`ξ_LD_Doppler_IntegratedGP`](@ref), [`int_on_mu_Doppler_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_LD_multipole`](@ref)
"""
function integrand_ξ_LD_Doppler_IntegratedGP(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1, D_s1, f_s1, ℋ_s1, ℛ_s1 = P1.comdist, P1.D, P1.f, P1.ℋ, P1.ℛ_LD
     s2, ℛ_s2 = P2.comdist, P2.ℛ_LD
     χ2, D2, a2, f2, ℋ2 = IP.comdist, IP.D, IP.a, IP.f, IP.ℋ
     Ω_M0 = cosmo.params.Ω_M0

     Δχ2_square = s1^2 + χ2^2 - 2 * s1 * χ2 * y
     Δχ2 = Δχ2_square > 0 ? √(Δχ2_square) : 0.0

     common = 3 * ℋ_s1 * f_s1 * D_s1 * ℋ0^2 * Ω_M0 * ℛ_s1
     #common = ℋ0^2 * Ω_M0 * D2 / (s2 * a2)
     #factor = Δχ2^2 * (χ2 * y - s1) * (s2 * (f2 - 1) * ℋ2 * ℛ_s2 + 1)

     #I00 = cosmo.tools.I00(Δχ2)
     #I20 = cosmo.tools.I20(Δχ2)
     #I40 = cosmo.tools.I40(Δχ2)
     #I02 = cosmo.tools.I02(Δχ2)

     #first = common * factor * (1 / 15 * I00 + 2 / 21 * I20 + 1 / 35 * I40 + I02)

     #new_J31 = -3 * χ2^3 * y * f0 * ℋ0 * (ℛ_s1 + 1) * (s2 * (f2 - 1) * ℋ2 * ℛ_s2 + 1)
     new_J31 = Δχ2^2 * D2 * (χ2 * y - s1) / (a2 * s2)  * (s2 * ℛ_s2 * ℋ2 * (f2 - 1) - 1)
     I13 = cosmo.tools.I13(Δχ2)

     second = common * new_J31 * I13

     #println("J00 = $new_J00, \t I00(Δχ) = $(I00)")
     #println("J02 = $new_J02, \t I20(Δχ) = $(I20)")
     #println("J31 = $new_J31, \t I13(Δχ) = $(I13)")
     #println("J22 = $new_J22, \t I22(Δχ) = $(I22)")

     return second
end


function integrand_ξ_LD_Doppler_IntegratedGP(
     χ2::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ2, cosmo)
     return integrand_ξ_LD_Doppler_IntegratedGP(IP, P1, P2, y, cosmo; kwargs...)
end


@doc raw"""
     ξ_LD_Doppler_IntegratedGP(s1, s2, y, cosmo::Cosmology;
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
the integrand function `integrand_ξ_LD_Doppler_IntegratedGP`.


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


See also: [`integrand_ξ_LD_Doppler_IntegratedGP`](@ref), [`int_on_mu_Doppler_IntegratedGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_LD_multipole`](@ref)
"""
function ξ_LD_Doppler_IntegratedGP(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100)

     χ2s = range(1e-6, 1.0, length = N_χs) .* s2

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξs = [
          en * GaPSE.integrand_ξ_LD_Doppler_IntegratedGP(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ2s, int_ξs)
     #println("res = $res")
     return res / en
end




##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LD_IntegratedGP_Doppler(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_LD_Doppler_IntegratedGP(s2, s1, y, cosmo; kwargs...)
end

