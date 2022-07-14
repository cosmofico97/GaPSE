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


"""
     integrand_ξ_GNC_Lensing_LocalGP(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the Lensing-LocalGP cross-correlation function 
``\\xi^{\\kappa \\phi} (s_1, s_2, \\cos{\\theta})``, i.e. the function 
``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = 
     \\frac{
          9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(s_2) (1 + \\mathcal{R}(s_2)) s_2
     }{4 a(s_2) s_1} 
     \\frac{D(\\chi_1)(s_1 - \\chi_1) }{a(\\chi_1)}
     \\left( J_{31} I^3_1(\\Delta\\chi_1) +  J_{22} I^2_2(\\Delta\\chi_1) \\right)
```

where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_1 = \\sqrt{\\chi_1^2 + s_2^2 - 2 \\chi_1 s_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
     J_{31} & = -2 y \\Delta\\chi_1^2 \\\\
     J_{22} & = \\chi_1 s_2 (1 - y^2)
\\end{align*}
```

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_GNC_Lensing_LocalGP`](@ref), [`int_on_mu_Lensing_LocalGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function integrand_ξ_GNC_Lensing_LocalGP(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1 = P1.comdist
     s2, D_s2, f_s2, a_s2, ℋ_s2, ℛ_s2 = P2.comdist, P2.D, P2.f, P2.a, P2.ℋ, P2.ℛ_GNC
     χ1, D1, a1 = IP.comdist, IP.D, IP.a
     s_b_s1, s_b_s2 = cosmo.params.s_b, cosmo.params.s_b
     𝑓_evo_s2 = cosmo.params.𝑓_evo
     Ω_M0 = cosmo.params.Ω_M0

     Δχ1_square = χ1^2 + s2^2 - 2 * χ1 * s2 * y
     Δχ1 = Δχ1_square > 0 ? √(Δχ1_square) : 0

     common = D_s2 * ℋ0^2 * Ω_M0 * s2 * D1 * (χ1 - s1) *  (5 * s_b_s1 - 2) * (
                   2 * f_s2 * a_s2 * ℋ_s2^2 * (𝑓_evo_s2 - 3)
                   + 3 * ℋ0^2 * Ω_M0 * (f_s2 + ℛ_s2 + 5 * s_b_s2 - 2)
              ) / (a1 * a_s2 * s1)
     factor = 2 * y * χ1^2 - χ1 * s2 * (y^2 + 3) + 2 * y * s2^2

     J20 = 1/2 * y * Δχ1^2

     I00 = cosmo.tools.I00(Δχ1)
     I20 = cosmo.tools.I20(Δχ1)
     I40 = cosmo.tools.I40(Δχ1)
     I02 = cosmo.tools.I02(Δχ1)

     return common * (
                factor * (1 / 60 * I00 + 1 / 42 * I20 + 1 / 140 * I40)
                + J20 * I02
           )
end


function integrand_ξ_GNC_Lensing_LocalGP(
     χ1::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ1, cosmo)
     return integrand_ξ_GNC_Lensing_LocalGP(IP, P1, P2, y, cosmo)
end


"""
     ξ_GNC_Lensing_LocalGP(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Int = 100):: Float64

Return the Lensing-LocalGP cross-correlation function 
``\\xi^{\\kappa \\phi} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{\\kappa \\phi} (s_1, s_2, \\cos{\\theta}) = 
     \\frac{
          9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(s_2) (1 + \\mathcal{R}(s_2)) s_2
     }{4 a(s_2) s_1} 
     \\int_0^{s_1} \\mathrm{d}\\chi_1 \\frac{D(\\chi_1)(s_1 - \\chi_1) }{a(\\chi_1)}
     \\left( J_{31} I^3_1(\\Delta\\chi_1) +  J_{22} I^2_2(\\Delta\\chi_1) \\right)
```

where ``\\mathcal{H} = a H``, 
``\\Delta\\chi_1 = \\sqrt{\\chi_1^2 + s_2^2 - 2 \\chi_1 s_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
     J_{31} & = -2 y \\Delta\\chi_1^2 \\\\
     J_{22} & = \\chi_1 s_2 (1 - y^2)
\\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNC_Lensing_LocalGP`.


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


See also: [`integrand_ξ_GNC_Lensing_LocalGP`](@ref), [`int_on_mu_Lensing_LocalGP`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function ξ_GNC_Lensing_LocalGP(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Int = 100)

     χ1s = s1 .* range(1e-6, 1, length = N_χs)

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ1s]

     int_ξs = [
          en * GaPSE.integrand_ξ_GNC_Lensing_LocalGP(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ1s, int_ξs)
     #println("res = $res")
     return res / en
end




##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_GNC_LocalGP_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNC_Lensing_LocalGP(s2, s1, y, cosmo; kwargs...)
end

