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
     integrand_ξ_GNCxLD_Lensing_Lensing(
          IP1::Point, IP2::Point,
          P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the Lensing auto-correlation function 
``\\xi^{\\kappa\\kappa} (s_1, s_2, \\cos{\\theta})``, i.e. the function 
``f(s_1, s_2, y, \\chi_1, \\chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \\chi_1, \\chi_2) = 
\\frac{1}{2}
\\frac{
     \\mathcal{H}_0^4 \\Omega_{ \\mathrm{M0}}^2 D_1 D_2 (\\chi_1 - s_1)(\\chi_2 - s_2)
}{
     s_1 s_2 a(\\chi_1) a(\\chi_2) }
(J_{00} \\, I^0_0(\\chi) + J_{02} \\, I^0_2(\\chi) + 
     J_{31} \\, I^3_1(\\chi) + J_{22} \\, I^2_2(\\chi))
```

where ``D_1 = D(\\chi_1)``, ``D_2 = D(\\chi_2)`` and so on, ``\\mathcal{H} = a H``, 
``\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2\\chi_1\\chi_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
    J_{00} & = - \\frac{3 \\chi_1^2 \\chi_2^2}{4 \\chi^4} (y^2 - 1) 
               (8 y (\\chi_1^2 + \\chi_2^2) - 9 \\chi_1 \\chi_2 y^2 - 7 \\chi_1 \\chi_2) \\\\
    J_{02} & = - \\frac{3 \\chi_1^2 \\chi_2^2}{2 \\chi^4} (y^2 - 1)
               (4 y (\\chi_1^2 + \\chi_2^2) - 3 \\chi_1 \\chi_2 y^2 - 5 \\chi_1 \\chi_2) \\\\
    J_{31} & = 9 y \\chi^2 \\\\
    J_{22} & = \\frac{9 \\chi_1 \\chi_2}{4 \\chi^4}
               [ 2 (\\chi_1^4 + \\chi_2^4) (7 y^2 - 3) 
                 - 16 y \\chi_1 \\chi_2 (\\chi_1^2 + \\chi_2^2) (y^2+1) 
               + \\chi_1^2 \\chi_2^2 (11 y^4 + 14 y^2 + 23)]
\\end{align*}
```

## Inputs

- `IP1::Point` and `IP2::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1` and `χ2` respectively.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments

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


See also: [`ξ_GNCxLD_Lensing_Lensing`](@ref), [`integrand_on_mu_Lensing`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
function integrand_ξ_GNCxLD_Lensing_Lensing(
     IP1::Point, IP2::Point,
     P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1 = P1.comdist
     s2 = P2.comdist
     χ1, D1, a1 = IP1.comdist, IP1.D, IP1.a
     χ2, D2, a2 = IP2.comdist, IP2.D, IP2.a
     s_b_s1 = cosmo.params.s_b
     Ω_M0 = cosmo.params.Ω_M0

     Δχ_square = χ1^2 + χ2^2 - 2 * χ1 * χ2 * y
     Δχ = Δχ_square > 0 ? √(Δχ_square) : 0

     denomin = s1 * s2 * a1 * a2
     factor = - ℋ0^4 * Ω_M0^2 * D1 * (s1 - χ1) * D2 * (s2 - χ2) * (5 * s_b_s1 - 2)

     χ1χ2 = χ1 * χ2

     new_J00 = - 3/4 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
     new_J02 = - 3/2 * χ1χ2^2 / Δχ^4 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
     new_J31 = 9 * y * Δχ^2
     new_J22 = 9/4 * χ1χ2 / Δχ^4 * (
          2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
          - 16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
          + χ1χ2^2 * (11y^4 + 14y^2 + 23)
     )

     I00 = cosmo.tools.I00(Δχ)
     I20 = cosmo.tools.I20(Δχ)
     I13 = cosmo.tools.I13(Δχ)
     I22 = cosmo.tools.I22(Δχ)

     res = new_J00 * I00 + new_J02 * I20 + new_J31 * I13 + new_J22 * I22

     return factor / denomin * res
end

function integrand_ξ_GNCxLD_Lensing_Lensing(
     χ1::Float64, χ2::Float64,
     s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP1, IP2 = Point(χ1, cosmo), Point(χ2, cosmo)
     return integrand_ξ_GNCxLD_Lensing_Lensing(IP1, IP2, P1, P2, y, cosmo; kwargs...)
end

#=
function func_Δχ_min(s1, s2, y; frac = 1e-4)
     Δs = s(s1, s2, y)
     return frac * Δs
end
=#

function ξ_GNCxLD_Lensing_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs_2::Int = 100)

     χ1s = P1.comdist .* range(1e-6, 1.0, length = N_χs_2)
     χ2s = P2.comdist .* range(1e-6, 1.0, length = N_χs_2 + 7)

     IP1s = [GaPSE.Point(x, cosmo) for x in χ1s]
     IP2s = [GaPSE.Point(x, cosmo) for x in χ2s]

     int_ξ_Lensings = [
          en * GaPSE.integrand_ξ_GNCxLD_Lensing_Lensing(IP1, IP2, P1, P2, y, cosmo)
          for IP1 in IP1s, IP2 in IP2s
     ]

     res = trapz((χ1s, χ2s), int_ξ_Lensings)
     #println("res = $res")
     return res / en
end


function ξ_GNCxLD_Lensing_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_GNCxLD_Lensing_Lensing(P1, P2, y, cosmo; kwargs...)
end



"""
     ξ_GNCxLD_Lensing_Lensing(P1::Point, P2::Point, y, cosmo::Cosmology;
          en::Float64 = 1e6,
          N_χs::Int = 100) :: Float64

     ξ_GNCxLD_Lensing_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...) = 
          ξ_GNCxLD_Lensing_Lensing(Point(s1, cosmo), Point(s2, cosmo), y, cosmo; kwargs...)

          
Return the Lensing auto-correlation function 
``\\xi^{\\kappa\\kappa} (s_1, s_2, \\cos{\\theta})`` concerning the perturbed
luminosity distance, defined as follows:
    
```math
\\xi^{\\kappa\\kappa} (s_1, s_2, \\cos{\\theta}) = 
\\int_0^{s_1} \\mathrm{d} \\chi_1 \\int_0^{s_2} \\mathrm{d} \\chi_2 
\\frac{1}{2}
\\frac{
     \\mathcal{H}_0^4 \\Omega_{ \\mathrm{M0}}^2 D_1 D_2 (\\chi_1 - s_1)(\\chi_2 - s_2)
}{
     s_1 s_2 a(\\chi_1) a(\\chi_2) }
(J_{00} \\, I^0_0(\\chi) + J_{02} \\, I^0_2(\\chi) + 
     J_{31} \\, I^3_1(\\chi) + J_{22} \\, I^2_2(\\chi))
```

where ``D_1 = D(\\chi_1)``, ``D_2 = D(\\chi_2)`` and so on, ``\\mathcal{H} = a H``, 
``\\chi = \\sqrt{\\chi_1^2 + \\chi_2^2 - 2\\chi_1\\chi_2\\cos{\\theta}}``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\\begin{align*}
    J_{00} & = - \\frac{3 \\chi_1^2 \\chi_2^2}{4 \\chi^4} (y^2 - 1) 
               (8 y (\\chi_1^2 + \\chi_2^2) - 9 \\chi_1 \\chi_2 y^2 - 7 \\chi_1 \\chi_2) \\\\
    J_{02} & = - \\frac{3 \\chi_1^2 \\chi_2^2}{2 \\chi^4} (y^2 - 1)
               (4 y (\\chi_1^2 + \\chi_2^2) - 3 \\chi_1 \\chi_2 y^2 - 5 \\chi_1 \\chi_2) \\\\
    J_{31} & = 9 y \\chi^2 \\\\
    J_{22} & = \\frac{9 \\chi_1 \\chi_2}{4 \\chi^4}
               [ 2 (\\chi_1^4 + \\chi_2^4) (7 y^2 - 3) 
                 - 16 y \\chi_1 \\chi_2 (\\chi_1^2 + \\chi_2^2) (y^2+1) 
               + \\chi_1^2 \\chi_2^2 (11 y^4 + 14 y^2 + 23)]
\\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_GNCxLD_Lensing_Lensing`.



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
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_GNCxLD_Lensing_Lensing`](@ref), [`integrand_on_mu_Lensing`](@ref)
[`integral_on_mu`](@ref), [`ξ_GNC_multipole`](@ref)
"""
ξ_GNCxLD_Lensing_Lensing



##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_Lensing_Lensing(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNCxLD_Lensing_Lensing(s2, s1, y, cosmo; kwargs...)
end

