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
     integrand_ξ_lensingdoppler(
          IP::Point, P1::Point, P2::Point,
          y, cosmo::Cosmology) :: Float64

Return the integrand of the lensing auto-correlation function 
``\xi^{v_{\parallel}\kappa} (s_1, s_2, \cos{\theta})``, i.e. the function 
``f(s_1, s_2, y, \chi_1, \chi_2)`` defined as follows:  

```math
f(s_1, s_2, y, \chi_1, \chi_2) = 
     \mathcal{H}_0^2 \Omega_{M0} D(s_2) f(s_2) \mathcal{H}(s_2) \mathcal{R}(s_2) 
     \frac{ D(\chi_1) (\chi_1 - s_1) }{a(\chi_1) s_1} 
     \left(
          J_{00} I^0_0(\chi) + J_{02} I^0_2(\chi) + J_{04} I^0_4(\chi) + J_{20} I^2_0(\chi)
     \right)
```

where ``\mathcal{H} = a H``, 
``\chi = \sqrt{\chi_1^2 + s_2^2 - 2\chi_1s_2\cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``) 
and the ``J`` coefficients are given by 

```math
\begin{align*}
     J_{00} & = \frac{1}{15}(\chi_1^2 y + \chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \\
     J_{02} & = \frac{1}{42 \chi^2} 
          (4 \chi_1^4 y + 4 \chi_1^3 (2 y^2 - 3) s_2 + \chi_1^2 y (11 - 23 y^2) s_2^2 + 
          \chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \\
     J_{04} & = \frac{1}{70 \chi^2}
          (2 \chi_1^4 y + 2 \chi_1^3 (2y^2 - 3) s_2 - \chi_1^2 y (y^2 + 5) s_2^2 + 
          \chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \\
     J_{20} & = y \chi^2
\end{align*}
```

## Inputs

- `IP::Point`: `Point` inside the integration limits, placed 
  at comoving distance `χ1`.

- `P1::Point` and `P2::Point`: extreme `Point` of the integration, placed 
  at comoving distance `s1` and `s2` respectively.

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


See also: [`ξ_lensingdoppler`](@ref), [`int_on_mu_lensingdoppler`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function integrand_ξ_lensingdoppler(
     IP::Point, P1::Point, P2::Point,
     y, cosmo::Cosmology)

     s1 = P1.comdist
     s2, D_s2, f_s2, ℋ_s2, ℛ_s2 = P2.comdist, P2.D, P2.f, P2.ℋ, P2.ℛ
     χ1, D1, a_χ1 = IP.comdist, IP.D, IP.a
     Ω_M0 = cosmo.params.Ω_M0

     Δχ = √(χ1^2 + s2^2 - 2 * χ1 * s2 * y)

     factor = D1 * (χ1 - s1) / (s1 * a_χ1)
     prefactor = ℋ0^2 * Ω_M0 * D_s2 * f_s2 * ℋ_s2 * ℛ_s2

     new_J00 = 1.0 / 15.0 * (χ1^2 * y + χ1 * (4 * y^2 - 3) * s2 - 2 * y * s2^2)
     new_J02 = 1.0 / (42 * Δχ^2) * (
          4 * χ1^4 * y + 4 * χ1^3 * (2 * y^2 - 3) * s2
          + χ1^2 * y * (11 - 23 * y^2) * s2^2 
          + χ1 * (23 * y^2 - 3) * s2^3 - 8 * y * s2^4)
     new_J04 = 1.0 / (70 * Δχ^2) *(
          2 * χ1^4 * y + 2 * χ1^3 * (2 * y^2 - 3) * s2
          - χ1^2 * y * (y^2 + 5) * s2^2 
          + χ1 * (y^2 + 9) * s2^3 - 4 * y * s2^4)

     new_J20 = y * Δχ^2

     I00 = cosmo.tools.I00(Δχ)
     I20 = cosmo.tools.I20(Δχ)
     I40 = cosmo.tools.I40(Δχ)
     I02 = cosmo.tools.I02(Δχ)

     #println("J00 = $new_J00, \t I00(Δχ) = $(I00)")
     #println("J02 = $new_J02, \t I20(Δχ) = $(I20)")
     #println("J31 = $new_J31, \t I13(Δχ) = $(I13)")
     #println("J22 = $new_J22, \t I22(Δχ) = $(I22)")

     parenth = (
          new_J00 * I00 + new_J02 * I20 +
          new_J04 * I40 + new_J20 * I02
     )


     res =  prefactor * factor * parenth
     #println("res = ", res, "\n")
     return res
end


function integrand_ξ_lensingdoppler(
     χ1::Float64, s1::Float64, s2::Float64,
     y, cosmo::Cosmology;
     kwargs...)

     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     IP = Point(χ1, cosmo)
     return integrand_ξ_lensing(IP, P1, P2, y, cosmo; kwargs...)
end


@doc raw"""
     ξ_lensingdoppler(s1, s2, y, cosmo::Cosmology;
          en::Float64 = 1e6, N_χs::Integer = 100):: Float64

Return the lensing auto-correlation function 
``\xi^{\kappa\kappa} (s_1, s_2, \cos{\theta})``, defined as follows:
    
```math
\xi^{v_{\parallel}\kappa} (s_1, s_2, \cos{\theta}) = 
     \mathcal{H}_0^2 \Omega_{M0} D(s_2) f(s_2) \mathcal{H}(s_2) \mathcal{R}(s_2) 
     \int_0^{s_1} \mathrm{d} \chi_1 
     \frac{ D(\chi_1) (\chi_1 - s_1) }{a(\chi_1) s_1} 
     \left(
          J_{00} I^0_0(\chi) + J_{02} I^0_2(\chi) + J_{04} I^0_4(\chi) + J_{20} I^2_0(\chi)
     \right)
```

where ``\mathcal{H} = a H``, 
``\chi = \sqrt{\chi_1^2 + s_2^2 - 2 \chi_1 s_2 \cos{\theta}}``, 
``y = \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2``) 
and the ``J`` coefficients are given by:

```math
\begin{align*}
     J_{00} & = \frac{1}{15}(\chi_1^2 y + \chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \\
     J_{02} & = \frac{1}{42 \chi^2} 
          (4 \chi_1^4 y + 4 \chi_1^3 (2 y^2 - 3) s_2 + \chi_1^2 y (11 - 23 y^2) s_2^2 + 
          \chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \\
     J_{04} & = \frac{1}{70 \chi^2}
          (2 \chi_1^4 y + 2 \chi_1^3 (2y^2 - 3) s_2 - \chi_1^2 y (y^2 + 5) s_2^2 + 
          \chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \\
     J_{20} & = y \chi^2
\end{align*}
```

The computation is made applying [`trapz`](@ref) (see the 
[Trapz](https://github.com/francescoalemanno/Trapz.jl) Julia package) to
the integrand function `integrand_ξ_lensingdoppler`.



## Inputs

- `s1` and `s2`: comovign distances where the function must be evaluated

- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `Δχ_min::Float64 = 1e-6` : when ``\Delta\chi = \sqrt{\chi_1^2 + \chi_2^2 - 2 \, \chi_1 \chi_2 y} \to 0^{+}``,
  some ``I_\ell^n`` term diverges, but the overall parenthesis has a known limit:

  ```math
     \lim_{\chi\to0^{+}} (J_{00} \, I^0_0(\chi) + J_{02} \, I^0_2(\chi) + 
          J_{31} \, I^3_1(\chi) + J_{22} \, I^2_2(\chi)) = 
          \frac{4}{15} \, (5 \, \sigma_2 + \frac{2}{3} \, σ_0 \,s_1^2 \, \chi_2^2)
  ```

  So, when it happens that ``\chi < \Delta\chi_\mathrm{min}``, the function considers this limit
  as the result of the parenthesis instead of calculating it in the normal way; it prevents
  computational divergences.

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.


See also: [`integrand_ξ_lensingdoppler`](@ref), [`int_on_mu_lensingdoppler`](@ref)
[`integral_on_mu`](@ref), [`ξ_multipole`](@ref)
"""
function ξ_lensingdoppler(s1, s2, y, cosmo::Cosmology;
     en::Float64 = 1e6, N_χs::Integer = 100)

     adim_χs = range(1e-6, 1.0, N_χs)
     χ1s = adim_χs .* s1

     P1, P2 = GaPSE.Point(s1, cosmo), GaPSE.Point(s2, cosmo)
     IPs = [GaPSE.Point(x, cosmo) for x in χ1s]

     int_ξs = [
          en * GaPSE.integrand_ξ_lensingdoppler(IP, P1, P2, y, cosmo)
          for IP in IPs
     ]

     res = trapz(χ1s, int_ξs)
     #println("res = $res")
     return res / en
end



##########################################################################################92



@doc raw"""
     int_on_mu_lensingdoppler(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true, 
          en::Float64 = 1e6,
          N_χs::Integer = 100) :: Float64

Return the integrand on ``\mu = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}`` 
of the lensing auto-correlation function, i.e.
the following function ``f(s_1, s, \mu)``:

```math
     f(s_1, s, \mu) = \xi^{\kappa\kappa} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) \,  \phi(s_2) \, F\left(\frac{s}{s_1}, \mu \right)
```
where ``y =  \cos{\theta} = \hat{\mathbf{s}}_1 \dot \hat{\mathbf{s}}_2`` and
``s = \sqrt{s_1^2 + s_2^2 - 2 \, s_1 \, s_2 \, y}``.

In case `use_windows` is set to `false`, the window functions ``\phi`` and ``F``
are removed, i.e is returned the following function ``f^{'}(s_1, s, \mu)``:

```math
     f^{'}(s_1, s, \mu) = \xi^{\kappa\kappa} (s_1, s_2, \cos{\theta}) 
          \, \mathcal{L}_L(\mu) 
```

The function ``\xi^{\kappa\kappa}(s_1, s_2, \cos{\theta})`` is calculated
from `ξ_lensing`; note that these is an internal conversion of coordiate sistems
from `(s1, s, μ)` to `(s1, s2, y)` thorugh the functions `y` and `s2`

## Inputs

- `s1`: the comoving distance where must be evaluated the integral

- `s`: the comoving distance from `s1` where must be evaluated the integral

- `μ`: the cosine between `s1` and `s` where must be evaluated the integral

- `cosmo::Cosmology`: cosmology to be used in this computation


## Optional arguments 

- `L::Integer = 0`: order of the Legendre polynomial to be used

- `en::Float64 = 1e6`: just a float number used in order to deal better 
  with small numbers;

- `use_windows::Bool = false`: tells if the integrand must consider the two
   window function ``\phi`` and ``F``

- `N_χs::Integer = 100`: number of points to be used for sampling the integral
  along the ranges `(0, s1)` (for `χ1`) and `(0, s1)` (for `χ2`); it has been checked that
  with `N_χs ≥ 50` the result is stable.

See also: [`integrand_ξ_lensingdoppler`](@ref), [`ξ_lensingdoppler`](@ref),
[`integral_on_mu`](@ref), [`map_integral_on_mu`](@ref),
[`spline_F`](@ref), [`ϕ`](@ref), [`Cosmology`](@ref), 
[`y`](@ref), [`s2`](@ref)
"""
function int_on_mu_lensingdoppler(s1, s, μ, cosmo::Cosmology;
          L::Integer = 0, 
          use_windows::Bool = true, 
          en::Float64 = 1e6,
          N_χs::Integer = 100)

     s2_value = s2(s1, s, μ)
     y_value = y(s1, s, μ)
     res = if use_windows == true
          ϕ_s2 = ϕ(s2_value; s_min = cosmo.s_min, s_max = cosmo.s_max)
          (ϕ_s2 > 0.0) || (return 0.0)
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensingdoppler(s1, s2_value, y_value, cosmo;
               en = en, N_χs = N_χs)
          #println("int = $int")
          int .* (ϕ_s2 * spline_F(s / s1, μ, cosmo.windowF) * Pl(μ, L))
     else
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensingdoppler(s1, s2_value, y_value, cosmo;
               en = en,  N_χs = N_χs)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return res
end
