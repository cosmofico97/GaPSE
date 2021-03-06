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


function ξ_GNCxLD_LocalGP_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology)
     s1, D1, f1, a1, ℛ1, ℋ1 = P1.comdist, P1.D, P1.f, P1.a, P1.ℛ_GNC, P1.ℋ
     s2, D2, a2, ℜ2 = P2.comdist, P2.D, P2.a, P2.ℛ_LD
     s_b1 = cosmo.params.s_b
     𝑓_evo1 = cosmo.params.𝑓_evo
     Ω_M0 = cosmo.params.Ω_M0

     Δs = s(s1, s2, y)

     factor = - 3 * Δs^4 * ℋ0^2 * Ω_M0 * D1 * D2 * (1 + ℜ2) / (4 * a1 * a2)
     parenth = 2 * f1 * ℋ1^2 * a1 * (𝑓_evo1 - 3) + 3 * ℋ0^2 * Ω_M0 * (f1 + ℛ1 + 5 * s_b1 - 2)

     I04_tilde = cosmo.tools.I04_tilde(Δs)

     res = factor * parenth * I04_tilde

     return res
end


function ξ_GNCxLD_LocalGP_LocalGP(s1, s2, y, cosmo::Cosmology)
     P1, P2 = Point(s1, cosmo), Point(s2, cosmo)
     return ξ_GNCxLD_LocalGP_LocalGP(P1, P2, y, cosmo)
end


"""
     ξ_GNCxLD_LocalGP_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64

     ξ_GNCxLD_LocalGP_LocalGP(s1, s2, y, cosmo::Cosmology) = 
          ξ_GNCxLD_LocalGP_LocalGP(Point(s1, cosmo), Point(s2, cosmo), y, cosmo::Cosmology)

Return the local gravitational potential auto-correlation function concerning the perturbed
luminosity distance, defined as follows:
```math
\\xi^{\\phi\\phi} (s_1, s_2, \\cos{\\theta}) = 
     \\frac{9 \\mathcal{H}_0^4 \\Omega_{M0}^2 D(s_1) D(s_2)s^4}{4 a(s_1) a(s_2)}
     (1 + \\mathcal{R}_1 + \\mathcal{R}_2 + \\mathcal{R}_1\\mathcal{R}_2)
     \\tilde{I}^4_0(s)
```
where ``D_1 = D(s_1)``, ``D_2 = D(s_2)`` and so on, ``\\mathcal{H} = a H``, 
``y = \\cos{\\theta} = \\hat{\\mathbf{s}}_1 \\cdot \\hat{\\mathbf{s}}_2`` and:
```math
\\tilde{I}^4_0 (s) &= \\int_0^\\infty \\frac{\\mathrm{d}q}{2\\pi^2} 
          q^2 \\, P(q) \\, \\frac{j_0(q s) - 1}{(q s)^4}
```

## Inputs

- `P1::Point` and `P2::Point`: `Point` where the CF has to be calculated; they contain all the 
  data of interest needed for this calculus (comoving distance, growth factor and so on)
     
- `y`: the cosine of the angle between the two points `P1` and `P2`

- `cosmo::Cosmology`: cosmology to be used in this computation

See also: [`Point`](@ref), [`Cosmology`](@ref)
"""
ξ_GNCxLD_LocalGP_LocalGP


##########################################################################################92

##########################################################################################92

##########################################################################################92



function ξ_LDxGNC_LocalGP_LocalGP(s1, s2, y, cosmo::Cosmology; kwargs...)
     ξ_GNCxLD_LocalGP_LocalGP(s2, s1, y, cosmo; kwargs...)
end
