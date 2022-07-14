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

function ξ_PP_Doppler_L0(P::Point, cosmo::Cosmology)
     Peff = Point(cosmo.s_eff, cosmo)
     D, f, ℋ, ℛ = Peff.D, Peff.f, Peff.ℋ, Peff.ℛ_LD
     s = P.comdist

     1.0 / 3.0 * f^2 * ℋ^2 * ℛ^2 * D^2 * s^2 * cosmo.tools.I02(s)
end

function ξ_PP_Doppler_L0(s1, cosmo::Cosmology)
     P1 = Point(s1, cosmo)
     return ξ_PP_Doppler_L0(P1, cosmo)
end

function ξ_PP_Doppler_L2(P::Point, cosmo::Cosmology)
     Peff = Point(cosmo.s_eff, cosmo)
     D, f, ℋ, ℛ = Peff.D, Peff.f, Peff.ℋ, Peff.ℛ_LD
     s = P.comdist
     -2.0 / 3.0 * f^2 * ℋ^2 * ℛ^2 * D^2 * s^2 * cosmo.tools.I22(s)
end

function ξ_PP_Doppler_L2(s1, cosmo::Cosmology)
     P1 = Point(s1, cosmo)
     return ξ_PP_Doppler_L2(P1, cosmo)
end

function ξ_PP_Doppler(s1, y, cosmo::Cosmology)
     return ξ_PP_Doppler_L0(s1, cosmo) + ξ_PP_Doppler_L2(s1, cosmo) * Pl(y, 2)
end



##########################################################################################92



function integrand_ξ_PPD_multipole(s, μ, cosmo::Cosmology;
     L::Int=0, use_windows::Bool=true)

     res = if use_windows == true
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = ξ_PP_Doppler(s, μ, cosmo)
          #println("int = $int")
          int .* (spline_integrF(s, μ, cosmo.windowFint)/cosmo.WFI_norm * Pl(μ, L))
     else
          #println("s1 = $s1 \\t s2 = $(s2(s1, s, μ)) \\t  y=$(y(s1, s, μ))")
          int = ξ_PP_Doppler(s, μ, cosmo)
          #println("int = $int")
          #println( "Pl(μ, L) = $(Pl(μ, L))")
          int .* Pl(μ, L)
     end

     #println("res = $res")
     return (2.0 * L + 1.0) / 2.0 * res
end

function ξ_PPD_multipole(
     s, cosmo::Cosmology;
     L::Int=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     μ_atol::Float64=0.0,
     μ_rtol::Float64=1e-2)

     orig_f(μ) = enhancer * integrand_ξ_PPD_multipole(s, μ, cosmo;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=μ_atol, rtol=μ_rtol)[1]

     return int / enhancer
end



##########################################################################################92



function map_ξ_PPD_multipole(cosmo::Cosmology,
     v_ss=nothing;
     pr::Bool=true,
     N_log::Int=1000,
     L::Int=0,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length=N_log) : v_ss
     xis = pr ? begin
          @showprogress "PP Doppler, L=$L: " [
               ξ_PPD_multipole(s, cosmo; L=L, kwargs...) for s in ss
          ]
     end : [
          ξ_PPD_multipole(s, cosmo; L=L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_PPD_multipole for PP_Doppler " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


##########################################################################################92



function print_map_ξ_PPD_multipole(
     cosmo::Cosmology,
     out::String,
     v_ss=nothing;
     kwargs...)

     t1 = time()
     vec = map_ξ_PPD_multipole(cosmo, v_ss; kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of the ξ multipole concerning the PP Doppler.")
          parameters_used(io, cosmo)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "# ")
          println(io, "# s [Mpc/h_0] \t \t xi")
          for (s, xi) in zip(vec[1], vec[2])
               println(io, "$s \t $xi")
          end
     end
end


