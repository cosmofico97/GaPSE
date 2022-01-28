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

function integrand_ξ_lensing(χ1, χ2, s1, s2, y; tol = 0.5, enhancer = 1)
     (s(s1, s2, y) >= tol) || (return 0.0)

     Δχ = √(χ1^2 + χ2^2 - 2 * χ1 * χ2 * y)

     D1, D2 = D(χ1), D(χ2)
     a_χ1, a_χ2 = 1 / (1 + z_of_s(χ1)), 1 / (1 + z_of_s(χ2))
     #psb1, psb2 = -1.0, -1.0
     #s_b1, s_b2 = s_b(s1), s_b(s2)
     #psb1, psb2 = 5 * s_b1 - 2, 5 * s_b2 - 2

     denomin = s1 * s2 * Δχ^4 * a_χ1 * a_χ2
     factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2)
     #factor = ℋ0^4 * Ω_M0^2 * D1 * (χ1 - s1) * D2 * (χ2 - s2) * psb1 * psb2

     if Δχ > 1e-6
          χ1χ2 = χ1 * χ2

          new_J00 = -0.75 * χ1χ2^2 * (y^2 - 1) * (8 * y * (χ1^2 + χ2^2) - χ1χ2 * (9 * y^2 + 7))
          new_J02 = -1.5 * χ1χ2^2 * (y^2 - 1) * (4 * y * (χ1^2 + χ2^2) - χ1χ2 * (3 * y^2 + 5))
          new_J31 = 9 * y * Δχ^6
          new_J22 = 2.25 * χ1χ2 * (
                         2 * (χ1^4 + χ2^4) * (7 * y^2 - 3)
                         -
                         16 * y * χ1χ2 * (y^2 + 1) * (χ1^2 + χ2^2)
                         +
                         χ1χ2^2 * (11y^4 + 14y^2 + 23)
                    )


          return enhancer * factor / denomin * (
               new_J00 * I00(Δχ) + new_J02 * I20(Δχ) +
               new_J31 * I13(Δχ) + new_J22 * I22(Δχ)
          )
     else
          lim = 4.0 / 15.0 * (5.0 * σ_2 + 2.0 / 3.0 * σ_0 * s1^2 * χ2^2)

          return 9.0 / 4.0 * enhancer * factor / denomin * lim
     end
end


function ξ_lensing(s1, s2, y; tol = 0.5, enhancer = 1, kwargs...)
     my_int(var) = integrand_ξ_lensing(var[1], var[2], s1, s2, y; tol = tol, enhancer = enhancer)
     a = [0.0, 0.0]
     b = [s1, s2]
     int = hcubature(my_int, a, b; kwargs...)
     #println(int)
     return int
end

function int_on_mu_lensing(s1, s, μ; L::Integer = 0, enhancer = 1, tol = 0.5, χ_atol = 1e-3, χ_rtol = 1e-3)
     if ϕ(s2(s1, s, μ)) > 0
          #println("s1 = $s1 \t s2 = $(s2(s1, s, μ)) \t  y=$(y(s1, s, μ))")
          int = ξ_lensing(s1, s2(s1, s, μ), y(s1, s, μ);
               tol = tol, rtol = χ_rtol, atol = χ_atol, enhancer = enhancer)
          #println("int = $int")
          return int #.* (spline_F(s / s1, μ) * Pl(μ, L))
     else
          return (0.0, 0.0)
     end
end

function integral_on_mu_lensing(s1, s; pr::Bool = true, L::Integer = 0,
     enhancer = 1, tol = 0.5, χ_atol = 1e-3, χ_rtol = 1e-3, kwargs...)

     f(μ) = int_on_mu_lensing(s1, s, μ; enhancer = enhancer, tol = tol,
          L = L, χ_atol = χ_atol, χ_rtol = χ_rtol)[1]
     int = quadgk(μ -> f(μ), -0.99, 0.99; kwargs...)
     println("s1 = $s1 \t s2 = $s \t int = $int")
     return int
end


function map_integral_on_mu_lensing(s1 = s_eff;
     L::Integer = 0, pr::Bool = true,
     χ_atol = 1e-3, χ_rtol = 1e-3,
     enhancer = 1e6, tol = 0.5,
     kwargs...)

     t1 = time()
     ss = 10 .^ range(-1, 3, length = 100)
     #ss = range(tol, 1000, length = 1000)
     f(s) = integral_on_mu_lensing(s1, s; pr = pr, L = L, enhancer = enhancer,
          χ_atol = χ_atol, χ_rtol = χ_rtol, tol = tol, kwargs...)
     vec = @showprogress [f(s) ./ enhancer for s in ss]
     xis, xis_err = [x[1] for x in vec], [x[2] for x in vec]
     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu_lensing [in s] = $(t2-t1)\n")
     return (ss, xis, xis_err)
end


function PS_lensing(; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool = true, kwargs...)

     if ϕ(s_eff) > 0
          t1 = time()
          ks, pks = xicalc(s -> 2 * π^2 * integral_on_mu_lensing(s_eff, s; pr = pr, L = L, kwargs...)[1], L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
          t2 = time()
          pr && println("\ntime needed for PS_lensing [in s] = $(t2-t1)\n")

          if iseven(L)
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-1)^(L / 2)) .* pks
          else
               return ks, ((2 * L + 1) / A(s_min, s_max, θ_MAX) * ϕ(s_eff) * (-im)^L) .* pks
          end
     else
          throw(ErrorException(
               "ϕ(s_eff) should be >0, but in s_eff=$s_eff " *
               "is ϕ(s_eff) = $(ϕ(s_eff))! \n"
          ))
     end
end



##########################################################################################92



function print_map_int_on_mu_lensing(out::String; L::Integer = 0,
     s1 = s_eff, pr::Bool = true, kwargs...)

     t1 = time()
     vec = map_integral_on_mu_lensing(s1; L = L, pr = pr, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of xi_lensing.")
          parameters_used(io)
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
          print(io, "# kwards passed: ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for (i, key) in enumerate(keys(kwargs))
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          println(io, "\ns \t xi \t xi_error")
          for (s, xi, xi_err) in zip(vec[1], vec[2], vec[3])
               println(io, "$s \t $xi \t $(xi_err)")
          end
     end
end
