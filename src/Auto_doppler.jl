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


function ξ_doppler(s1, s2, y; enhancer = 1, tol = 1)
     delta_s = s(s1, s2, y)
     (delta_s > tol) || (return 0.0)

     D1, D2 = D(s1), D(s2)

     f1, ℋ1 = f(s1), ℋ(s1)
     f2, ℋ2 = f(s2), ℋ(s2)
     #f1, ℋ1, ℋ1_p, s_b1 = f(s1), ℋ(s1), ℋ_p(s1), s_b(s1)
     #f2, ℋ2, ℋ2_p, s_b2 = f(s2), ℋ(s2), ℋ_p(s2), s_b(s2)
     #ℛ1, ℛ2 = ℛ(s1, ℋ1, ℋ1_p, s_b1), ℛ(s2, ℋ2, ℋ2_p, s_b2)
     ℛ1, ℛ2 = 1 - 1 / (ℋ1 * s1), 1 - 1 / (ℋ2 * s2)

     prefac = D1 * D2 * f1 * f2 * ℛ1 * ℛ2 * ℋ1 * ℋ2
     c1 = 3 * s1 * s2 - 2 * y * (s1^2 + s2^2) + s1 * s2 * y^2

     parenth = I00(delta_s) / 45.0 + I20(delta_s) / 31.5 + I40(delta_s) / 105.0

     first = prefac * (c1 * parenth + I02(delta_s) * y * delta_s^2 / 3.0)

     return 0.5 * enhancer * first
end


function integrand_on_mu_doppler(s1, s, μ; L::Integer = 0, enhancer = 1, tol = 1)
     (ϕ(s2(s1, s, μ)) > 0) || (return 0.0)
     val = ξ_doppler(s1, s2(s1, s, μ), y(s1, s, μ), enhancer = enhancer, tol = tol)
     return val * spline_F(s / s1, μ) * Pl(μ, L)
end


function integral_on_mu_doppler(s1, s; L::Integer = 0, enhancer = 1, tol = 1, kwargs...)
     f(μ) = integrand_on_mu_doppler(s1, s, μ; L = L, enhancer = enhancer, tol = tol)
     return quadgk(f, -1, 1; kwargs...)
end

function map_integral_on_mu_doppler(s1 = s_eff;
     L::Integer = 0, pr::Bool = true, enhancer = 1e6, tol = 1, kwargs...)

     t1 = time()
     ss = 10 .^ range(log10(tol), 3, length = 100)
     f(s) = integral_on_mu_doppler(s1, s; L = L, enhancer = enhancer, tol = tol, kwargs...)
     vec = @showprogress [f(s) ./ enhancer for s in ss]
     xis, xis_err = [x[1] for x in vec], [x[2] for x in vec]
     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu_doppler [in s] = $(t2-t1)\n")
     return (ss, xis, xis_err)
end




function print_map_int_on_mu_doppler(out::String; L::Integer = 0,
     s1 = s_eff, pr::Bool = true, kwargs...)

     t1 = time()
     vec = map_integral_on_mu_doppler(s1; L = L, pr = pr, kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of xi_doppler.")
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



##########################################################################################92



# mean time of evaluation: 141 seconds!
# too long to be used
function PS_doppler_exact(; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,
     L::Integer = 0, pr::Bool=true, kwargs...)

     if ϕ(s_eff) > 0
          t1 = time()
          ks, pks = xicalc(s -> 2 * π^2 * integral_on_mu_doppler(s_eff, s; L=L, kwargs...)[1], L, 0;
               N = N, kmin = int_s_min, kmax = int_s_max, r0 = 1 / int_s_max)
          t2 = time()
          pr && println("\ntime needed for PS_doppler [in s] = $(t2-t1)\n")

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

