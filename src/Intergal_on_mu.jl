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


function integral_on_mu(s1, s, integrand::Function;
     L::Integer = 0,
     enhancer::Float64 = 1e6,
     use_windows::Bool = true,
     μ_atol::Float64 = 1e-3,
     μ_rtol::Float64 = 1e-3,
     kwargs...
)

     f(μ) = integrand(s1, s, μ; enhancer = enhancer, L = L,
          use_windows = use_windows, kwargs...)[1]

     #println("s1 = $s1 \t s = $s")
     int = quadgk(μ -> f(μ), -1.0, 1.0; rtol = μ_rtol, atol = μ_atol)
     #println("s1 = $s1 \t s2 = $s \t int = $int")
     return int ./ enhancer
end


function integral_on_mu(s1, s, name::String; kwargs...)
     error = "$name is not a valid GR effect name.\n" *
             "Valid GR effect names are the following:\n" *
             string(keys(dict_gr_mu) .* " , "...)

     @assert (name ∈ keys(dict_gr_mu)) error
     integral_on_mu(s1, s, dict_gr_mu[name]; kwargs...)
end


function map_integral_on_mu(
     effect::Union{String, Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing,
     s1::Float64 = s_eff; 
     pr::Bool = true, enhancer = 1e6, kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(-1, 3, length = 100) : v_ss
     f(s) = integral_on_mu(s1, s, effect; enhancer = enhancer, kwargs...)
     vec = @showprogress [f(s) for s in ss]
     xis, xis_err = [x[1] for x in vec], [x[2] for x in vec]
     t2 = time()
     pr && println("\ntime needed for map_integral_on_mu for $effect [in s] = $(t2-t1)")
     return (ss, xis, xis_err)
end


function print_map_int_on_mu(out::String,
     effect::Union{String, Function},
     v_ss::Union{Vector{Float64},Nothing} = nothing,
     s1::Float64 = s_eff;
     kwargs...)

     t1 = time()
     vec = map_integral_on_mu(effect, v_ss, s1; kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# This is an integration map on mu of the $effect GR effect.")
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


