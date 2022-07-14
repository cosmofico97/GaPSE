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
# along with  If not, see <http://www.gnu.org/licenses/>.
#


struct TF
     left_value::Float64
     left::Float64

     spline::Dierckx.Spline1D

     r_si::Float64
     r_b::Float64
     r_a::Float64
     right::Float64


     function TF(ks, Tks)
          left_val = sum(Tks[1:10]) / 10
          r_si, r_b, r_a = power_law_from_data(
               ks, Tks, [-2.0, 1.0], ks[end-15], ks[end]; con=false)
          spline = Spline1D(ks, Tks; bc="error")

          new(left_val, ks[2], spline, r_si, r_b, r_a, ks[end])
     end

end;


function (tf::TF)(k)
     if k < tf.left
          return tf.left_value
     elseif k > tf.right
          return power_law(k, tf.r_si, tf.r_b, tf.r_a)
     else
          return tf.spline(k)
     end
end;


##########################################################################################92


function α_bias(k, tf::TF; bf=1.0, D=1.0, Ω_M0=0.29992)
     return 1.5 * bf * Ω_M0 * (100 / 299792.458)^2 / (0.779017 * D * k^2 * tf(k))
end


struct IntegralIPSalpha
     l_si::Float64
     l_b::Float64
     l_a::Float64
     left::Float64

     spline::Dierckx.Spline1D

     r_si::Float64
     r_b::Float64
     r_a::Float64
     right::Float64

     function IntegralIPSalpha(tf::TF, cosmo::Cosmology, l, n=0;
          D=nothing, bf=1.0,
          N::Int=1024, kmin=1e-6, kmax=1e4, s0=1e-4,
          fit_left_min=nothing, fit_left_max=nothing, p0_left=nothing,
          fit_right_min=nothing, fit_right_max=nothing, p0_right=nothing)

          DD = isnothing(D) ? cosmo.D_of_s(cosmo.s_eff) : D
          Ω_M00 = cosmo.params.Ω_M0

          rs, xis = xicalc(k -> cosmo.IPS(k) * α_bias(k, tf; bf=bf, D=DD, Ω_M0=Ω_M00), l, n;
               N=N, kmin=kmin, kmax=kmax, r0=s0)

          fit_left_MIN = !isnothing(fit_left_min) ? fit_left_min : begin
               l ≈ 0.0 ? 5e-2 : l ≈ 2.0 ? 5e-1 : rs[2]
          end
          fit_left_MAX = !isnothing(fit_left_max) ? fit_left_max : begin
               l ≈ 0.0 ? 1e-1 : l ≈ 2.0 ? 1e0 : rs[16]
          end
          p_0_left = isnothing(p0_left) ? [-1.0, 1.0] : p0_left
          l_si, l_b, l_a = power_law_from_data(
               rs, xis, p_0_left, fit_left_MIN, fit_left_MAX; con=false)

          fit_right_MIN = isnothing(fit_right_min) ? rs[length(rs)-16] : fit_right_min
          fit_right_MAX = isnothing(fit_right_max) ? rs[length(rs)-1] : fit_right_max
          p_0_right = isnothing(p0_right) ? [-4.0, 1.0] : p0_right
          r_si, r_b, r_a = power_law_from_data(
               rs, xis, p_0_right, fit_right_MIN, fit_right_MAX; con=false)

          ind_left = findfirst(x -> x > fit_left_MIN, rs) - 1
          ind_right = findfirst(x -> x >= fit_right_MAX, rs)
          new_rs = vcat(rs[ind_left:ind_right])
          new_Js = vcat(xis[ind_left:ind_right])
          spline = Spline1D(new_rs, new_Js; bc="error")

          #println("\nleft = $l_si , $l_b , $l_a, $fit_left_min")
          #println("right = $r_si , $r_b , $r_a, $fit_right_MAX\n")

          new(l_si, l_b, l_a, fit_left_MIN, spline, r_si, r_b, r_a, fit_right_MAX)
     end
end;


function (Jl::IntegralIPSalpha)(x)
     if x < Jl.left
          return power_law(x, Jl.l_si, Jl.l_b, Jl.l_a)
     elseif x > Jl.right
          #warning("i am going too right! ")
          return power_law(x, Jl.r_si, Jl.r_b, Jl.r_a)
     else
          return Jl.spline(x)
     end
end;



##########################################################################################92



struct CosmoPNG
     tf::TF
     file_TF::String

     J0::IntegralIPSalpha
     J2::IntegralIPSalpha

     params::Dict{Symbol,T1} where {T1}

     function CosmoPNG(
               cosmo::Cosmology,
               file_TF::String;
               comments::Bool = true, D = nothing, bf = 1.0,
               flm0 = 5e-2, flM0 = 1e-1, flm2 = 5e-1, flM2 = 1e0,
               kmin = 1e-6, kmax=1e4, N::Int = 1024
          )

          table = readdlm(file_TF; comments=comments)
          ks_tf = convert(Vector{Float64}, table[:, 1])
          pks_tf = convert(Vector{Float64}, table[:, 2])

          DD = isnothing(D) ? cosmo.D_of_s(cosmo.s_eff) : D
          tf = TF(ks_tf, pks_tf)
     
          J0 = IntegralIPSalpha(tf, cosmo, 0, 0; D=DD, bf=bf, kmin=kmin, kmax=kmax,
               fit_left_min=flm0, fit_left_max=flM0, N=N)
          J2 = IntegralIPSalpha(tf, cosmo, 2, 0; D=DD, bf=bf, kmin=kmin, kmax=kmax,
               fit_left_min=flm2, fit_left_max=flM2, N=N)
     
          params = Dict(:D=>D, :bf=>bf, :flm0=>flm0, :flM0=>flM0, 
               :flm2=>flm2, :flM2=>flM2, :kmin=>kmin, :kmax=>kmax, :N=>N)

          new(tf, file_TF, J0, J2, params)
     end
end


##########################################################################################92



function ξ_S_L0(P::Point, cosmo::Cosmology, cosmopng::CosmoPNG)
     b = cosmo.params.b
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     2.0 * (b + f / 3.0) * D^2 * cosmopng.J0(s)
end

function ξ_S_L0(s1, cosmo::Cosmology, cosmopng::CosmoPNG)
     P1 = Point(s1, cosmo)
     return ξ_S_L0(P1, cosmo, cosmopng)
end

function ξ_S_L2(P::Point, cosmo::Cosmology, cosmopng::CosmoPNG)
     s = P.comdist

     Peff = Point(cosmo.s_eff, cosmo)
     D, f = Peff.D, Peff.f

     -4.0 / 3.0 * f * D^2 * cosmopng.J2(s)
end

function ξ_S_L2(s1, cosmo::Cosmology, cosmopng::CosmoPNG)
     P1 = Point(s1, cosmo)
     return ξ_S_L2(P1, cosmo, cosmopng)
end


function ξ_S(s1, y, cosmo::Cosmology, cosmopng::CosmoPNG)
     ξ_S_L0(s1, cosmo, cosmopng) + ξ_S_L2(s1, cosmo, cosmopng) * Pl(y, 2)
end



##########################################################################################92



function integrand_ξ_S_multipole(s, μ, cosmo::Cosmology, cosmopng::CosmoPNG;
     L::Int=0, use_windows::Bool=true)

     res = if use_windows == true
          ξ_S(s, μ, cosmo, cosmopng) .* (
               spline_integrF(s, μ, cosmo.windowFint) / cosmo.WFI_norm * Pl(μ, L)
               )
     else
          ξ_S(s, μ, cosmo, cosmopng) .* Pl(μ, L)
     end

     return (2.0 * L + 1.0) / 2.0 * res
end

function ξ_S_multipole(
     s, cosmo::Cosmology, cosmopng::CosmoPNG;
     L::Int=0,
     use_windows::Bool=true,
     enhancer::Float64=1e6,
     μ_atol::Float64=0.0,
     μ_rtol::Float64=1e-2)

     orig_f(μ) = enhancer * integrand_ξ_S_multipole(s, μ, cosmo, cosmopng;
          L=L, use_windows=use_windows)

     int = quadgk(μ -> orig_f(μ), -1.0, 1.0; atol=μ_atol, rtol=μ_rtol)[1]

     return int / enhancer
end

function map_ξ_S_multipole(cosmo::Cosmology, cosmopng::CosmoPNG,
     v_ss=nothing;
     pr::Bool=true,
     N_log::Int=1000,
     L::Int=0,
     kwargs...)

     t1 = time()
     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length=N_log) : v_ss
     xis = pr ? begin
          @showprogress "ξ_S, L=$L: " [
               ξ_S_multipole(s, cosmo, cosmopng; L=L, kwargs...) for s in ss
          ]
     end : [
          ξ_S_multipole(s, cosmo, cosmopng; L=L, kwargs...) for s in ss
     ]

     t2 = time()
     pr && println("\ntime needed for map_ξ_S_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) ")
     return (ss, xis)
end


function print_map_ξ_S_multipole(
     cosmo::Cosmology, cosmopng::CosmoPNG,
     out::String,
     v_ss=nothing;
     kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     t1 = time()
     vec = map_ξ_S_multipole(cosmo, cosmopng, v_ss; kwargs...)
     t2 = time()

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "#\n# This is an integration map on mu of the ξ multipole of S,\n" *
                      "# which is the difference between the Power Spectrum with PNG and with f_NL = 0.")
          println(io, "# Transfer function read from the file: $(cosmopng.file_TF)")
          print(io, "# Parameters used for the considered CosmoPNG: ")
          print(io, "\n")
          for key in keys(cosmopng.params)
               println(io, "# \t\t$(key) = $(cosmopng.params[key])")
          end
          println(io, "#")

          parameters_used(io, cosmo; logo=false)
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


