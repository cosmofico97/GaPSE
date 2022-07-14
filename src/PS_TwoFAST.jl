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



function TwoFAST_PS_multipole(f_in;
     int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
     L::Int=0, N::Int=1024, pr::Bool=true,
     k0::Union{Nothing,Float64}=nothing,
     right::Union{Float64,Nothing}=nothing)

     k_0 = isnothing(k0) ? 1.0 / int_s_max : k0

     t1 = time()
     ks, pks = xicalc(s -> 2 * π^2 * f_in(s), L, 0;
          N=N, kmin=int_s_min, kmax=int_s_max, r0=k_0)
     t2 = time()
     pr && println("\ntime needed for this Power Spectrum computation [in s] = $(t2-t1)\n")

     if isnothing(right)
          if iseven(L)
               return ks, (1 / A_prime * (-1)^(L / 2)) .* pks
          else
               return ks, (1 / A_prime * (-im)^L) .* pks
          end
     else
          if iseven(L)
               return ks[ks.<right], (1 / A_prime * (-1)^(L / 2)) .* pks[ks.<right]
          else
               return ks[ks.<right], (1 / A_prime * (-im)^L) .* pks[ks.<right]
          end
     end
end


function TwoFAST_PS_multipole(SS, FS;
     int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
     epl::Bool=true, pr::Bool=true, L::Int=0,
     N_left::Int=12, N_right::Int=12,
     p0_left=[-2.0, 1.0], p0_right=[-2.0, 1.0],
     cut_first_n::Int = 0, cut_last_n::Int=0)

     @assert length(SS) == length(FS) "length(SS) == length(FS) must hold!"
     @assert cut_first_n ≥ 0 "cut_first_n ≥ 0 must hold!"
     @assert cut_last_n ≥ 0 "cut_last_n ≥ 0 must hold!"
     @assert cut_first_n + cut_last_n < length(FS) "cut_first_n + cut_last_n < length(ss) must hold!"

     a, b = 1 + cut_first_n, length(SS) - cut_last_n
     ss, fs = SS[a:b], FS[a:b]

     N = length(ss)
     k_0 = 1.0 / max(ss...)
     right = 1.0 / min(ss...)

     f_in, INT_s_min, INT_s_max =
          if epl == true
               if all(fs[begin:begin+5] .≈ 0.0) || all(fs[end-5:end] .≈ 0.0)
                    spl = Spline1D(ss, fs; bc="error")
                    f(x) = ((x ≤ ss[1]) || (x ≥ ss[end])) ? 0.0 : spl(x)
                    f, int_s_min, int_s_max
               else
                    EPLs(ss, fs, p0_left, p0_right;
                         N_left=N_left, N_right=N_right), int_s_min, int_s_max
               end
          else
               Spline1D(ss, fs; bc="error"), min(ss...), max(ss...)
          end

     return TwoFAST_PS_multipole(f_in; int_s_min=INT_s_min, int_s_max=INT_s_max,
          k0=k_0, right=right, N=N, pr=pr, L=L)
end



function TwoFAST_all_PS_multipole(input::String,
          group::String=VALID_GROUPS[end];
          L::Int=0, pr::Bool=true, kwargs...)

     check_group(group; valid_groups=VALID_GROUPS)
     check_fileisingroup(input, group)

     pr && begin
          print("\nI'm computing the PS_multipole from the file \"$input\"")
          if group == "GNC"
               println("for the Galaxy Number Counts.")
          elseif group == "LD"
               println("for the Luminosity Distance perturbations.")
          elseif group == "GNCxLD"
               println("for the cross correlations between " *
                       "Galaxy Number Counts and Luminosity Distance perturbations.")
          elseif group == "LDxGNC"
               println("for the cross correlations between " *
                       "Luminosity Distance perturbations and Galaxy Number Counts.")
          else
               println("(no specific group considered).")
          end
     end

     sps = group == "GNC" ? "GNC GR effects" :
           group == "LD" ? "LD GR effects" :
           group == "GNCxLD" ? "GNCxLD GR effects" :
           group == "LDxGNC" ? "LDxGNC GR effects" :
           "generic file"

     ks, VEC = begin
          table = readdlm(input; comments=true)
          xs = convert(Vector{Float64}, table[:, 1])
          all_YS = [convert(Vector{Float64}, col)
                    for col in eachcol(table[:, 2:end])]
          res = pr ? begin
               @showprogress sps * ", L=$L: " [
                    TwoFAST_PS_multipole(xs, ys; L=L, pr=false, kwargs...)
                    for ys in all_YS]
          end : begin
               [TwoFAST_PS_multipole(xs, ys; L=L, pr=false, kwargs...)
                for ys in all_YS]
          end

          res[1][1], [res[i][2] for i in 1:length(res)]
     end

     return ks, VEC
end
