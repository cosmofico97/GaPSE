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


function FFTLog_PS_multipole(ss, xis;
     pr::Bool=true,
     L::Int=0, ν::Union{Float64,Nothing}=nothing,
     n_extrap_low::Int=500,
     n_extrap_high::Int=500, n_pad::Int=500,
     cut_first_n::Int=0, cut_last_n::Int=0)

     @assert length(ss) == length(xis) "length(ss) == length(xis) must hold!"
     @assert cut_first_n ≥ 0 "cut_first_n ≥ 0 must hold!"
     @assert cut_last_n ≥ 0 "cut_last_n ≥ 0 must hold!"
     @assert cut_first_n + cut_last_n < length(ss) "cut_first_n + cut_last_n < length(ss) must hold!"

     a, b = 1 + cut_first_n, length(ss) - cut_last_n

     SS = ss[a:b]
     XIS = xis[a:b] .* SS .^ 3

     new_ν = if isnothing(ν)
          #l_si, l_b, l_a = GaPSE.power_law_from_data(SS, XIS, [1.0, 1.0], SS[begin], SS[begin+10])
          #r_si, r_b, r_a = GaPSE.power_law_from_data(SS, XIS, [1.0, 1.0], SS[end-10], SS[end])
          #1.9 + (l_si + r_si)/2
          1.5
     else
          ν
     end
     #NEW_SS, NEW_XIS = 

     t1 = time()
     plan = SingleBesselPlan(x=SS,
          n_extrap_low=n_extrap_low, ν=new_ν,
          n_extrap_high=n_extrap_high, n_pad=n_pad)
     prepare_FFTLog!(plan, [L])
     pks = reshape(evaluate_FFTLog(plan, XIS), (:,))
     ks = reshape(get_y(plan), (:,))
     t2 = time()
     pr && println("\ntime needed for this Power Spectrum computation [in s] = $(t2-t1)\n")


     if iseven(L)
          return ks, (1 / A_prime * (-1)^(L / 2)) .* pks
     else
          return ks, (1 / A_prime * (-im)^L) .* pks
     end
end



function FFTLog_all_PS_multipole(input::String,
     group::String=VALID_GROUPS[end];
     L::Int=0, pr::Bool=true,
     ν::Union{Float64,Nothing,Vector{Float64}}=nothing,
     n_extrap_low::Int=500,
     n_extrap_high::Int=500, n_pad::Int=500,
     cut_first_n::Int=0, cut_last_n::Int=0)

     @assert cut_first_n ≥ 0 "cut_first_n ≥ 0 must hold!"
     @assert cut_last_n ≥ 0 "cut_last_n ≥ 0 must hold!"


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

     table = readdlm(input; comments=true)
     ss = convert(Vector{Float64}, table[:, 1])
     all_xis = [convert(Vector{Float64}, col)
                for col in eachcol(table[:, 2:end])]

     used_νs = if typeof(ν) == Vector{Float64}
          @assert length(ν) == length(all_xis) "length(ν) == length(all_xis) must hold!"
          ν
     else
          [ν for i in 1:length(all_xis)]
     end

     a, b = 1 + cut_first_n, length(ss) - cut_last_n
     #=

     SS = ss[a:b]
     ALL_XIS = [xis[a:b] .* SS .^ 3 for xis in all_xis]


     plan = SingleBesselPlan(x=SS,
     n_extrap_low=n_extrap_low, ν=ν,
     n_extrap_high=n_extrap_high, n_pad=n_pad)
     prepare_FFTLog!(plan, [L])
     ALL_pks = pr ? begin
     @showprogress sps * ", L=$L: " [reshape(evaluate_FFTLog(plan, XIS), (:,))
     for XIS in ALL_XIS] 
     end : begin
     [reshape(evaluate_FFTLog(plan, XIS), (:,))
     for XIS in ALL_XIS] 
     end
     =#

     plan = SingleBesselPlan(x=ss[a:b],
          n_extrap_low=n_extrap_low, ν=1.1,
          n_extrap_high=n_extrap_high, n_pad=n_pad)
     prepare_FFTLog!(plan, [L])
     ALL_pks = pr ? begin
          @showprogress sps * ", L=$L: " [
               begin
                    _, pks = FFTLog_PS_multipole(ss, xis;
                         pr=false, L=L, ν=specific_ν,
                         n_extrap_low=n_extrap_low,
                         n_extrap_high=n_extrap_high, n_pad=n_pad,
                         cut_first_n=cut_first_n, cut_last_n=cut_last_n)
                    pks
               end
               for (xis, specific_ν) in zip(all_xis, used_νs)]
     end : begin
          [
               begin
                    _, pks = FFTLog_PS_multipole(ss, xis;
                         pr=false, L=L, ν=specific_ν,
                         n_extrap_low=n_extrap_low,
                         n_extrap_high=n_extrap_high, n_pad=n_pad,
                         cut_first_n=cut_first_n, cut_last_n=cut_last_n)
                    pks
               end
               for (xis, specific_ν) in zip(all_xis, used_νs)
          ]
     end

     ks = reshape(get_y(plan), (:,))

     return ks, ALL_pks


     #=
     if iseven(L)
          return ks, [(1 / A_prime * (-1)^(L / 2)) .* pks for pks in ALL_pks]
     else
          return ks, [(1 / A_prime * (-im)^L) .* pks for pks in ALL_pks]
     end
     =#
end

