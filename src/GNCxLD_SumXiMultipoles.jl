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
     sum_ξ_GNCxLD_multipole(s1, s, cosmo::Cosmology; 
          kwargs...) ::Tuple{Float64, Vector{Float64}}

Evaluate the multipole of order `L` of all the GR effects TPCF multipoles and their sum
in `s1` and a distance `s` from it for the input Cosmology `cosmo`.

It makes a for-loop on the `GaPSE.GR_EFFECTS_GNCxLD` strings, calling `ξ_GNCxLD_multipole` for
each of them. They are currently:

`$(string(GaPSE.GR_EFFECTS_GNCxLD .* " , "...))`

## Optional arguments

- `kwards...` : all these keyword arguments will be passed to `ξ_GNCxLD_multipole`

## Returns

A tuple containing:
- the sum of all the ξ multipoles as first element
- a `Vector{Float64}` with all the values of each ξ; they are ordered
  following `GR_EFFECTS_GNCxLD`


See also: [`integrand_ξ_GNCxLD_multipole`](@ref), [`ξ_GNCxLD_multipole`](@ref),
[`map_sum_ξ_GNCxLD_multipole`](@ref), [`print_map_sum_ξ_GNCxLD_multipole`](@ref),
[`Cosmology`](@ref), [`GR_EFFECTS_GNCxLD`](@ref)
"""
function sum_ξ_GNCxLD_multipole(s1, s, cosmo::Cosmology; kwargs...)
     ALL = [ξ_GNCxLD_multipole(s1, s, effect, cosmo; specif_kwargs_GNCxLD(effect, kwargs)...)
            for effect in GaPSE.GR_EFFECTS_GNCxLD]

     return sum(ALL), ALL
end



##########################################################################################92



"""
     map_sum_ξ_GNCxLD_multipole(
          cosmo::Cosmology,
          v_ss = nothing;
          s1 = nothing,
          N_log::Int = 1000,
          kwargs...) ::Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{Float64}}}

Evaluate the multipole of order `L` of all the GR effects TPCF multipoles and their sum
in `s1` , for all the `s` values stored inside `v_ss` and for the input Cosmology `cosmo`.
If `v_ss = nothing`, it is set `v_ss = 10 .^ range(0, 3, length = N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.

It makes a for-loop on the `GaPSE.GR_EFFECTS_GNCxLD` strings, calling `map_ξ_GNCxLD_multipole` for
each of them. They are currently:

`$(string(GaPSE.GR_EFFECTS_GNCxLD .* " , "...))`

## Optional arguments

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `Cosmology`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `v_ss`, i.e. `range(0, 3, N_log)`; it is ignored if `v_ss ≠ nothing` 

- `kwards...` : all these keyword arguments will be passed to `map_ξ_GNCxLD_multipole`

## Returns

A tuple containing:
- the vector `v_ss` itself as first element;
- the  `Vector{Float64}` of the sum of all the ξ multipoles as second one
- a `Vector{Vector{Float64}}` with all the values of each ξ; they are ordered
  following `GR_EFFECTS_GNCxLD`


See also: [`map_ξ_GNCxLD_multipole`](@ref),
[`sum_ξ_GNCxLD_multipole`](@ref), [`print_map_sum_ξ_GNCxLD_multipole`](@ref),
[`Cosmology`](@ref), [`GR_EFFECTS_GNCxLD`](@ref)
"""
function map_sum_ξ_GNCxLD_multipole(
     cosmo::Cosmology,
     v_ss = nothing;
     N_log::Int = 1000,
     sum_xi::Bool = true,
     kwargs...)

     ss = isnothing(v_ss) ? 10 .^ range(0, 3, length = N_log) : v_ss

     ALL = [
          begin
               _, xis = map_ξ_GNCxLD_multipole(cosmo, effect, ss; 
                    sum_xi = sum_xi, specif_kwargs_GNCxLD(effect, kwargs)...)
               xis
          end for effect in GaPSE.GR_EFFECTS_GNCxLD
     ]

     return ss, sum(ALL), ALL
end



##########################################################################################92



"""
     print_map_sum_ξ_GNCxLD_multipole(
          cosmo::Cosmology,
          out::String,
          v_ss = nothing;
          s1 = nothing,
          N_log::Int = 1000,
          L::Int = 0,
          single::Bool = true,
          kwargs...) 

Evaluate the multipole of order `L` of all the GR effects TPCF multipoles and their sum
in `s1`, for all the `s` values stored inside `v_ss` and for the input Cosmology `cosmo`; 
saves the results inside the file `out`.
If `v_ss = nothing`, it is set `v_ss = 10 .^ range(0, 3, length = N_log)`.
If `s1 = nothing`, it is set `s1 = cosmo.s_eff`.

It makes a for-loop on the `GaPSE.GR_EFFECTS_GNCxLD` strings, calling `map_ξ_GNCxLD_multipole` for
each of them. They are currently:

`$(string(GaPSE.GR_EFFECTS_GNCxLD .* " , "...))`

## Optional arguments

- `s1 = nothing` : comoving distance from the observer where the TPCF should be evaluated;
  if `s1 = nothing`, it is automatically set `s1 = cosmo.s_eff` from the given input `Cosmology`.

- `N_log::Int = 1000` : number of points to be used in the default logaritmically-spaced 
  range for `v_ss`, i.e. `range(0, 3, N_log)`; it is ignored if `v_ss ≠ nothing` 

- `L::Int = 0`: order of the Legendre polynomial to be used

- `single::Bool = true` : if `true`, all the CFs are printed inside the file of the sum, in a 
  table with 18 columns (first one for `v_ss`, second for their sum an the next 16 for each effect).
  Otherwise, a new directory "all_standalones_CFs" is created (in the same path given in `out`) and 
  they are separately saved in files there placed.

- `kwards...` : all these keyword arguments will be passed to `map_ξ_GNCxLD_multipole`


See also: [`map_ξ_GNCxLD_multipole`](@ref),
[`sum_ξ_GNCxLD_multipole`](@ref), [`map_sum_ξ_GNCxLD_multipole`](@ref),
[`Cosmology`](@ref), [`GR_EFFECTS_GNCxLD`](@ref)
"""
function print_map_sum_ξ_GNCxLD_multipole(
     cosmo::Cosmology,
     out::String,
     v_ss = nothing;
     s1 = nothing,
     L::Int = 0,
     single::Bool = true,
     pr::Bool = true,
     sum_xi::Bool = true,
     kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     dir = length(split(out, "/")) == 1 ? "all_GNCxLD_standalones_CFs/" :
           join(split(out, "/")[begin:end-1] .* "/") * "all_GNCxLD_standalones_CFs/"

     if single == false
          isdir(dir) || mkdir(dir)
     end

     s_1 = isnothing(s1) ? cosmo.s_eff : s1
     t1 = time()
     ss, xis, ALL = map_sum_ξ_GNCxLD_multipole(cosmo, v_ss;
          L = L, s1 = s_1, pr = pr, sum_xi = sum_xi, kwargs...)
     t2 = time()

     (pr) && println("\ntime needed for map_sum_ξ_GNCxLD_multipole " *
                   "[in s] = $(@sprintf("%.5f", t2-t1)) \n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "#\n# This is an integration map on mu of the sum" *
                      " of all the ξ_GNCxLD_L=$L multipole GR effects.")
          !(single == true) ||
               println(io, "# In input was set \"single = $single\", " *
                           "so, together with their sum, all the CFs are here reported.\n#")
          !(single == false) ||
               println(io, "# In input was set \"single = $single\", so here is showed only \n" *
                           "# the sum of all the CFs, and them are saved separately in the \n" *
                           "# following directory: $dir \n#")

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
          isnothing(s1) ||
               println(io, "#\n# NOTE: the computation is done not in s1 = s_eff, \n" *
                           "#\t because you specified in input s1 = $s1 Mpc/h_0!")
          println(io, "# ")

          if (single == false)
               println(io, "# s [Mpc/h_0] \t \t xi_SUM")
               for (s, xi) in zip(ss, xis)
                    println(io, "$s \t $xi")
               end

          elseif single == true
               println(io, "# 1: s [Mpc/h_0] \t 2: xi_SUM \t " *
                    join([string(i) for i in 3:length(GR_EFFECTS_GNCxLD)+2] .* 
                         ": xi_" .* GR_EFFECTS_GNCxLD .* " \t "))
               for (i, (s, xi)) in enumerate(zip(ss, xis))
                    println(io, "$s \t $xi \t " *
                                join(["$(v[i]) \t " for v in ALL]))
               end
          else
               throw(ErrorException("how did you arrive here?"))
          end
     end

     if single == false
          for (effect, vec) in zip(GaPSE.GR_EFFECTS_GNCxLD, ALL)
               open(dir * "xi_GNCxLD_" * effect * "_L$L" * ".txt", "w") do io
                    println(io, "# This is an integration map on mu of the ξ multipole $effect GR effect.")
                    parameters_used(io, cosmo)
                    #println(io, "# computational time needed (in s) : $(@sprintf("%.4f", t2-t1))")
                    print(io, "# kwards passed: ")

                    if isempty(kwargs)
                         println(io, "none")
                    else
                         print(io, "\n")
                         for key in keys(kwargs)
                              println(io, "# \t\t$(key) = $(kwargs[key])")
                         end
                    end
                    isnothing(s1) || println(io, "#\n# NOTE: the computation is done not in " *
                                                  "s1 = s_eff, because you specified in input s1 = $s1 !")
                    println(io, "# ")
                    println(io, "# s [Mpc/h_0] \t \t xi")
                    for (s, xi) in zip(ss, vec)
                         println(io, "$s \t $xi")
                    end

               end
          end
     end
end

