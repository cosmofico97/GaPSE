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

function PS_multipole(ss, fs; twofast::Bool = false, 
     cut_first_n::Int = 0, cut_last_n::Int=0, kwargs...)

     @assert length(ss) == length(fs) "length(ss) == length(fs) must hold!"
     @assert cut_first_n ≥ 0 "cut_first_n ≥ 0 must hold!"
     @assert cut_last_n ≥ 0 "cut_last_n ≥ 0 must hold!"
     @assert cut_first_n + cut_last_n < length(ss) "cut_first_n + cut_last_n < length(ss) must hold!"

     a, b = 1 + cut_first_n, length(ss) - cut_last_n

     if twofast == true
          return TwoFAST_PS_multipole(ss[a:b], fs[a:b]; kwargs...)
     else
          return FFTLog_PS_multipole(ss[a:b], fs[a:b]; kwargs...)
     end
end


function PS_multipole(input::String; kwargs...)

     xi_table = readdlm(input, comments=true)
     ss = convert(Vector{Float64}, xi_table[:, 1])
     fs = convert(Vector{Float64}, xi_table[:, 2])

     return PS_multipole(ss, fs; kwargs...)
end




"""
     PS_multipole(input::String; kwargs...
          ) ::Tuple{Vector{Float64}, Vector{Float64}}

     PS_multipole(ss, fs;
          int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
          epl::Bool=true,
          N_left::Int=12, N_right::Int=12,
          p0_left=[-2.0, 1.0], p0_right=[-2.0, 1.0],
          kwargs...)

     PS_multipole(f_in;
          int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
          L::Int=0, pr::Bool=true, N::Int = 1024,
          k0::Union{Nothing,Float64}=nothing,
          right::Union{Float64,Nothing}=nothing)

Return the `L`-order PS multipole through the
following Fast Fourier Transform and the effective redshift approximation:

```math
P_L(k) = \\frac{2 L + 1}{A^{'}} (-i)^L \\, \\phi(s_\\mathrm{eff}) \\int_0^\\infty 
        \\mathrm{d} s \\; s^2 \\, j_L(ks) \\, f_\\mathrm{in}(s) \\; ,
        \\quad \\; A^{'} = \\frac{1}{4\\,\\pi}
```

The former method takes the name of the file where the
TPCF multipole in exam is saved in, opens that file and pass the `xs` and  `ys` vector
of values to the second method. All the keyword- arguments given to the first method are directly
transferred to the second one.

In the second method, you have to pass the `xs` and `ys` values of the TPCF you want to
exam. Internally, if `epl==true` the data are fitted with two power laws at the edges, creating a
`EPLs` struct that is passed to the third method.
The great advantage is that the integration can be extended over the limits imposed by the vector themself,
increasing by far the precision of the PS computation.

In the last method, you have to give its function/spline `f_in`.

It's recommended to use either the first or the second method.

The analytical expression previously showed can be easily obtained from the 
standard one:
```math
\\begin{split}
    P_L(k) = &\\frac{2 L + 1}{A} (-i)^L \\, 
        \\int_0^\\infty \\mathrm{d} s_1 \\; s_1^2 
        \\int_0^\\infty \\mathrm{d} s \\; s^2 
        \\int_{-1}^{+1} \\mathrm{d} \\mu \\;
        j_L(ks) \\, \\xi(s_1, s, \\mu) \\, \\phi(s_1) \\, \\phi(s_2) \\,
        \\mathcal{L}_L(\\mu) F\\left(\\frac{s}{s_1}, \\mu \\right) \\\\
        &\\mathrm{with} \\; \\;s_2 = s_2(s_1, s, μ) = \\sqrt{s_1^2 + s^2 + 2s_1s\\mu}
        \\; 
        , \\quad A(s_\\mathrm{max}, s_\\mathrm{min}, \\theta_\\mathrm{max}) 
        \\frac{
          V(s_\\mathrm{max}, s_\\mathrm{min}, \\theta_\\mathrm{max})
          }{4 \\, \\pi^2}
\\end{split}
```
with the definition
```math
f_\\mathrm{in}(s_1, s) =  \\int_{-1}^{+1} \\mathrm{d} \\mu \\;
        \\xi(s_1, s, \\mu) \\, \\phi(s_2) \\,
        \\mathcal{L}_L(\\mu) \\, F\\left(\\frac{s}{s_1}, \\mu \\right)
```
and the application of the effective redshift approximation.

The computation is made through the `xicalc` function of the 
[TwoFAST](https://github.com/hsgg/TwoFAST.jl) Julia package.
Note that in the computation the integration range ``0\\leq s \\leq \\infty`` 
is reduced to `int_s_min ≤ s ≤ int_s_max``

## Optional arguments

- `int_s_min::Float64 = 1e-1` and `int_s_max::Float64 = 1e3` : extremes of
  integration. Do not worry if the input TPCF is not defined in all the 
  integration range: it will be fitted with two power laws at its extremes, 
  thanks to the function `EPLs`

- `N_left::Int = 12` and `N_right::Int = 12` : number of points from left
  right edges to be used for the power law fitting in `EPLs`. They matters only
  if in the given input file ξ is not defined until the extremes of integration
  `int_s_min` and `int_s_max`

- `L::Int = 0`: order of the Legendre polynomial to be used; note that 
  the multipole order `L` must be the same of the input TPCF in exam! 
  Otherwise, the results would have no sense at all!

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `N::Int = 1024` : number of points to be returned by `xicalc`


## Returns

A `Tuple{Vector{Float64}, Vector{Float64}}` with:
- the `k` values vector as first element;
- the correspoding PS `pk` values vector as second one.


See also: [`V_survey`](@ref), [`A`](@ref), [`A_prime`](@ref),
[`EPLs`](@ref),  [`print_PS_multipole`](@ref)
"""
PS_multipole



##########################################################################################92


function print_PS_multipole(input::String, out::String;
     L::Int=0, pr::Bool=true, kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     pr && println("""\nI'm computing the PS_multipole from the file "$input" """)

     time_1 = time()
     vec = PS_multipole(input; L=L, pr=pr, kwargs...)
     time_2 = time()

     N = length(vec[1])

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "# Power Spectrum Multipole computation of the file: \"$input\"")
          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
          print(io, "# kwards passed to \"print_PS_multipole\": ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          println(io, "# ")
          println(io, "# k [h_0/Mpc] \t \t  P [(Mpc/h_0)^3]")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t $pk")
          end
     end
end

function print_PS_multipole(ss, fs, out::String;
     L::Int=0, pr::Bool=true, kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     pr && println("""\nI'm computing the PS_multipole from the two input vectors.""")

     time_1 = time()
     vec = PS_multipole(ss, fs; L=L, pr=pr, kwargs...)
     time_2 = time()

     N = length(vec[1])

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          println(io, "# Power Spectrum Multipole computation from two input vectors.")
          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
          print(io, "# kwards passed to \"print_PS_multipole\": ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          println(io, "# ")
          println(io, "# k [h_0/Mpc] \t \t  P [(Mpc/h_0)^3]")
          for (k, pk) in zip(vec[1], vec[2])
               println(io, "$k \t " * GaPSE.number_to_string(pk))
          end
     end
end


"""
     print_PS_multipole(input::String, out::String;
          L::Int = 0, N::Int = 1024,
          pr::Bool = true, kwargs...)

Takes in input a filename `input` where is stored a TPCF multipole,
calculate the `L`-order PS multipole through the
following Fast Fourier Transform and the effective redshift approximation

```math
P_L(k) = \\frac{2 L + 1}{A^{'}} (-i)^L \\, \\phi(s_\\mathrm{eff}) \\int_0^\\infty 
        \\mathrm{d} s \\; s^2 \\, j_L(ks) \\, f_\\mathrm{in}(s) \\; ,
        \\quad \\; A^{'} = \\frac{1}{4\\,\\pi}
```

and save it in the file `out`, together with the options used for the computation.

## Optional arguments

- `L::Int = 0`: order of the Legendre polynomial to be used; note that 
  the multipole order `L` must be the same of the input TPCF in exam! 
  Otherwise, the results would have no sense at all!

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `N::Int = 1024` : number of points to be returned by `xicalc`

- `kwargs...` : other keyword arguments that will be passed to `PS_multipole`

See also: [`V_survey`](@ref), [`A`](@ref), [`A_prime`](@ref),
[`EPLs`](@ref), [`PS_multipole`](@ref)
"""
print_PS_multipole



##########################################################################################92


function all_PS_multipole(input::String,
     group::String=VALID_GROUPS[end];
     twofast::Bool = false,
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

     time_1 = time()

     ks, VEC = twofast ? TwoFAST_all_PS_multipole(input,
               group; L=L, pr=false, kwargs...) : 
          FFTLog_all_PS_multipole(input,
               group; L=L, pr=false, kwargs...)

     time_2 = time()

     pr && println("\ntime needed for all the Power Spectra computation [in s] = $(time_2-time_1)\n")

     return ks, VEC
end


function print_all_PS_multipole(input::String, out::String,
     group::String=VALID_GROUPS[end]; twofast::Bool = false,
     L::Int=0, pr::Bool=true, kwargs...)

     check_parent_directory(out)
     check_namefile(out)
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

     time_1 = time()

     ks, VEC = all_PS_multipole(input, group; 
          twofast = twofast, L = L, pr = false, kwargs...) 
     
     time_2 = time()

     N = length(ks)

     pr && println("\ntime needed for all the Power Spectra computation [in s] = $(time_2-time_1)\n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, BRAND)

          print(io, "#\n# Power Spectra Multipole computation ")
          if group == "GNC"
               println(io, "for the Galaxy Number Counts GR effect" *
                           "\n# from the file: $input")
          elseif group == "LD"
               println(io, "for the Luminosity Distance perturbations GR effect" *
                           "\n# from the file: $input")
          elseif group == "GNCxLD"
               println(io, "for the cross correlations between \n#" *
                           "Galaxy Number Counts and Luminosity Distance perturbations " *
                           "from the file:\n# $input")
          elseif group == "LDxGNC"
               println(io, "for the cross correlations between \n#" *
                           "Luminosity Distance perturbations and Galaxy Number Counts " *
                           "from the file:\n# $input")
          else
               println(io, "without any specific group considered" *
                           "\n# from the file: $input")
          end

          println(io, "#\n# For this PS_multipole computation we set: ")
          println(io, "# \t #points used in Fourier transform N = $N")
          println(io, "# \t multipole degree in consideration L = $L")
          println(io, "# overall computational time needed (in s) : $(@sprintf("%.4f", time_2-time_1))")
          print(io, "# kwards passed to \"print_all_PS_multipole\": ")

          if isempty(kwargs)
               println(io, "none")
          else
               print(io, "\n")
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end
          println(io, "# ")
          println(io, "# (all the following Power Spectra are measured in (Mpc/h_0)^3)")

          if group ≠ VALID_GROUPS[end]
               effs = group == "GNC" ? GR_EFFECTS_GNC :
                      group == "LD" ? GR_EFFECTS_LD :
                      group == "GNCxLD" ? GR_EFFECTS_GNCxLD :
                      group == "LDxGNC" ? GR_EFFECTS_LDxGNC :
                      throw(ErrorException("how did you arrive here???"))

               println(io, "# 1: k [h_0/Mpc] \t 2: P_SUM \t " *
                           join([string(i) for i in 3:length(effs)+2] .*
                                ": P_" .* effs .* " \t "))
          else
               println(io, "# 1: k [h_0/Mpc] \t" *
                           join([string(i) for i in 2:length(VEC)+1] .*
                                ": P_" .*
                                [string(i) for i in 2:length(VEC)+1] .* " \t ")
               )
          end

          for (i, k) in enumerate(ks)
               println(io, "$k \t" * join([GaPSE.number_to_string(v[i])*" \t " for v in VEC]))
          end
     end
end
