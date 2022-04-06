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

function PS_multipole(f_in;
     int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
     L::Integer=0, N::Integer=1024, pr::Bool=true,
     k0::Union{Nothing,Float64}=nothing,
     right::Union{Float64,Nothing}=nothing)

     k_0 = isnothing(k0) ? 1.0 / int_s_max : k0

     t1 = time()
     ks, pks = xicalc(s -> 2 * π^2 * f_in(s), L, 0;
          N=N, kmin=int_s_min, kmax=int_s_max, r0=k_0)
     t2 = time()
     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(t2-t1)\n")

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


function PS_multipole(ss, fs;
     int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
     epl::Bool=true,
     N_left::Integer=12, N_right::Integer=12,
     p0_left=[-2.0, 1.0], p0_right=[-2.0, 1.0],
     kwargs...)

     @assert length(ss) == length(fs) "xs and ys must have same length"
     k_0 = 1.0 / max(ss...)
     right = 1.0 / min(ss...)

     f_in, INT_s_min, INT_s_max =
          if epl == true
               if all(fs[begin:begin+5] .≈ 0.0) && all(fs[end-5:end] .≈ 0.0)
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

     return PS_multipole(f_in; int_s_min=INT_s_min, int_s_max=INT_s_max,
          k0=k_0, right=right, kwargs...)
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
          N_left::Integer=12, N_right::Integer=12,
          p0_left=[-2.0, 1.0], p0_right=[-2.0, 1.0],
          kwargs...)

     PS_multipole(f_in;
          int_s_min::Float64=1e-1, int_s_max::Float64=1e3,
          L::Integer=0, N::Integer=1024, pr::Bool=true,
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

- `N_left::Integer = 12` and `N_right::Integer = 12` : number of points from left
  right edges to be used for the power law fitting in `EPLs`. They matters only
  if in the given input file ξ is not defined until the extremes of integration
  `int_s_min` and `int_s_max`

- `L::Integer = 0`: order of the Legendre polynomial to be used; note that 
  the multipole order `L` must be the same of the input TPCF in exam! 
  Otherwise, the results would have no sense at all!

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `N::Integer = 1024` : number of points to be returned by `xicalc`


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
     L::Integer=0, N::Integer=1024,
     pr::Bool=true, kwargs...)

     pr && println("\nI'm computiong the PS_multipole from the file $input")

     time_1 = time()
     vec = PS_multipole(input; N=N, L=L, pr=pr, kwargs...)
     time_2 = time()

     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(time_2-time_1)\n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
          println(io, "# Power Spectrum Multipole computation of the file: $input")
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
     L::Integer=0, N::Integer=1024,
     pr::Bool=true, kwargs...)

     pr && println("\nI'm computiong the PS_multipole from the two input vectors.")

     time_1 = time()
     vec = PS_multipole(ss, fs; N=N, L=L, pr=pr, kwargs...)
     time_2 = time()

     pr && println("\ntime needed for Power Spectrum  computation [in s] = $(time_2-time_1)\n")

     isfile(out) && run(`rm $out`)
     open(out, "w") do io
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
               println(io, "$k \t $pk")
          end
     end
end


"""
     print_PS_multipole(input::String, out::String;
          L::Integer = 0, N::Integer = 1024,
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

- `L::Integer = 0`: order of the Legendre polynomial to be used; note that 
  the multipole order `L` must be the same of the input TPCF in exam! 
  Otherwise, the results would have no sense at all!

- `pr::Bool = true` : do you want the progress bar showed on screen, in order to 
  check the time needed for the computation? (`true` recommended)

- `N::Integer = 1024` : number of points to be returned by `xicalc`

- `kwargs...` : other keyword arguments that will be passed to `PS_multipole`

See also: [`V_survey`](@ref), [`A`](@ref), [`A_prime`](@ref),
[`EPLs`](@ref), [`PS_multipole`](@ref)
"""
print_PS_multipole
