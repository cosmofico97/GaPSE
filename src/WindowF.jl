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
     DEFAULT_FMAP_OPTS_hcub = Dict(
          :θ_max => π / 2.0::Float64, 
          :tolerance => 1e-8::Float64, 
          :rtol => 1e-2::Float64, 
          :atol => 1e-3::Float64,
          :pr => true::Bool,
     )

The default values to be used for the `F` function when you
want to perform the computation with `hcubature`.

See also: [`integrand_F`](@ref), [`F_hcub`](@ref), [`map_F`](@ref)
"""
const DEFAULT_FMAP_OPTS_hcub = Dict(
     :θ_max => π / 2.0::Float64, 
     :tolerance => 1e-8::Float64, 
     :rtol => 1e-2::Float64, 
     :atol => 1e-3::Float64,
     :pr => true::Bool,
)


"""
     DEFAULT_FMAP_OPTS_trap = Dict(
          :θ_max => π / 2.0::Float64, 
          :tolerance => 1e-8::Float64, 
          :N => 300::Int64, 
          :en => 1.0::Float64,
          :pr => true::Bool,
     )


The default values to be used for the `F` function when you
want to perform the computation with `trap`.

See also: [`integrand_F`](@ref), [`F_trap`](@ref), [`map_F`](@ref)
"""
const DEFAULT_FMAP_OPTS_trap = Dict(
     :θ_max => π / 2.0::Float64, 
     :tolerance => 1e-8::Float64, 
     :N => 300::Int64, 
     :en => 1.0::Float64,
     :pr => true::Bool,
)


##########################################################################################92



"""
    integrand_F(θ_1, θ, x, μ, θ_max; tolerance=1e-8) :: Float64

Return the integrand of the function ``F(x,\\mu; \\theta_\\mathrm{max})``, i.e the 
function ``f(x,\\mu, \\theta, \\theta_1; \\theta_\\mathrm{max})``:

```math
\\begin{split}
f(x,\\mu, \\theta, \\theta_1; \\theta_\\mathrm{max}) = 
    \\; &\\Theta\\left( \\frac
        {x \\cos \\theta + \\cos \\theta_1}{\\sqrt{x^1+2+2x\\mu}} - 
        \\cos(\\theta_\\mathrm{max}) 
        \\right) 
    \\; \\times \\; \\Theta(\\mu-\\cos(\\theta+\\theta_1)) \\; \\, \\times \\\\
    & \\quad \\Theta(\\cos(\\theta - \\theta_1)-\\mu) \\; \\times \\;
    \\frac{4\\pi \\sin\\theta\\sin\\theta_1}
        {\\sqrt{(\\sin\\theta\\sin\\theta_1)^2-(\\cos\\theta\\cos\\theta_1-\\mu)^2}}
\\end{split}
```
```math
\\begin{equation}
F(x,\\mu; \\theta_\\mathrm{max}) = \\int_0^{\\theta_\\mathrm{max}} 
        \\mathrm{d}\\theta_1 \\int_0^\\pi \\mathrm{d} \\theta 
        \\; f(x,\\mu, \\theta, \\theta_1; \\theta_\\mathrm{max})
\\end{equation}
```

`tolerance` is a parameter needed in case of small negative denominator: the Heaviside
theta function mathematically prevent that 
``\\mathrm{den}=(\\sin\\theta\\sin\\theta_1)^2-(\\cos\\theta\\cos\\theta_1-\\mu)^2``
becomes negative, but computationally might happen that ``\\mathrm{den}`` results as a
very small negative number (for instance `-1.2368946523-18`); in this case `tolerance`
solve the problem, returning 0 if ``0<-\\mathrm{den}< \\mathrm{tolerance}``

See also: [`F`](@ref), [`map_F`](@ref)
"""
function integrand_F(θ_1, θ, x, μ, θ_max; tolerance = 1e-8)
     if (x * cos(θ) + cos(θ_1)) / √(x^2 + 1 + 2 * x * μ) - cos(θ_max) > 0 &&
        (μ - cos(θ + θ_1)) > 0 &&
        (cos(θ - θ_1) - μ) > 0
          den = (sin(θ) * sin(θ_1))^2 - (cos(θ) * cos(θ_1) - μ)^2
          if den > 0
               return 4.0 * π * sin(θ_1) * sin(θ) / √den
          elseif -den < tolerance
               return 0.0
          else
               throw(ErrorException("negative denominator, greater than $(tolerance): $(den)"))
          end
     else
          return 0.0
     end
end;



"""
     F_hcub(x, μ; θ_max = π/2, tolerance = 1e-8, 
          atol = 1e-2, rtol = 1e-5, 
          kwargs...) ::Tuple{Float64, Float64}

Computes with `hcubature` the value of ``F(x,\\mu; \\theta_\\mathrm{max})``, 
defined as follows:

```math
\\begin{split}
F(x,\\mu; \\theta_\\mathrm{max}) = & \\;4\\pi 
    \\int_0^{\\theta_\\mathrm{max}} \\mathrm{d}\\theta_1 \\int_0^\\pi \\mathrm{d} \\theta \\; 
    \\, \\Theta\\left(\\frac
        {x \\cos \\theta + \\cos \\theta_1}{\\sqrt{x^1+2+2x\\mu}} - 
        \\cos(\\theta_\\mathrm{max}) 
        \\right) 
    \\, \\Theta(\\mu-\\cos(\\theta+\\theta_1)) \\\\
    &\\Theta(\\cos(\\theta - \\theta_1)-\\mu) \\;
    \\frac{\\sin\\theta\\sin\\theta_1}
        {\\sqrt{(\\sin\\theta\\sin\\theta_1)^2-(\\cos\\theta\\cos\\theta_1-\\mu)^2}}
\\end{split}
```

`tolerance` is a parameter needed in case of small negative denominator: the Heaviside
theta function mathematically prevent that 
``\\mathrm{den}=(\\sin\\theta\\sin\\theta_1)^2-(\\cos\\theta\\cos\\theta_1-\\mu)^2``
becomes negative, but computationally might happen that ``\\mathrm{den}`` results as a
very small negative number (for instance `-1.2368946523-18`); in this case `tolerance`
solve the problem, returning 0 if ``0<-\\mathrm{den}< \\mathrm{tolerance}``.

The double integral is performed with [`hcubature`](@ref) function from the Julia
Package [`HCubature`](@ref); `rtol`, `atol` and all the `kwargs` insert into `F` 
are directly transferred to `hcubature`. 

The output of this function is a `Tuple{Float64, Float64}`, containing respectively
the value of the integral ad its error.

PAY ATTENTION: do not set too small `atol` and `rtol`, or the computation
can easily become overwhelming! 

NOTE: for computational efficiency and stability, it is highly recommended to use 
the other function`F_trap`, based on the trapezoidal rule, in order to compute
this function.

See also:  [`F_trap`](@ref), [`map_F`](@ref), [`integrand_F`](@ref), 
[`check_compatible_dicts`](@ref)
"""
function F_hcub(x, μ; θ_max = π/2, tolerance = 1e-8, atol = 1e-2, rtol = 1e-5, kwargs...)
     my_int(var) = integrand_F(var[1], var[2], x, μ, θ_max; tolerance = tolerance)
     a = [0.0, 0.0]
     b = [θ_max, π]
     return hcubature(my_int, a, b; rtol = rtol, atol = atol, kwargs...)
end;



"""
     F_trap(x, μ; θ_max = π/2, tolerance = 1e-8, 
          atol = 1e-2, rtol = 1e-5, 
          kwargs...) ::Float64

Computes with `trap` the value of ``F(x,\\mu; \\theta_\\mathrm{max})``, 
defined as follows:

```math
\\begin{split}
F(x,\\mu; \\theta_\\mathrm{max}) = & \\;4\\pi 
    \\int_0^{\\theta_\\mathrm{max}} \\mathrm{d}\\theta_1 \\int_0^\\pi \\mathrm{d} \\theta \\; 
    \\, \\Theta\\left(\\frac
        {x \\cos \\theta + \\cos \\theta_1}{\\sqrt{x^1+2+2x\\mu}} - 
        \\cos(\\theta_\\mathrm{max}) 
        \\right) 
    \\, \\Theta(\\mu-\\cos(\\theta+\\theta_1)) \\\\
    &\\Theta(\\cos(\\theta - \\theta_1)-\\mu) \\;
    \\frac{\\sin\\theta\\sin\\theta_1}
        {\\sqrt{(\\sin\\theta\\sin\\theta_1)^2-(\\cos\\theta\\cos\\theta_1-\\mu)^2}}
\\end{split}
```

`tolerance` is a parameter needed in case of small negative denominator: the Heaviside
theta function mathematically prevent that 
``\\mathrm{den}=(\\sin\\theta\\sin\\theta_1)^2-(\\cos\\theta\\cos\\theta_1-\\mu)^2``
becomes negative, but computationally might happen that ``\\mathrm{den}`` results as a
very small negative number (for instance `-1.2368946523-18`); in this case `tolerance`
solve the problem, returning 0 if ``0<-\\mathrm{den}< \\mathrm{tolerance}``.

The double integral is performed with [`trapz`](@ref) function from the Julia
Package [`Trapz`](@ref), that is based on the trapezoidal rule. `N` is the number
of point to be used to sample INDIPENDENTLY `θ_1` and `θ`, so consider that there
is a `N^2` time dependence.
It's recommended to set `100 < N < 1000`.

NOTE: there is another function, called `F_hcub`, that performs this calculus. 
Nevertheless, for computational efficiency and stability, it is highly recommended to use 
to use this one.

See also:  [`F_hcub`](@ref), [`map_F`](@ref), [`integrand_F`](@ref), 
[`check_compatible_dicts`](@ref)
"""
function F_trap(x, μ; θ_max = π/2, N::Int=300, en=1.0, 
        tolerance = 1e-13)
    
     χ1s = θ_max .* range(0.0, 1.0, length=N)
     χ2s = π .* range(0.0, 1.0, length=N+7)

     zs = [
          en * integrand_F(χ1, χ2, x, μ, θ_max; tolerance = tolerance)
          for χ1 in χ1s, χ2 in χ2s
     ]

     res = trapz((χ1s, χ2s), zs)
    
     return res / en
end;


function print_map_F(out::String, x_step::Float64 = 0.01, μ_step::Float64 = 0.01;
     trap::Bool = true, x1 = 0, x2 = 3, μ1 = -1, μ2 = 1,
     Fmap_opts::Dict = Dict{Symbol,Any}(), kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     @assert x1 >= 0.0 "The lower limit of x must be >0, not $(x1)!"
     @assert x2 > x1 "The upper limit of x must be > than the lower one, not $(x2)<=$(x1)!"
     @assert x2 <= 10.0 "The upper limit of x must be <10, not $(x2)!"
     @assert μ1 >= -1.0 "The lower limit of μ must be >=-1, not $(μ1)!"
     @assert μ2 <= 1.0 "The upper limit of μ must be <=1, not $(μ2)!"
     @assert μ2 > μ1 "The upper limit of μ must be > than the lower one, not $(μ2)<=$(μ1)!"
     @assert 0 < x_step < 1 "The integration step of x must be 0<x_step<1, not $(x_step)!"
     @assert 0 < μ_step < 1 "The integration step of μ must be 0<μ_step<1, not $(μ_step)!"

     @assert typeof(Fmap_opts) <: Dict{Symbol,T1} where {T1}
          "the keys of the Fmap_opts dict have to be Symbols (like :k_min, :N, ...)"
     trap ? check_compatible_dicts(DEFAULT_FMAP_OPTS_trap, Fmap_opts, "Fmap_opts") :
               check_compatible_dicts(DEFAULT_FMAP_OPTS_hcub, Fmap_opts, "Fmap_opts")
     Fmap_dict = trap ? merge(DEFAULT_FMAP_OPTS_trap, Fmap_opts) : merge(DEFAULT_FMAP_OPTS_hcub, Fmap_opts)

     μs = μ1:μ_step:μ2
     xs = x1:x_step:x2
     xs_grid = [x for x = xs for μ = μs]
     μs_grid = [μ for x = xs for μ = μs]

     time_1 = time()

     new_F(x, μ) = trap ?
          F_trap(x, μ; θ_max = Fmap_dict[:θ_max], 
               tolerance = Fmap_dict[:tolerance], N = Fmap_dict[:N], 
               en = Fmap_dict[:en]) :
          F_hcub(x, μ; θ_max = Fmap_dict[:θ_max], 
          tolerance = Fmap_dict[:tolerance], atol = Fmap_dict[:atol], 
          rtol = Fmap_dict[:rtol], kwargs...)    

     Fs_grid = Fmap_dict[:pr] ? begin
          @showprogress "window F evaluation: " map(new_F, xs_grid, μs_grid)
     end : map(new_F, xs_grid, μs_grid)

     time_2 = time()

     #run(`rm $(out)`)
     open(out, "w") do io
          println(io, BRAND)
          println(io, "# This is an integration map on μ and x of the window function F(x, μ)")
          println(io, "# For its analytical definition, check the code.\n#")
          println(io, "# It was set trap = $trap, so the computation was performed with "*
                         begin trap ? "trap()" : "hcubature()" end )
          println(io, "# Parameters used in this integration map:")
          println(io, "# x_min = $(x1) \t x_max = $(x2) \t x_step = $(x_step)")
          println(io, "# mu_min = $(μ1) \t mu_max = $(μ2) \t mu_step = $(μ_step)")
          println(io, "# computational time (in s) : $(@sprintf("%.3f", time_2-time_1))")
          println(io, "# kwards passed: ")

          for key in keys(Fmap_dict)
               println(io, "# \t\t$(key) = $(Fmap_dict[key])")
          end

          if !isempty(kwargs)
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          if trap == true
               println(io, "#\n# x \t mu \t F")
               for (x, μ, F) in zip(xs_grid, μs_grid, Fs_grid)
                    println(io, "$x\t $μ \t $F")
               end
          else
               println(io, "#\n# x \t mu \t F \t F_error")
               for (x, μ, F) in zip(xs_grid, μs_grid, Fs_grid)
                    println(io, "$x\t $μ \t $(F[1]) \t $(F[2])")
               end
          end
     end
end



function print_map_F(out::String, xs::Vector{Float64}, μs::Vector{Float64}; 
     trap::Bool = true, Fmap_opts::Dict = Dict{Symbol,Any}(), kwargs...)

     check_parent_directory(out)
     check_namefile(out)

     @assert all(xs .>= 0.0) "All xs must be >=0.0!"
     @assert all([xs[i+1] > xs[i] for i in 1:(length(xs)-1)]) "xs must be a float vector of increasing values!"
     @assert all(xs .<= 10.0) "All xs must be <=10.0!"
     @assert all(μs .>= -1.0) "All μs must be >=-1.0!"
     @assert all([μs[i+1] > μs[i] for i in 1:(length(μs)-1)]) "μs must be a float vector of increasing values!"
     @assert all(μs .<= 1.0) "All μs must be <=1.0!"

     @assert typeof(Fmap_opts) <: Dict{Symbol,T1} where {T1}
          "the keys of the Fmap_opts dict have to be Symbols (like :k_min, :N, ...)"
     trap ? check_compatible_dicts(DEFAULT_FMAP_OPTS_trap, Fmap_opts, "Fmap_opts") :
               check_compatible_dicts(DEFAULT_FMAP_OPTS_hcub, Fmap_opts, "Fmap_opts")
     Fmap_dict = trap ? merge(DEFAULT_FMAP_OPTS_trap, Fmap_opts) : merge(DEFAULT_FMAP_OPTS_hcub, Fmap_opts)


     xs_grid = [x for x = xs for μ = μs]
     μs_grid = [μ for x = xs for μ = μs]

     time_1 = time()

     new_F(x, μ) = trap ?
          F_trap(x, μ; θ_max = Fmap_dict[:θ_max], 
               tolerance = Fmap_dict[:tolerance], N = Fmap_dict[:N], 
               en = Fmap_dict[:en]) :
          F_hcub(x, μ; θ_max = Fmap_dict[:θ_max], 
          tolerance = Fmap_dict[:tolerance], atol = Fmap_dict[:atol], 
          rtol = Fmap_dict[:rtol], kwargs...) 

     Fs_grid = Fmap_dict[:pr] ? begin
          @showprogress "window F evaluation: " map(new_F, xs_grid, μs_grid)
     end : map(new_F, xs_grid, μs_grid)

     time_2 = time()

     #run(`rm $(out)`)
     open(out, "w") do io
          println(io, BRAND)
          println(io, "# This is an integration map on μ and x of the window function F(x, μ)")
          println(io, "# It was set trap = $trap, so the computation was performed with "*
                         begin trap ? "trap()" : "hcubature()" end )
          println(io, "# For its analytical definition, check the code.\n#")
          println(io, "# Parameters used in this integration map:")
          println(io, "# computational time (in s) : $(@sprintf("%.3f", time_2-time_1))")
          println(io, "# kwards passed: ")

          for key in keys(Fmap_dict)
               println(io, "# \t\t$(key) = $(Fmap_dict[key])")
          end

          if !isempty(kwargs)
               for key in keys(kwargs)
                    println(io, "# \t\t$(key) = $(kwargs[key])")
               end
          end

          if trap == true
               println(io, "#\n# x \t mu \t F")
               for (x, μ, F) in zip(xs_grid, μs_grid, Fs_grid)
                    println(io, "$x\t $μ \t $F")
               end
          else
               println(io, "#\n# x \t mu \t F \t F_error")
               for (x, μ, F) in zip(xs_grid, μs_grid, Fs_grid)
                    println(io, "$x\t $μ \t $(F[1]) \t $(F[2])")
               end
          end
     end
end




"""
     print_map_F(out::String, x_step::Float64 = 0.01, μ_step::Float64 = 0.01;
          trap::Bool = true, x1 = 0, x2 = 3, μ1 = -1, μ2 = 1, 
          Fmap_opts::Dict = Dict{Symbol,Any}(), 
          kwargs...)

     print_map_F(out::String, xs::Vector{Float64}, μs::Vector{Float64};
          trap::Bool = true, Fmap_opts::Dict = Dict{Symbol,Any}(),
          kwargs...)

Evaluate the window function ``F(x,\\mu; \\theta_\\mathrm{max})`` in a rectangual grid 
of ``\\mu`` and ``x`` values, and print the results in the `out` file.

In the first method you have to specify manually, both for `x` and `μ`, start
(`x1` and `μ1`), stop (`x2` and `μ2`), and step (`x_step` and `μ_step`).
In the second one, you need to pass the values you want to calculate 
the function in, through the vectors `xs` and `μs`.

The bool variable `trap` tells if you want to perform the computation of `F` with
`F_trap` (if `trap==true`) or with `F_hcub` (if `trap==false`).
Both for computational efficiency and stability, it's highly recommended to use
the former (i.e. the default one). 

`Fmap_opts` is instead the way you should exploit in order to pass to `F_trap`/`F_hcub`
other options you are interested in.
You may pass only the key and the value you are focused on,
and all the other default ones will be considered.

For example, if you set `trap == false` and:

`Fmap_opts = Dict(:tolerance => 1e-5, :θ_max => 2.0)`

then the dictionary with all the options that will be passed to `F` will be:

`Fmap_dict = merge(DEFAULT_FMAP_OPTS_hcub, Fmap_opts) = 
     :θ_max => 2.0,           # CHANGED VALUE
     :tolerance => 1e-5,      # CHANGED VALUE
     :rtol => 1e-2,           # default
     :atol => 1e-3,           # default
     :pr => true,             # default
)`

Check the documentation of `DEFAULT_FMAP_OPTS_hcub` and `DEFAULT_FMAP_OPTS_trap`
for more information about these default values.


See also: [`DEFAULT_FMAP_OPTS_hcub`](@ref), [`DEFAULT_FMAP_OPTS_trap`](@ref)
[`F_trap`](@ref), [`F_hcub`](@ref), [`integrand_F`](@ref)
"""
print_map_F



##########################################################################################92



"""
    WindowF(
        xs::Vector{Float64}
        μs::Vector{Float64}
        Fs::Matrix{Float64}
        )

Struct containing xs, μs and Fs values of the window function ``F(x, μ)``.
`xs` and `μs` are 1D vectors containing each value only once, while 
`Fs` values are contained in a matrix of size `(length(xs), length(μs))`, so:
- along a fixed column the changing value is `x`
- along a fixed row the changing value is `μ`

## Constructors

`WindowF(file::String)` : read the F map from the file `file`. Such a file might
be produced by `print_map_F`, check its docstring. 

It does not matter if the pattern is

```data
# xs      μs      Fs
0.0       -1.0       ...
0.0       -0.9       ...
0.0       -0.8       ...
...       ...      ...
```

or 

```data
# xs      μs      Fs
0.0       -1.0       ...
0.1       -1.0       ...
0.2       -1.0       ...
...       ...      ...
```

because the constructor will recognise it. What does matter is the columns order:
`xs` first, then `μs` and finally `Fs`.

See also: [`print_map_F`](@ref), [`F_trap`](@ref), [`spline_F`](@ref)
"""
struct WindowF
     xs::Vector{Float64}
     μs::Vector{Float64}
     Fs::Matrix{Float64}


     function WindowF(file::String)
          data = readdlm(file, comments = true)
          xs, μs, Fs = data[:, 1], data[:, 2], data[:, 3]
          @assert size(xs) == size(μs) == size(Fs) "xs, μs and Fs must have the same length!"

          new_xs = unique(xs)
          new_μs = unique(μs)
          new_Fs =
               if xs[2] == xs[1] && μs[2] ≠ μs[1]
                    transpose(reshape(Fs, (length(new_μs), length(new_xs))))
               elseif xs[2] ≠ xs[1] && μs[2] == μs[1]
                    reshape(Fs, (length(new_xs), length(new_μs)))
               else
                    throw(ErrorException("What kind of convenction for the file $file" *
                                         " are you using? I do not recognise it."))
               end
          new(new_xs, new_μs, new_Fs)
     end
end


"""
     spline_F(x, μ, str::WindowF)) ::Float64

Return the 2-dim spline value of ``F`` in the given `(x,μ)`, where
``F`` is defined in the input `WindowF`.
The spline is obtained through the `interpolate` function of the 
[`GridInterpolations`](https://github.com/sisl/GridInterpolations.jl) Julia
package.

See also: [`WindowF`](@ref)
"""
function spline_F(x, μ, str::WindowF)
     grid = GridInterpolations.RectangleGrid(str.xs, str.μs)
     GridInterpolations.interpolate(grid, reshape(str.Fs, (:, 1)), [x, μ])
end

#=
map_F_data = readdlm(FILE_F_MAP, comments = true)
map_F_data_dict = Dict([name => map_F_data[2:end, i] for (i, name) in enumerate(NAMES_F_MAP)]...)

_xs = unique(map_F_data_dict["x"])
_μs = unique(map_F_data_dict["mu"])
_Fs = map_F_data_dict["F"]


# for my F map convenction
my_F_grid = GridInterpolations.RectangleGrid(_μs, _xs)
spline_F(x, μ) = GridInterpolations.interpolate(my_F_grid, _Fs, [μ, x])

# for the opposite convenction for F map
#other_F_grid = GridInterpolations.RectangleGrid(_xs, _μs)
#spline_F(x, μ) = GridInterpolations.interpolate(other_F_grid, _Fs, [x, μ])
=#

