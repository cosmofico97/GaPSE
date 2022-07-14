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
     const DEFAULT_IPS_OPTS = Dict(
          :fit_left_min => 1e-6::Float64, 
          :fit_left_max => 3e-6::Float64,
          :fit_right_min => 1e1::Float64, 
          :fit_right_max => 2e1::Float64,
          )

The default values to be stored in `CosmoParams` concerning the 
Input Power Spectrum. In the `Cosmology` that will have such `CosmoParams` as input,
they will be used in its `InputPS`.

See also: [`CosmoParams`](@ref), [`Cosmology`](@ref), [`InputPS`](@ref)
"""
const DEFAULT_IPS_OPTS = Dict(
     :fit_left_min => 1e-6::Float64,
     :fit_left_max => 3e-6::Float64,
     :fit_right_min => 1e1::Float64,
     :fit_right_max => 2e1::Float64,
)


"""
     const DEFAULT_IPSTOOLS_OPTS = Dict(
          :N => 1024::Int,
          :fit_min => 0.05::Float64,
          :fit_max => 0.5::Float64,
          :con => true::Bool,
          :k_min => 1e-6::Float64,
          :k_max => 10.0::Float64,
     )

The default values to be stored in `CosmoParams` concerning the 
Input Power Spectrum Tools. In the `Cosmology` that will have such `CosmoParams` as input,
they will be used in its `IPSTools`.

See also: [`CosmoParams`](@ref), [`Cosmology`](@ref), [`IPSTools`](@ref)
"""
const DEFAULT_IPSTOOLS_OPTS = Dict(
     :fit_min => 0.05::Float64,
     :fit_max => 0.5::Float64,
     :N => 1024::Int,
     :con => true::Bool,
     :k_min => 1e-6::Float64,
     :k_max => 10.0::Float64,
)



"""
     const DEFAULT_WFI_OPTS = Dict(
          :llim=> 0.0::Float64,
          :rlim=> Inf::Float64,
          :N => 200::Int64,
          :trap => true::Bool,
          :rtol => 1e-2::Float64,
          :atol => 0.0::Float64,
          :ss_start => 0.0::Float64,
          :ss_step => 21.768735478453323::Float64,
          :ss_stop => 0.0::Float64,
          )

The default values to be stored in `CosmoParams` concerning the 
Integrated Iwndow Function F. In the `Cosmology` that will have such `CosmoParams` as input,
they will be used in its `WindowFIntegrated`.

See also: [`CosmoParams`](@ref), [`Cosmology`](@ref), [`WindowFIntegrated`](@ref),
"""
const DEFAULT_WFI_OPTS = Dict(
     :llim=> 0.0::Float64,
     :rlim=> Inf::Float64,
     :N => 200::Int64,
     :trap => true::Bool,
     :rtol => 1e-2::Float64,
     :atol => 0.0::Float64,
     :ss_start => 0.0::Float64,
     :ss_step => 21.768735478453323::Float64,
     :ss_stop => 0.0::Float64,
)



##########################################################################################92



"""
     CosmoParams(
          z_min::Float64
          z_max::Float64
          θ_max::Float64

          Ω_b::Float64
          Ω_cdm::Float64
          Ω_M0::Float64
          h_0::Float64

          b::Float64
          s_b::Float64
          𝑓_evo::Float64

          s_lim::Float64

          IPS::Dict{Symbol,T1} where {T1}
          IPSTools::Dict{Symbol,T2} where {T2}
          WFI::Dict{Symbol,T3} where {T3}
     )


Struct that contains all the parameters and options that are 
matter of concerns for the `Cosmology` we are interested in.

## Arguments

- `z_min::Float64` and `z_max::Float64` : the minimum and maximum redshifts of the
  survey we want to study.

- `θ_max::Float64` : Angular maximum value of the survey. It is
  implicitly assumed an azimutal simmetry of the survey.

- `Ω_b::Float64`, `Ω_cdm::Float64` and `Ω_M0::Float64` : barionic, cold-dark-matter and
  total matter density parameters.

- `h_0::Float64` : today's Hubble adimensional parameter (`H_0 = h_0 * 100 km/(s * Mpc)`).

- `b::Float64` : galaxy bias.

- `s_b::Float64` : magnification bias, i.e. the slope of the luminosity function at the luminosity threshold.

- `𝑓_evo::Float64` : evolution bias.

- `s_lim::Float64` : the lower-bound value for the function `func_ℛ_LD`; it is necessary, because
  `ℛ_LD` blows up for ``s \\rightarrow 0^{+}``. Consequently, if the `func_ℛ_LD` input value is 
  `0 ≤ s < s_lim`, the returned value is always `func_ℛ_LD(s_lim)`.

- `IPS::Dict{Symbol,T1} where {T1}` : dictionary concerning all the options that should be 
  passed to `InputPS` in the contruction of a `Cosmology`. The allowed keys, with their default
  values, are stored in  `DEFAULT_IPS_OPTS`, and are the following:
     - `:fit_left_min => 1e-6` and `:fit_left_max => 3e-6` : the limits (min and max) where the PS
       must be fitted with a (pure) power law, for small wavenumbers. 
     - `:fit_right_min => 1e1` and `:fit_right_max => 2e1` : the limits (min and max) where the PS
       must be fitted with a (pure) power law, for high wavenumbers. 

- `IPSTools::Dict{Symbol,T2} where {T2}` : dictionary concerning all the options that should be 
  passed to `IPSTools` in the contruction of a `Cosmology`. The allowed keys, with their default
  values, are stored in  `DEFAULT_IPSTOOLS_OPTS`, and are the following:
     - `:fit_min => 0.05` and `:fit_max => 0.5` : the limits (min and max) 
       where the integral ``I_\\ell^n`` in `Cosmology` must be fitted with a power law, 
       for small distances. This operation is necessary, because `xicalc`, in this context, 
       gives wrong results for too small input distance `s`; nevertheless, all these ``I_\\ell^n`` 
       integrals have fixed power-law trends for ``s \\rightarrow 0``, so this approach gives
       good results.
     - `:N => 1024` : number of points to be used in the Sperical Bessel Fourier Transform made
       by `xicalc` in `IPSTools`.
     - `:k_min => 1e-6` and `:k_max => 10.0` : extremes of integration for the `σ_i`
       integrals in `IPSTools`.
     - `:con => true` : do you want that the fit of all the ``I_\\ell^n`` in `IPSTools` for 
       the LEFT edge is not a simple power-law ``y = f(x) = b \\, x^s``, but also consider 
       a constant ``a``, such that ``y = f(x) = a + b \\, x^s``?


## Constructors

`CosmoParams(z_min, z_max, θ_max;
          Ω_b = 0.0489, Ω_cdm = 0.251020, h_0 = 0.70, s_lim = 1e-2,
          IPS_opts::Dict = Dict{Symbol,Any}(),
          IPSTools_opts::Dict = Dict{Symbol,Any}()
     )`
The associations are trivials, with `Ω_M0 = Ω_cdm + Ω_b`.
For the two dictionary, you may pass only the key and the value you are interested in,
and all the other default ones will be considered.
For example, if you set:

`IPSTools_opts = Dict(:N => 150, :con => false, :k_max => 30.0)`

then the dictionary with all the options that will be passed to `IPSTools` will be:

`IPSTools = merge(DEFAULT_IPSTOOLS_OPTS, IPSTools_opts) = 
     :fit_min => 0.05,   # default
     :fit_max => 0.5,    # default
     :N => 150,          # CHANGED VALUE
     :con => false,      # CHANGED VALUE
     :k_min => 1e-6,     # default
     :k_max => 30.0,     # CHANGED VALUE
)`

and similar for `IPS_opts`.


See also: [`Cosmology`](@ref), [`IPSTools`](@ref),  [`InputPS`](@ref), 
[`func_ℛ_LD`](@ref), [`DEFAULT_IPSTOOLS_OPTS`](@ref), [`DEFAULT_IPS_OPTS`](@ref),
[`DEFAULT_WFI_OPTS`](@ref), [`check_compatible_dicts`](@ref)
"""
struct CosmoParams
     z_min::Float64
     z_max::Float64
     θ_max::Float64

     Ω_b::Float64
     Ω_cdm::Float64
     Ω_M0::Float64
     h_0::Float64

     b::Float64
     s_b::Float64
     𝑓_evo::Float64

     s_lim::Float64

     IPS::Dict{Symbol,T1} where {T1}
     IPSTools::Dict{Symbol,T2} where {T2}
     WFI::Dict{Symbol,T3} where {T3}

     function CosmoParams(z_min, z_max, θ_max;
          Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.70, s_lim=1e-2,
          b=1.0, s_b=0.0, 𝑓_evo=0.0,
          IPS_opts::Dict=Dict{Symbol,Any}(),
          IPSTools_opts::Dict=Dict{Symbol,Any}(),
          WFI_opts::Dict=Dict{Symbol,Any}(),
     )
          str(n, a, b) = "the keys of the $n dict have to be Symbols (like :$a, :$b, ...)"
          
          @assert typeof(IPS_opts) <: Dict{Symbol,T1} where {T1} str("IPS_opts", "k_min", "N")

          @assert typeof(IPSTools_opts) <: Dict{Symbol,T2} where {T2} str("IPSTools_opts", "k_min", "N")

          @assert typeof(WFI_opts) <: Dict{Symbol,T3} where {T3} str("WFI_opts", "r_lim", "N")

          check_compatible_dicts(DEFAULT_IPS_OPTS, IPS_opts, "IPS_opts")
          check_compatible_dicts(DEFAULT_IPSTOOLS_OPTS, IPSTools_opts, "IPSTools_opts")
          check_compatible_dicts(DEFAULT_WFI_OPTS, WFI_opts, "WFI_opts")

          IPS = merge(DEFAULT_IPS_OPTS, IPS_opts)
          IPSTools = merge(DEFAULT_IPSTOOLS_OPTS, IPSTools_opts)
          WFI = merge(DEFAULT_WFI_OPTS, WFI_opts)

          @assert 0.0 < z_min < z_max " 0.0 < z_min < z_max must hold!"
          @assert 0.0 ≤ θ_max ≤ π / 2.0 " 0.0 ≤ θ_max ≤ π/2.0 must hold!"
          @assert 0.0 ≤ Ω_b ≤ 1.0 " 0.0 ≤ Ω_b ≤ 1.0 must hold!"
          @assert 0.0 ≤ Ω_cdm ≤ 1.0 " 0.0 ≤ Ω_cdm ≤ 1.0 must hold!"
          @assert 0.0 < h_0 ≤ 1.0 " 0.0 < h_0 ≤ 1.0 must hold!"
          @assert 0.0 < s_lim < 10.0 "0.0 < s_lim < 10.0 must hold!"

          @assert 0.0 < IPS[:fit_left_min] < IPS[:fit_left_max] < 1e-1
          " 0 < fit_left_min < fit_left_max < 0.1 must hold!"
          @assert 0.5 < IPS[:fit_right_min] < IPS[:fit_right_max] < 1e6
          " 0.5 < fit_right_min < fit_right_max < 1e6 must hold!"

          @assert 0.0 ≤ IPSTools[:k_min] < IPSTools[:k_max] " 0.0 ≤ k_min < k_max must hold!"
          @assert IPSTools[:N] > 7 " N > 7 must hold!"
          @assert 1e-2 ≤ IPSTools[:fit_min] < IPSTools[:fit_max] < 10.0 " 1e-2 " *
                                                                        "≤ fit_min < fit_max < 10.0 must hold!"
                                                            
          @assert b > 0.0 " b > 0 must hold!"

          @assert 0.0 ≤ WFI[:llim] < WFI[:rlim] " 0.0 ≤ llim < rlim must hold!"
          @assert WFI[:N] > 10 " N > 10 must hold!"
          @assert 0 < WFI[:rtol] < 1 " 0 < rtol < 1 must hold!" 
          @assert 0 ≤ WFI[:atol] < 1 " 0 ≤ atol < 1 must hold!" 
          @assert WFI[:ss_start] ≥ 0.0 " ss_start ≥ 0.0 must hold!"
          @assert WFI[:ss_step] > 0 " ss_step > 0 must hold!"
          @assert iszero(WFI[:ss_stop]) || WFI[:ss_stop] > WFI[:ss_start] "ss_stop == 0 or ss_stop > ss_start must hold!"

          new(z_min, z_max, θ_max, 
               Ω_b, Ω_cdm, Ω_cdm + Ω_b, h_0, 
               b, s_b, 𝑓_evo, 
               s_lim,
               IPS, IPSTools, WFI)
     end
end

