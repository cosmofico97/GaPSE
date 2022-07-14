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

#=
function my_Plm(l::Int, m::Int, x)
     @assert l ≥ 0 " l ≥ 0 must hold, while $l does not!"
     @assert l ≥ abs(m) " l ≥ |m| must hold, while $l and $m do not!"
     if m ≥ 0
          return Plm(l, m, x)
     else
          return (-1)^m * factorial(l - m) / factorial(l + m) * Plm(l, -m, x)
     end
end;



"""
     Aabcd(θ_max, a::Int, b::Int, c::Int, d::Int; kwargs...)

Calculates the following expression:

```math
\\begin{align*}
     A_{a b}^{c d}  
     &= \\int \\mathrm{d}\\theta \\sin \\left( \\theta \\right) W \\left( \\theta \\right) 
          P_a^c \\left( \\cos \\theta \\right)  P_b^d \\left( \\cos \\theta \\right) \\
     &= \\int_0^{\\theta_{\\mathrm{max}} \\mathrm{d}\\theta \\sin \\left( \\theta \\right) 
          P_a^c \\left( \\cos \\theta \\right)  P_b^d \\left( \\cos \\theta \\right) \\
     &= \\int_{\\cos \\theta_{\\mathrm{max}}}^1 \\mathrm{d}x 
          P_a^c \\left( x \\right)  P_b^d \\left( x \\right)
\\end{align*}
```

The integral is performed thanks to `quadgk`; all the `kwargs...` are passed to that function.
"""
function Aabcd(θ_max, a::Int, b::Int, c::Int, d::Int; kwargs...)
     quadgk(x -> my_Plm(a, c, x) * my_Plm(b, d, x), cos(θ_max), 1; kwargs...)[1]
end


function vec_index_Aabcd(l::Int)
     @assert l ≥ 0 " l ≥ 0 must hold!"
     matr = [[0, 0, 0, 0] for i in 1:(l+1)^4]
     i = 1
     for a in 0:l, b in 0:l
          for c in range(-a, a), d in range(-b, b)
               matr[i] = [a, b, c, d]
               i += 1
          end
     end
     matr
end



##########################################################################################92



function Wnnll(θ_max::Float64, n1::Int, n2::Int, l3::Int, l4::Int;
     rtol=1e-2, kwargs...)

     vec_indexes_A = vec_index_Aabcd(max(n1, n2, l3, l4))
     vec_Aabcd = [Aabcd(θ_max, indexes_A...; rtol=rtol, kwargs...) for indexes_A in vec_indexes_A]
     dict_index_Aabcd = Dict(x => y for (x, y) in zip(vec_indexes_A, vec_Aabcd))

     16 * π^4 * sum(
          [
          begin
               #println("n1 = $(n1) \t l3 = $(l3) \t m3 = $(m3)");
               #println("n1-m3 = $(n1-m3) \t l3-m3 = $(l3-m3) \t n2-m4 = $(n2-m4) \t l4-m4 = $(l4-m4)");
               #println("n1+m3 = $(n1+m3) \t l3+m3 = $(l3+m3) \t n2+m4 = $(n2+m4) \t l4+m4 = $(l4+m4)");
               factorial(n1 - m3) / factorial(n1 + m3) * factorial(l3 - m3) / factorial(l3 + m3) *
               factorial(n2 - m4) / factorial(n2 + m4) * factorial(l4 - m4) / factorial(l4 + m4) *
               dict_index_Aabcd[[n1, l3, m3, m3]] *
               dict_index_Aabcd[[n2, l3, m4, m3]] *
               dict_index_Aabcd[[l4, n1, m4, m3]] *
               dict_index_Aabcd[[n2, l4, m4, m4]]
          end
          for m3 in -l3:l3, m4 in -l4:l4
     ]
     )
end

function Wnnll(dict_i_A::Dict{Vector{Int64},Float64},
     n1::Int, n2::Int, l3::Int, l4::Int)
     16 * π^4 * sum(
          [
          begin
               #println("n1 = $(n1) \t l3 = $(l3) \t m3 = $(m3)");
               #println("n1-m3 = $(n1-m3) \t l3-m3 = $(l3-m3) \t n2-m4 = $(n2-m4) \t l4-m4 = $(l4-m4)");
               #println("n1+m3 = $(n1+m3) \t l3+m3 = $(l3+m3) \t n2+m4 = $(n2+m4) \t l4+m4 = $(l4+m4)");
               factorial(n1 - m3) / factorial(n1 + m3) * factorial(l3 - m3) / factorial(l3 + m3) *
               factorial(n2 - m4) / factorial(n2 + m4) * factorial(l4 - m4) / factorial(l4 + m4) *
               dict_i_A[[n1, l3, m3, m3]] *
               dict_i_A[[n2, l3, m4, m3]] *
               dict_i_A[[l4, n1, m4, m3]] *
               dict_i_A[[n2, l4, m4, m4]]
          end
          for m3 in -l3:l3, m4 in -l4:l4
     ]
     )
end

function Wnnll(θ_max::Float64, indexes::Vector{Vector{Int64}};
     rtol=1e-2, kwargs...)

     l_max = max([max(x...) for x in indexes]...)
     vec_indexes_A = vec_index_Aabcd(l_max)
     vec_Aabcd = [Aabcd(θ_max, indexes_A...; rtol=rtol, kwargs...) for indexes_A in vec_indexes_A]
     dict_index_Aabcd = Dict(x => y for (x, y) in zip(vec_indexes_A, vec_Aabcd))

     [
          begin
               n1, n2, l3, l4 = x[1], x[2], x[3], x[4]
               16 * π^4 * sum(
                    [
                    begin
                         #println("n1 = $(n1) \t l3 = $(l3) \t m3 = $(m3)");
                         #println("n1-m3 = $(n1-m3) \t l3-m3 = $(l3-m3) \t n2-m4 = $(n2-m4) \t l4-m4 = $(l4-m4)");
                         #println("n1+m3 = $(n1+m3) \t l3+m3 = $(l3+m3) \t n2+m4 = $(n2+m4) \t l4+m4 = $(l4+m4)");
                         factorial(n1 - m3) / factorial(n1 + m3) * factorial(l3 - m3) / factorial(l3 + m3) *
                         factorial(n2 - m4) / factorial(n2 + m4) * factorial(l4 - m4) / factorial(l4 + m4) *
                         dict_index_Aabcd[[n1, l3, m3, m3]] *
                         dict_index_Aabcd[[n2, l3, m4, m3]] *
                         dict_index_Aabcd[[l4, n1, m4, m3]] *
                         dict_index_Aabcd[[n2, l4, m4, m4]]
                    end
                    for m3 in -l3:l3, m4 in -l4:l4
               ]
               )
          end
          for x in indexes
     ]
end

function Wnnll(dict_i_A::Dict{Vector{Int64},Float64},
     indexes::Vector{Vector{Int64}})

     [
          begin
               n1, n2, l3, l4 = x[1], x[2], x[3], x[4]
               16 * π^4 * sum(
                    [
                    begin
                         #println("n1 = $(n1) \t l3 = $(l3) \t m3 = $(m3)");
                         #println("n1-m3 = $(n1-m3) \t l3-m3 = $(l3-m3) \t n2-m4 = $(n2-m4) \t l4-m4 = $(l4-m4)");
                         #println("n1+m3 = $(n1+m3) \t l3+m3 = $(l3+m3) \t n2+m4 = $(n2+m4) \t l4+m4 = $(l4+m4)");
                         factorial(n1 - m3) / factorial(n1 + m3) * factorial(l3 - m3) / factorial(l3 + m3) *
                         factorial(n2 - m4) / factorial(n2 + m4) * factorial(l4 - m4) / factorial(l4 + m4) *
                         dict_i_A[[n1, l3, m3, m3]] *
                         dict_i_A[[n2, l3, m4, m3]] *
                         dict_i_A[[l4, n1, m4, m3]] *
                         dict_i_A[[n2, l4, m4, m4]]
                    end
                    for m3 in -l3:l3, m4 in -l4:l4
               ]
               )
          end
          for x in indexes
     ]
end



"""
     Wnnll(θ_max, n1::Int, n2::Int, l3::Int, l4::Int; kwargs...)
     Wnnll(dict_i_A::Dict{Vector{Int64},Float64},
          n1::Int, n2::Int, l3::Int, l4::Int)
     Wnnll(θ_max::Float64, indexes::Vector{Vector{Int64}};
          rtol=1e-2, kwargs...)
     Wnnll(dict_i_A::Dict{Vector{Int64},Float64},
          indexes::Vector{Vector{Int64}})

These 4 methods compute the following function:

```math
     W_{n_1 n_2}^{\\ell_3 \\ell_4}  = 
          16 \\pi^4 \\sum_{m_3 m_4} 
          \\frac{\\left( n_1 -m_3 \\right)! }{ \\left(n_1 +m_3 \\right)!} 
          \\frac{\\left( \\ell_3 -m_3 \\right)! }{ \\left(\\ell_3 +m_3 \\right)!}
          \\frac{\\left( n_2 -m_4 \\right)! }{ \\left(n_2 +m_4 \\right)!} 
          \\frac{\\left( \\ell_4 -m_4 \\right)! }{ \\left(\\ell_4 +m_4 \\right)!} 
          A_{n_1 \\ell_3}^{m_3 m_3} A_{n_2 \\ell_3}^{m_4 m_3} 
          A_{\\ell_4 n_1 }^{m_4 m_3} A_{n_2 \\ell_4}^{m_4 m_4}
```
where 
```math
     A_{a b}^{c d} = 
          \\int_{\\cos \\theta_{\\mathrm{max}}}^1 \\mathrm{d}x 
          P_a^c \\left( x \\right)  P_b^d \\left( x \\right)
```
The ``A_{a b}^{c d}`` are computed through `Aabcd`.

See also: [`Aabcd`](@ref)
"""
Wnnll
=#
