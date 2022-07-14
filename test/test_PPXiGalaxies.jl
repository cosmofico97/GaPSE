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




@testset "test map_ξ_PPGalaxies_multipole" begin
     RTOL = 1e-3
     kwargs_xis_PP = Dict(
          :pr => false, :enhancer => 1e8,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :N_log => 100,
     )

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_withF_L$L" * ".txt"

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               calc_ss, calc_xis = GaPSE.map_ξ_PPGalaxies_multipole(COSMO,
                    10 .^ range(0, 3, length=300); use_windows=true,
                    L=L, kwargs_xis_PP...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
          end

          @testset "quadrupole" begin
               L = 2
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_withF_L$L" * ".txt"

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               calc_ss, calc_xis = GaPSE.map_ξ_PPGalaxies_multipole(COSMO,
                    10 .^ range(0, 3, length=300); use_windows=true,
                    L=L, kwargs_xis_PP...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
          end

          @testset "hexadecapole" begin
               L = 4
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_withF_L$L" * ".txt"

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               calc_ss, calc_xis = GaPSE.map_ξ_PPGalaxies_multipole(COSMO,
                    10 .^ range(0, 3, length=300); use_windows=true,
                    L=L, kwargs_xis_PP...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_noF_L$L" * ".txt"

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               calc_ss, calc_xis = GaPSE.map_ξ_PPGalaxies_multipole(COSMO,
                    10 .^ range(0, 3, length=300); use_windows=false,
                    L=L, kwargs_xis_PP...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
          end

          @testset "quadrupole" begin
               L = 2
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_noF_L$L" * ".txt"

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               calc_ss, calc_xis = GaPSE.map_ξ_PPGalaxies_multipole(COSMO,
                    10 .^ range(0, 3, length=300); use_windows=false,
                    L=L, kwargs_xis_PP...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
          end

          @testset "hexadecapole" begin
               L = 4
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_noF_L$L" * ".txt"

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               calc_ss, calc_xis = GaPSE.map_ξ_PPGalaxies_multipole(COSMO,
                    10 .^ range(0, 3, length=300); use_windows=false,
                    L=L, kwargs_xis_PP...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])
          end
     end
end



##########################################################################################92




@testset "test print_map_ξ_PPGalaxies_multipole" begin
     RTOL = 1e-3
     kwargs_xis_PP = Dict(
          :pr => false, :enhancer => 1e8,
          :μ_atol => 0.0, :μ_rtol => 1e-2,
          :N_log => 100,
     )

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               name = "calc_xi_ppgalaxies_withF_L$L" * ".txt"
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_withF_L$L" * ".txt"

               isfile(name) && rm(name)

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_map_ξ_PPGalaxies_multipole(COSMO, name,
                    10 .^ range(0, 3, length=300); use_windows=true,
                    L=L, kwargs_xis_PP...)

               calc_table = readdlm(true_xi; comments=true)
               calc_ss = convert(Vector{Float64}, calc_table[:, 1])
               calc_xis = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

               rm(name)
          end

          @testset "quadrupole" begin
               L = 2
               name = "calc_xi_ppgalaxies_withF_L$L" * ".txt"
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_withF_L$L" * ".txt"

               isfile(name) && rm(name)

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_map_ξ_PPGalaxies_multipole(COSMO, name,
                    10 .^ range(0, 3, length=300); use_windows=true,
                    L=L, kwargs_xis_PP...)

               calc_table = readdlm(true_xi; comments=true)
               calc_ss = convert(Vector{Float64}, calc_table[:, 1])
               calc_xis = convert(Vector{Float64}, calc_table[:, 2])

               println()

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

               rm(name)
          end

          @testset "hexadecapole" begin
               L = 4
               name = "calc_xi_ppgalaxies_withF_L$L" * ".txt"
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_withF_L$L" * ".txt"

               isfile(name) && rm(name)

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_map_ξ_PPGalaxies_multipole(COSMO, name,
                    10 .^ range(0, 3, length=300); use_windows=true,
                    L=L, kwargs_xis_PP...)

               calc_table = readdlm(true_xi; comments=true)
               calc_ss = convert(Vector{Float64}, calc_table[:, 1])
               calc_xis = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

               rm(name)
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               name = "calc_xi_ppgalaxies_noF_L$L" * ".txt"
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_noF_L$L" * ".txt"

               isfile(name) && rm(name)

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_map_ξ_PPGalaxies_multipole(COSMO, name,
                    10 .^ range(0, 3, length=300); use_windows=false,
                    L=L, kwargs_xis_PP...)

               calc_table = readdlm(true_xi; comments=true)
               calc_ss = convert(Vector{Float64}, calc_table[:, 1])
               calc_xis = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

               rm(name)
          end

          @testset "quadrupole" begin
               L = 2
               name = "calc_xi_ppgalaxies_noF_L$L" * ".txt"
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_noF_L$L" * ".txt"

               isfile(name) && rm(name)

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_map_ξ_PPGalaxies_multipole(COSMO, name,
                    10 .^ range(0, 3, length=300); use_windows=false,
                    L=L, kwargs_xis_PP...)

               calc_table = readdlm(true_xi; comments=true)
               calc_ss = convert(Vector{Float64}, calc_table[:, 1])
               calc_xis = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

               rm(name)
          end

          @testset "hexadecapole" begin
               L = 4
               name = "calc_xi_ppgalaxies_noF_L$L" * ".txt"
               true_xi = "datatest/PPXiGalaxies/xi_ppgalaxies_noF_L$L" * ".txt"

               isfile(name) && rm(name)

               table = readdlm(true_xi; comments=true)
               ss = convert(Vector{Float64}, table[:, 1])
               xis = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_map_ξ_PPGalaxies_multipole(COSMO, name,
                    10 .^ range(0, 3, length=300); use_windows=false,
                    L=L, kwargs_xis_PP...)

               calc_table = readdlm(true_xi; comments=true)
               calc_ss = convert(Vector{Float64}, calc_table[:, 1])
               calc_xis = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ss, calc_ss)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(xis, calc_xis)])

               rm(name)
          end
     end
end

