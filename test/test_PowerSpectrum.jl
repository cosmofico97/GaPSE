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



@testset "test TwoFAST PS_multipole" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:epl => true, :pr => false, :twofast => true,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :int_s_min => 1e-4, :int_s_max => 1e4)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_TwoFAST.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_TwoFAST.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_TwoFAST.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_TwoFAST.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
     end
end


@testset "test TwoFAST print_PS_multipole first method" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:epl => true, :pr => false, :twofast => true,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :int_s_min => 1e-4, :int_s_max => 1e4)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)
               "datatest/Power_Spectrum/x"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end
end

@testset "test TwoFAST print_PS_multipole second method" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:epl => true, :pr => false, :twofast => true,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :int_s_min => 1e-4, :int_s_max => 1e4)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])


               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_TwoFAST.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_TwoFAST.txt"

               isfile(out_file) && rm(out_file)
               "datatest/Power_Spectrum/x"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end
end

@testset "test TwoFAST print_all_PS_multipole" begin
     RTOL = 1e-2

     @testset "zeros" begin
          name = "test_check_fileisingroup.txt"
          isfile(name) && rm(name)
          open(name, "w") do io
               println(io, "# boh")
               for i in 1:10
                    println(io, "$i $(2*i) $(3*i)")
               end
          end

          @test_throws AssertionError GaPSE.check_fileisingroup(name, "prova")
          @test_throws AssertionError GaPSE.check_fileisingroup(name, "LD")
          rm(name)

          open(name, "w") do io
               println(io, "# boh")
               for i in 1:10
                    println(io, "$i $(2*i) $(3*i)")
               end
               println(io, "13 14 ")
          end
          @test_throws AssertionError GaPSE.check_fileisingroup(name, "generic")
          rm(name)
     end

     kwargs_ps = Dict(:epl => true, :pr => false, :twofast => true,
          :N_left => 12, :N_right => 12,
          :p0_left => [-2.0, 1.0], :p0_right => [-2.0, 1.0],
          :int_s_min => 1e-4, :int_s_max => 1e4)

     @testset "test PS LD" begin
          L = 0
          group = "LD"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_TwoFAST.txt"
          out_file = "calc_pk_L$L" * "_TwoFAST.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS GNC" begin
          L = 0
          group = "GNC"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_TwoFAST.txt"
          out_file = "calc_pk_L$L" * "_TwoFAST.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS GNCxLD" begin
          L = 0
          group = "GNCxLD"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_TwoFAST.txt"
          out_file = "calc_pk_L$L" * "_TwoFAST.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS LDxGNC" begin
          L = 0
          group = "LDxGNC"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_TwoFAST.txt"
          out_file = "calc_pk_L$L" * "_TwoFAST.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS generic" begin
          L = 0
          group = "generic"
          input = "datatest/Power_Spectrum/map_sum_xi_LD" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_LD" * "_L$L" * "_noF_TwoFAST.txt"
          out_file = "calc_pk_L$L" * "_TwoFAST.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end
end

####################################################################################################92 




@testset "test FFTLog PS_multipole" begin
     RTOL = 1e-3
    kwargs_ps = Dict(:pr => false, :twofast => false,
         :ν => 1.5, :n_extrap_low => 500,
         :n_extrap_high => 500, :n_pad => 500)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_FFTLog.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_FFTLog.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_FFTLog.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_FFTLog.txt"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               calc_ks, calc_pks = GaPSE.PS_multipole(input; L=L, kwargs_ps...)

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
     end
end


@testset "test FFTLog print_PS_multipole first method" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:pr => false, :twofast => false,
          :ν => 1.5, :n_extrap_low => 500,
          :n_extrap_high => 500, :n_pad => 500)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)
               "datatest/Power_Spectrum/x"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               GaPSE.print_PS_multipole(input, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end
end

@testset "test FFTLog print_PS_multipole second method" begin
     RTOL = 1e-3
     kwargs_ps = Dict(:pr => false, :twofast => false,
          :ν => 1.5, :n_extrap_low => 500,
          :n_extrap_high => 500, :n_pad => 500)

     @testset "with F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])


               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_withF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_withF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_withF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end

     @testset "without F" begin
          @testset "monopole" begin
               L = 0
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end

          @testset "quadrupole" begin
               L = 2
               input = "datatest/Power_Spectrum/xi_LD_auto_doppler_noF_L$L" * ".txt"
               true_pk = "datatest/Power_Spectrum/ps_LD_auto_doppler_noF_L$L" * "_FFTLog.txt"
               out_file = "calc_pk_auto_doppler_noF_L$L" * "_FFTLog.txt"

               isfile(out_file) && rm(out_file)
               "datatest/Power_Spectrum/x"

               table = readdlm(true_pk; comments=true)
               ks = convert(Vector{Float64}, table[:, 1])
               pks = convert(Vector{Float64}, table[:, 2])

               in_table = readdlm(input; comments=true)
               in_ss = convert(Vector{Float64}, in_table[:, 1])
               in_xis = convert(Vector{Float64}, in_table[:, 2])
               GaPSE.print_PS_multipole(in_ss, in_xis, out_file; L=L, kwargs_ps...)
               calc_table = readdlm(out_file; comments=true)
               calc_ks = convert(Vector{Float64}, calc_table[:, 1])
               calc_pks = convert(Vector{Float64}, calc_table[:, 2])

               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])

               rm(out_file)
          end
     end
end

@testset "test FFTLog print_all_PS_multipole" begin
     RTOL = 1e-2

     @testset "zeros" begin
          name = "test_check_fileisingroup.txt"
          isfile(name) && rm(name)
          open(name, "w") do io
               println(io, "# boh")
               for i in 1:10
                    println(io, "$i $(2*i) $(3*i)")
               end
          end

          @test_throws AssertionError GaPSE.check_fileisingroup(name, "prova")
          @test_throws AssertionError GaPSE.check_fileisingroup(name, "LD")
          rm(name)

          open(name, "w") do io
               println(io, "# boh")
               for i in 1:10
                    println(io, "$i $(2*i) $(3*i)")
               end
               println(io, "13 14 ")
          end
          @test_throws AssertionError GaPSE.check_fileisingroup(name, "generic")
          rm(name)
     end

     kwargs_ps = Dict(:pr => false, :twofast => false,
          :ν => 1.5, :n_extrap_low => 500,
          :n_extrap_high => 500, :n_pad => 500)

     @testset "test PS LD" begin
          L = 0
          group = "LD"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_FFTLog.txt"
          out_file = "calc_pk_L$L" * "_FFTLog.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS GNC" begin
          L = 0
          group = "GNC"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_FFTLog.txt"
          out_file = "calc_pk_L$L" * "_FFTLog.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS GNCxLD" begin
          L = 0
          group = "GNCxLD"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_FFTLog.txt"
          out_file = "calc_pk_L$L" * "_FFTLog.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS LDxGNC" begin
          L = 0
          group = "LDxGNC"
          input = "datatest/Power_Spectrum/map_sum_xi_$group" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_$group" * "_L$L" * "_noF_FFTLog.txt"
          out_file = "calc_pk_L$L" * "_FFTLog.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end

     @testset "test PS generic" begin
          L = 0
          group = "generic"
          input = "datatest/Power_Spectrum/map_sum_xi_LD" * "_L$L" * "_noF.txt"
          true_pk = "datatest/Power_Spectrum/map_sum_ps_LD" * "_L$L" * "_noF_FFTLog.txt"
          out_file = "calc_pk_L$L" * "_FFTLog.txt"

          isfile(out_file) && rm(out_file)

          table = readdlm(true_pk; comments=true)
          ks = convert(Vector{Float64}, table[:, 1])
          all_pks = [convert(Vector{Float64}, col)
                     for col in eachcol(table[:, 2:end])]


          in_table = readdlm(input; comments=true)
          in_ss = convert(Vector{Float64}, in_table[:, 1])
          in_xis = convert(Vector{Float64}, in_table[:, 2])
          GaPSE.print_all_PS_multipole(input, out_file, group; L=L, kwargs_ps...)
          calc_table = readdlm(out_file; comments=true)
          calc_ks = convert(Vector{Float64}, calc_table[:, 1])
          calc_all_pks = [convert(Vector{Float64}, col)
                          for col in eachcol(calc_table[:, 2:end])]

          @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(ks, calc_ks)])
          for (calc_pks, pks) in zip(all_pks, calc_all_pks)
               @test all([isapprox(t, c; rtol=RTOL) for (t, c) in zip(pks, calc_pks)])
          end
          rm(out_file)
     end
end
