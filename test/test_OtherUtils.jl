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


@testset "test warning" begin
     @test "WARNING: ciao\n" == @capture_out GaPSE.warning("ciao")
end


@testset "test check_compatible_dicts" begin
     REF = Dict(
          :N => 1024::Integer, :fit_min => 0.05::Float64,
          :fit_max => 0.5::Float64, :con => true::Bool, :name => "file_name.txt"::String)

     var = 3.0
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:ma => 1.0, :con => true), "x")
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(var => 1.0))
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict("con" => true))

     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:con => 1.0), "x")
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:con => true, :N => 3.14))
     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:ma => 1.0), "x")

     @test_throws AssertionError GaPSE.check_compatible_dicts(REF, Dict(:fit_min => 1.0 + 3 * im), "x")

     @test isnothing(GaPSE.check_compatible_dicts(REF, Dict(:fit_min => 1, :fit_max => 17.43), "x"))
     @test isnothing(GaPSE.check_compatible_dicts(REF, Dict(:fit_min => 1e-2, :fit_max => 17)))
     @test isnothing(GaPSE.check_compatible_dicts(REF, Dict()))
end


@testset "test my_println_vec" begin
     vec = [x for x in 1:0.1:4]

     out_1 = @capture_out begin
          GaPSE.my_println_vec(vec, "vector"; N=8)
     end

     @test out_1 == "vector = [\n" *
                    "1.0 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , \n" *
                    "1.8 , 1.9 , 2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , \n" *
                    "2.6 , 2.7 , 2.8 , 2.9 , 3.0 , 3.1 , 3.2 , 3.3 , \n" *
                    "3.4 , 3.5 , 3.6 , 3.7 , 3.8 , 3.9 , 4.0 , \n" *
                    "];\n"

     out_2 = @capture_out begin
          GaPSE.my_println_vec(vec, "vector"; N=3)
     end
     @test out_2 == "vector = [\n" *
                    "1.0 , 1.1 , 1.2 , \n" *
                    "1.3 , 1.4 , 1.5 , \n" *
                    "1.6 , 1.7 , 1.8 , \n" *
                    "1.9 , 2.0 , 2.1 , \n" *
                    "2.2 , 2.3 , 2.4 , \n" *
                    "2.5 , 2.6 , 2.7 , \n" *
                    "2.8 , 2.9 , 3.0 , \n" *
                    "3.1 , 3.2 , 3.3 , \n" *
                    "3.4 , 3.5 , 3.6 , \n" *
                    "3.7 , 3.8 , 3.9 , \n" *
                    "4.0 , \n" *
                    "];\n"
end

@testset "test parent_directory" begin
     @test GaPSE.parent_directory("/Users/matteofoglieni/Downloads/file.txt") == "/Users/matteofoglieni/Downloads/"
     @test GaPSE.parent_directory("/Users/matteofoglieni/Downloads/") == "/Users/matteofoglieni/"
     @test GaPSE.parent_directory("/Users/matteofoglieni/Downloads") == "/Users/matteofoglieni/"
     @test GaPSE.parent_directory("/Users/matteofoglieni/") == "/Users/"
     @test GaPSE.parent_directory("/Users/matteofoglieni") == "/Users/"
     @test GaPSE.parent_directory("/Users/") == "/"
     @test GaPSE.parent_directory("/Users") == "/"

     @test GaPSE.parent_directory("matteofoglieni/Downloads/file.txt") == "matteofoglieni/Downloads/"
     @test GaPSE.parent_directory("matteofoglieni/Downloads/") == "matteofoglieni/"
     @test GaPSE.parent_directory("matteofoglieni/Downloads") == "matteofoglieni/"
     @test GaPSE.parent_directory("matteofoglieni/") == "./"
     @test GaPSE.parent_directory("matteofoglieni") == "./"
     @test GaPSE.parent_directory("file.txt") == "./"
end

@testset "test check_parent_directory" begin
     @test_throws AssertionError GaPSE.check_parent_directory("/Users/matteofoglieni/notadirectory/file.txt")
     @test_throws AssertionError GaPSE.check_parent_directory("notadirectory/file.txt")

     #@test isnothing(GaPSE.check_parent_directory("/Users/matteofoglieni/notadirectory/"))
     @test isnothing(GaPSE.check_parent_directory("notadirectory/"))
     #@test isnothing(GaPSE.check_parent_directory("/Users/"))
     @test isnothing(GaPSE.check_parent_directory("/"))
end

@testset "test return_namefile" begin
     @test_throws AssertionError GaPSE.return_namefile("/Users/matteofoglieni/Downloads/")
     @test_throws AssertionError GaPSE.return_namefile("Downloads/")
     @test_throws AssertionError GaPSE.return_namefile("./Downloads/")

     @test GaPSE.return_namefile("file") == "file"
     @test GaPSE.return_namefile("file.boh") == "file.boh"
     @test GaPSE.return_namefile("/Users/matteo.foglieni/ciao.file") == "ciao.file"
     @test GaPSE.return_namefile("matteo.foglieni/ciao.file.boh") == "ciao.file.boh"
     @test GaPSE.return_namefile("/Users/matteofoglieni/Downloads/file.txt") == "file.txt"
     @test GaPSE.return_namefile("./file.txt") == "file.txt"
     @test GaPSE.return_namefile("./file.dat") == "file.dat"
     @test GaPSE.return_namefile("file.txt") == "file.txt"
     @test GaPSE.return_namefile("file.dat") == "file.dat"
end


@testset "test check_namefile" begin
     @test_throws AssertionError GaPSE.check_namefile("/Users/matteofoglieni/Downloads/file")
     @test_throws AssertionError GaPSE.check_namefile("/Users/matteofoglieni/Downloads/file.boh")
     @test_throws AssertionError GaPSE.check_namefile("/Users/matteo.foglieni/ciao.file")
     @test_throws AssertionError GaPSE.check_namefile("matteo.foglieni/ciao.file.boh")
     @test_throws AssertionError GaPSE.check_namefile("file")
     @test_throws AssertionError GaPSE.check_namefile("file.boh")
     @test_throws AssertionError GaPSE.check_namefile("./file")
     @test_throws AssertionError GaPSE.check_namefile("./file.boh")

     @test isnothing(GaPSE.check_namefile("/Users/matteofoglieni/Downloads/file.txt"))
     @test isnothing(GaPSE.check_namefile("./file.txt"))
     @test isnothing(GaPSE.check_namefile("./file.dat"))
     @test isnothing(GaPSE.check_namefile("file.txt"))
     @test isnothing(GaPSE.check_namefile("file.dat"))
end

@testset "test check_group" begin
     @test_throws AssertionError GaPSE.check_group("prova")
     @test_throws AssertionError GaPSE.check_group("ld")
     @test_throws AssertionError GaPSE.check_group("generical")

     @test isnothing(GaPSE.check_group("LD"))
     @test isnothing(GaPSE.check_group("GNC"))
     @test isnothing(GaPSE.check_group("GNCxLD"))
     @test isnothing(GaPSE.check_group("LDxGNC"))
     @test isnothing(GaPSE.check_group("generic"))
end

@testset "test check_fileisingroup" begin
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
     @test isnothing(GaPSE.check_fileisingroup(name, "generic"))
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


##########################################################################################92


@testset "test number_to_string" begin
     @test GaPSE.number_to_string(3) == "3"
     @test GaPSE.number_to_string(-2.15) == "-2.15"
     @test GaPSE.number_to_string(2.15 * im) == "2.15im"
     @test GaPSE.number_to_string(0.0 + 2.15 * im) == "2.15im"
     @test GaPSE.number_to_string(-0.0 - 2.15 * im) == "-2.15im"
     @test GaPSE.number_to_string(-3.1415 - 2.15 * im) == "-3.1415-2.15im"
     @test GaPSE.number_to_string(3.1415 +   2.15 * im) == "3.1415+2.15im"
end
