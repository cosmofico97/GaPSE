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

@testset "test xi auto_lensing L = 0" begin
     table = readdlm("datatest/monopoles/xi_lensing_L0.txt"; comments = true)
     ss = convert(Vector{Float64}, table[:, 1])
     xis = convert(Vector{Float64}, table[:, 2])

     calc_xis = [GaPSE.ξ_multipole(COSMO.s_eff, s, "auto_lensing", COSMO;
          L = 0, N_μs = 20, N_χs = 30, Δχ_min = 1e-4, enhancer = 1e10, use_windows = false) for s in ss]

     @test all([isapprox(xi, calc_xi, rtol = 1e-3) for (xi, calc_xi) in
                zip(xis[ss.>0.5], calc_xis[ss.>0.5])])
end
