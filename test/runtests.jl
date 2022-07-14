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


using GaPSE, Test
using Dierckx, DelimitedFiles, QuadGK, Suppressor


const FILE_F_MAP = "datatest/F_REFERENCE_pi2.txt"
const FILE_PS = "datatest/file_pk.txt"
const FILE_IF_MAP = "datatest/IntegrF_REFERENCE_pi2.txt"
const FILE_ILN = "datatest/table_Iln.txt"
const FILE_BACKGROUND = "datatest/WideA_ZA_background.dat"

const Z_MIN = 0.05
const Z_MAX = 0.20
const S_MIN = 148.1920001343431
const S_MAX = 571.7022420911966
const S_EFF = 435.37470960794167
const Z_EFF = 0.15045636097417317
const VOLUME = 3.845366167993746e8
const HUBBLE_0 = 1e5 / 299792458.0

include("TEST_DATA.jl")

##########################################################################################92


@testset "test MathUtils" begin
     include("test_MathUtils.jl")
end

@testset "test OtherUtils" begin
     include("test_OtherUtils.jl")
end

@testset "test CosmoUtils" begin
     include("test_CosmoUtils.jl")
end

@testset "test IPSTools" begin
     include("test_IPSTools.jl")
end

@testset "test BackgroundData" begin
     include("test_BackgroundData.jl")
end

@testset "test CosmoParams" begin
     include("test_CosmoParams.jl")
end

@testset "test Cosmology" begin
     include("test_Cosmology.jl")
end

@testset "test WindowF" begin
     include("test_WindowF.jl")
end

@testset "test WindowFIntegrated" begin
     include("test_WindowFIntegrated.jl")
end

@testset "test Dicts" begin
     include("test_Dicts.jl")
end

@testset "test PowerSpectrum" begin
     include("test_PowerSpectrum.jl")
end



################################### COSMOLOGY IMPLEMENTATION #############################92


const PARAMS = GaPSE.CosmoParams(Z_MIN, Z_MAX, π / 2.0;
     Ω_b=0.0489, Ω_cdm=0.251020, h_0=0.7, s_lim=1e-2,
     IPS_opts=Dict(
          :fit_left_min => 1e-6, :fit_left_max => 3e-6,
          :fit_right_min => 1e1, :fit_right_max => 2e1,
     ),
     IPSTools_opts=Dict(
          :N => 1024, :fit_min => 0.05, :fit_max => 0.5,
          :con => true, :k_min => 1e-8, :k_max => 10.0,
     ),
     WFI_opts=Dict(:llim => 0.0, :rlim => Inf, :N => 1000,
          :trap => true, :rtol => 1e-2, :atol => 0.0,
          :ss_start => 0.0, :ss_step => 21.768735478453323,
          :ss_stop => 0.0)
)

const COSMO = GaPSE.Cosmology(PARAMS, FILE_BACKGROUND, FILE_PS, FILE_F_MAP, FILE_IF_MAP)

#=
common_kwargs = Dict(
     :pr => false,
     :use_windows => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_log => 100,
);

spec_effect = [
     "auto_lensing", "auto_integratedgp",
     "lensing_doppler", "doppler_lensing",
     "doppler_integratedgp", "integratedgp_doppler",
     "lensing_localgp", "localgp_lensing",
     "lensing_integratedgp", "integratedgp_lensing",
     "localgp_integratedgp", "integratedgp_localgp",
];

specific_kwargs = [effect ∈ spec_effect ? Dict(
     :en => 1e12, :N_χs => 30) : nothing for effect in GaPSE.GR_EFFECTS_LD]

joint_kwargs = [isnothing(spec) ? common_kwargs : merge(common_kwargs, spec)
                for spec in specific_kwargs];
=#

dict_L_dir = Dict(0 => "monopoles", 1 => "dipoles", 2 => "tripoles",
     3 => "quadrupoles", 4 => "pentapoles")

KWARGS_LD = Dict(
     :pr => true,
     :use_windows => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_χs => 100, :N_χs_2 => 30,
     :N_log => 100,
);

KWARGS_GNC = Dict(
     :pr => true,
     :use_windows => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_χs => 100, :N_χs_2 => 30,
     :N_log => 100,
);

KWARGS_GNCxLD = Dict(
     :pr => true,
     :use_windows => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_χs => 100, :N_χs_2 => 30,
     :N_log => 100,
);

KWARGS_LDxGNC = Dict(
     :pr => true,
     :use_windows => false,
     :enhancer => 1e8, :N_μs => 30,
     :μ_atol => 0.0, :μ_rtol => 1e-2,
     :N_χs => 100, :N_χs_2 => 30,
     :N_log => 100,
);


################################### TEST PLANE-PARALLEL APPROXIMATIONS#### ###############92

@testset "test PPDoppler" begin
     include("test_PPDoppler.jl")
end

@testset "test XiMatter" begin
     include("test_XiMatter.jl")
end

@testset "test PPXiGalaxies" begin
     include("test_PPXiGalaxies.jl")
end


################################### TEST LUMINOSITY DISTANCE PERTURBATIONS ###############92


@testset "test LD_XiMultipoles" begin
     include("test_LD_XiMultipoles.jl")
end

@testset "test LD_SumXiMultipoles" begin
     include("test_LD_SumXiMultipoles.jl")
end


##############################


table_xis_LD_L0_noF = readdlm("datatest/LD_SumXiMultipoles/map_sum_xi_LD_L0_noF.txt",
     comments=true);
ss_LD_L0_noF = convert(Vector{Float64}, table_xis_LD_L0_noF[:, 1]);
xis_sum_LD_L0_noF = convert(Vector{Float64}, table_xis_LD_L0_noF[:, 2]);
all_xis_LD_L0_noF = [convert(Vector{Float64}, table_xis_LD_L0_noF[:, i]) for i in 3:18];

@testset "test LD_AutoDoppler" begin
     include("test_LD_AutoCorrelations/test_LD_AutoDoppler.jl")
end

@testset "test LD_AutoIntegratedGP" begin
     include("test_LD_AutoCorrelations/test_LD_AutoIntegratedGP.jl")
end

@testset "test LD_AutoLocalGP" begin
     include("test_LD_AutoCorrelations/test_LD_AutoLocalGP.jl")
end

@testset "test LD_AutoLensing" begin
     include("test_LD_AutoCorrelations/test_LD_AutoLensing.jl")
end


##############################



@testset "test LD_DopplerLensing" begin
     include("test_LD_CrossCorrelations/test_LD_DopplerLensing.jl")
end

@testset "test LD_DopplerLocalGP" begin
     include("test_LD_CrossCorrelations/test_LD_DopplerLocalGP.jl")
end

@testset "test LD_DopplerIntegratedGP" begin
     include("test_LD_CrossCorrelations/test_LD_DopplerIntegratedGP.jl")
end

@testset "test LD_LensingIntegratedGP" begin
     include("test_LD_CrossCorrelations/test_LD_LensingIntegratedGP.jl")
end

@testset "test LD_LensingLocalGP" begin
     include("test_LD_CrossCorrelations/test_LD_LensingLocalGP.jl")
end

@testset "test LD_LocalGPIntegratedGP" begin
     include("test_LD_CrossCorrelations/test_LD_LocalGPIntegratedGP.jl")
end


################################### TEST RELATIVISTIC GALAXY NUMBER COUNTS ###############92


@testset "test GNC_XiMultipoles" begin
     include("test_GNC_XiMultipoles.jl")
end

@testset "test GNC_SumXiMultipoles" begin
     include("test_GNC_SumXiMultipoles.jl")
end


##############################

table_xis_GNC_L0_noF = readdlm("datatest/GNC_SumXiMultipoles/map_sum_xi_GNC_L0_noF.txt",
     comments=true);
ss_GNC_L0_noF = convert(Vector{Float64}, table_xis_GNC_L0_noF[:, 1]);
xis_sum_GNC_L0_noF = convert(Vector{Float64}, table_xis_GNC_L0_noF[:, 2]);
all_xis_GNC_L0_noF = [convert(Vector{Float64}, table_xis_GNC_L0_noF[:, i]) for i in 3:27];


@testset "test GNC_AutoNewton" begin
     include("test_GNC_AutoCorrelations/test_GNC_AutoNewton.jl")
end

@testset "test GNC_AutoDoppler" begin
     include("test_GNC_AutoCorrelations/test_GNC_AutoDoppler.jl")
end

@testset "test GNC_AutoIntegratedGP" begin
     include("test_GNC_AutoCorrelations/test_GNC_AutoIntegratedGP.jl")
end

@testset "test GNC_AutoLocalGP" begin
     include("test_GNC_AutoCorrelations/test_GNC_AutoLocalGP.jl")
end

@testset "test GNC_AutoLensing" begin
     include("test_GNC_AutoCorrelations/test_GNC_AutoLensing.jl")
end


##############################

@testset "test GNC_NewtonDoppler" begin
     include("test_GNC_CrossCorrelations/test_GNC_NewtonDoppler.jl")
end

@testset "test GNC_NewtonLensing" begin
     include("test_GNC_CrossCorrelations/test_GNC_NewtonDoppler.jl")
end

@testset "test GNC_NewtonLocalGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_NewtonLocalGP.jl")
end

@testset "test GNC_NewtonIntegratedGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_NewtonIntegratedGP.jl")
end

@testset "test GNC_DopplerLensing" begin
     include("test_GNC_CrossCorrelations/test_GNC_DopplerLensing.jl")
end

@testset "test GNC_DopplerLocalGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_DopplerLocalGP.jl")
end

@testset "test GNC_DopplerIntegratedGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_DopplerIntegratedGP.jl")
end

@testset "test GNC_LensingIntegratedGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_LensingIntegratedGP.jl")
end

@testset "test GNC_LensingLocalGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_LensingLocalGP.jl")
end

@testset "test GNC_LocalGPIntegratedGP" begin
     include("test_GNC_CrossCorrelations/test_GNC_LocalGPIntegratedGP.jl")
end


##### TEST RELATIVISTIC GALAXY NUMBER COUNTS X LUMINOSITY DISTANCE PERT. and viceversa ###92


@testset "test GNCxLD_XiMultipoles" begin
     include("test_GNCxLD_XiMultipoles.jl")
end

@testset "test LDxGNC_XiMultipoles" begin
     include("test_LDxGNC_XiMultipoles.jl")
end

@testset "test GNCxLD_SumXiMultipoles" begin
     include("test_GNCxLD_SumXiMultipoles.jl")
end

@testset "test LDxGNC_SumXiMultipoles" begin
     include("test_LDxGNC_SumXiMultipoles.jl")
end

##############################

table_xis_GNCxLD_L0_noF = readdlm("datatest/GNCxLD_SumXiMultipoles/map_sum_xi_GNCxLD_L0_noF.txt",
     comments=true);
ss_GNCxLD_L0_noF = convert(Vector{Float64}, table_xis_GNCxLD_L0_noF[:, 1]);
xis_sum_GNCxLD_L0_noF = convert(Vector{Float64}, table_xis_GNCxLD_L0_noF[:, 2]);
all_xis_GNCxLD_L0_noF = [convert(Vector{Float64}, table_xis_GNCxLD_L0_noF[:, i]) for i in 3:22];

table_xis_LDxGNC_L0_noF = readdlm("datatest/GNCxLD_SumXiMultipoles/map_sum_xi_LDxGNC_L0_noF.txt",
     comments=true);
ss_LDxGNC_L0_noF = convert(Vector{Float64}, table_xis_LDxGNC_L0_noF[:, 1]);
xis_sum_LDxGNC_L0_noF = convert(Vector{Float64}, table_xis_LDxGNC_L0_noF[:, 2]);
all_xis_LDxGNC_L0_noF = [convert(Vector{Float64}, table_xis_LDxGNC_L0_noF[:, i]) for i in 3:22];


@testset "test GNCxLD_NewtonDoppler" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_NewtonDoppler.jl")
end

@testset "test GNCxLD_NewtonLensing" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_NewtonLensing.jl")
end

@testset "test GNCxLD_NewtonLocalGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_NewtonLocalGP.jl")
end

@testset "test GNCxLD_NewtonIntegratedGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_NewtonIntegratedGP.jl")
end

############

@testset "test GNCxLD_DopplerDoppler" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_DopplerDoppler.jl")
end

@testset "test GNCxLD_DopplerLensing" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_DopplerLensing.jl")
end

@testset "test GNCxLD_DopplerLocalGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_DopplerLocalGP.jl")
end

@testset "test GNCxLD_DopplerIntegratedGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_DopplerIntegratedGP.jl")
end

############

@testset "test GNCxLD_LensingDoppler" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LensingDoppler.jl")
end

@testset "test GNCxLD_LensingLensing" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LensingLensing.jl")
end

@testset "test GNCxLD_LensingLocalGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LensingLocalGP.jl")
end

@testset "test GNCxLD_LensingIntegratedGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LensingIntegratedGP.jl")
end

############

@testset "test GNCxLD_LocalGPDoppler" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LocalGPDoppler.jl")
end

@testset "test GNCxLD_LocalGPLensing" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LocalGPLensing.jl")
end

@testset "test GNCxLD_LocalGPLocalGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LocalGPLocalGP.jl")
end

@testset "test GNCxLD_LocalGPIntegratedGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_LocalGPIntegratedGP.jl")
end

############

@testset "test GNCxLD_IntegratedGPDoppler" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_IntegratedGPDoppler.jl")
end

@testset "test GNCxLD_IntegratedGPLensing" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_IntegratedGPLensing.jl")
end

@testset "test GNCxLD_IntegratedGPLocalGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_IntegratedGPLocalGP.jl")
end

@testset "test GNCxLD_IntegratedGPIntegratedGP" begin
     include("test_GNCxLD_CrossCorrelations/test_GNCxLD_IntegratedGPIntegratedGP.jl")
end
