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

module GaPSE

using TwoFAST # Licence: MIT "Expat" (o GPL ?)

using Dierckx # Licence: BSD
using HCubature, QuadGK, WignerSymbols # Licence: MIT "Expat"
using LegendrePolynomials, AssociatedLegendrePolynomials # Licence: MIT "Expat"
using SpecialFunctions, Trapz, LsqFit # Licence: MIT
using GridInterpolations  # Licence: MIT "Expat"
using ProgressMeter, Printf  # Licence: MIT "Expat"

using Test, Documenter, DelimitedFiles  # Licence: MIT "Expat"


const BRAND = """
        ###############
        #    GaPSE    #
        ############### \n#"""

####### DO NOT MODIFY THESE NAMES AND/OR THEIR ORDER IF YOU ARE NOT SURE! ###############92

const NAMES_F_MAP = ["x", "mu", "F", "F_error"]
const NAMES_PS = ["k (h/Mpc)", "P (Mpc/h)^3"]
const NAMES_BACKGROUND = ["z", "proper time [Gyr]", "conf. time [Mpc]", "H [1/Mpc]",
     "comov. dist.", "ang.diam.dist.", "lum. dist.", "comov.snd.hrz.",
     "(.)rho_g", "(.)rho_b", "(.)rho_cdm", "(.)rho_lambda", "(.)rho_ur",
     "(.)rho_crit", "gr.fac. D", "gr.fac. f"]
const VALID_GROUPS = ["LD", "GNC", "GNCxLD", "LDxGNC", "generic"]
const LENGTH_VALID_GROUPS = [18, 27, 22, 22, nothing]

include("OtherUtils.jl")
include("MathUtils.jl")
include("Wllnn.jl")
include("WindowF.jl")
include("WindowFIntegrated.jl")
include("BackgroundData.jl")
include("CosmoParams.jl")
include("CosmoUtils.jl")
include("IPSTools.jl")
include("XiMatter.jl")
include("Cosmology.jl")
include("PPXiMatter.jl")


##################### PERTURBED LUMINOSITY DISTANCE #############################


include("LD_AutoCorrelations/LD_AutoDoppler.jl")
include("LD_AutoCorrelations/LD_AutoLensing.jl")
include("LD_AutoCorrelations/LD_AutoLocalGP.jl")
include("LD_AutoCorrelations/LD_AutoIntegratedGP.jl")

include("LD_CrossCorrelations/LD_LensingDoppler.jl")
include("LD_CrossCorrelations/LD_DopplerLocalGP.jl")
include("LD_CrossCorrelations/LD_DopplerIntegratedGP.jl")
include("LD_CrossCorrelations/LD_LensingLocalGP.jl")
include("LD_CrossCorrelations/LD_LensingIntegratedGP.jl")
include("LD_CrossCorrelations/LD_LocalGPIntegratedGP.jl")

include("PPDoppler.jl")


####################################### GALAXY NUMBER COUNTS #############################92


include("GNC_AutoCorrelations/GNC_AutoNewtonian.jl")
include("GNC_AutoCorrelations/GNC_AutoDoppler.jl")
include("GNC_AutoCorrelations/GNC_AutoLensing.jl")
include("GNC_AutoCorrelations/GNC_AutoLocalGP.jl")
include("GNC_AutoCorrelations/GNC_AutoIntegratedGP.jl")

include("GNC_CrossCorrelations/GNC_NewtonianDoppler.jl")
include("GNC_CrossCorrelations/GNC_NewtonianLensing.jl")
include("GNC_CrossCorrelations/GNC_NewtonianLocalGP.jl")
include("GNC_CrossCorrelations/GNC_NewtonianIntegratedGP.jl")
include("GNC_CrossCorrelations/GNC_LensingDoppler.jl")
include("GNC_CrossCorrelations/GNC_DopplerLocalGP.jl")
include("GNC_CrossCorrelations/GNC_DopplerIntegratedGP.jl")
include("GNC_CrossCorrelations/GNC_LensingLocalGP.jl")
include("GNC_CrossCorrelations/GNC_LensingIntegratedGP.jl")
include("GNC_CrossCorrelations/GNC_LocalGPIntegratedGP.jl")


##################### GALAXY NUMBER COUNTS x LUMINOSITY DISTANCE #########################92


include("GNCxLD_CrossCorrelations/GNCxLD_NewtonianDoppler.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_NewtonianLensing.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_NewtonianLocalGP.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_NewtonianIntegratedGP.jl")

include("GNCxLD_CrossCorrelations/GNCxLD_DopplerDoppler.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_DopplerLensing.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_DopplerLocalGP.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_DopplerIntegratedGP.jl")

include("GNCxLD_CrossCorrelations/GNCxLD_LensingDoppler.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_LensingLensing.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_LensingLocalGP.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_LensingIntegratedGP.jl")

include("GNCxLD_CrossCorrelations/GNCxLD_LocalGPDoppler.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_LocalGPLensing.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_LocalGPLocalGP.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_LocalGPIntegratedGP.jl")

include("GNCxLD_CrossCorrelations/GNCxLD_IntegratedGPDoppler.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_IntegratedGPLensing.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_IntegratedGPLocalGP.jl")
include("GNCxLD_CrossCorrelations/GNCxLD_IntegratedGPIntegratedGP.jl")


##########################################################################################92


include("Dicts.jl")


include("LD_XiMultipoles.jl")
include("LD_SumXiMultipoles.jl")
include("GNC_XiMultipoles.jl")
include("GNC_SumXiMultipoles.jl")
include("GNCxLD_XiMultipoles.jl")
include("GNCxLD_SumXiMultipoles.jl")
include("LDxGNC_XiMultipoles.jl")
include("LDxGNC_SumXiMultipoles.jl")
include("PowerSpectrum.jl")


function parameters_used(io::IO, cosmo::Cosmology; logo::Bool=true)
     logo && println(io, BRAND)
     println(io, "# The Cosmology considered had the following paremeters:\n# ")
     println(io, "# - Matter Power Spectrum input file: \"$(cosmo.file_ips)\"")
     println(io, "# - Background data input file: \"$(cosmo.file_data)\"")
     println(io, "# - F window function input file: \"$(cosmo.file_windowF)\"")
     println(io, "# - Integrated F window function input file: \"$(cosmo.file_IWF)\"")
     println(io, "#")

     println(io, "# - Basic CosmoParams considered: ")
     println(io, "#\t z_min = $(cosmo.params.z_min) \t z_max = $(cosmo.params.z_max)")
     println(io, "#\t ??_max = $(cosmo.params.??_max) [rad] \t h_0 = $(cosmo.params.h_0)")
     println(io, "#\t ??_b = $(cosmo.params.??_b) \t " *
                 "??_cdm = $(cosmo.params.??_cdm) \t ??_M0 = $(cosmo.params.??_M0)")
     println(io, "#\t b = $(cosmo.params.b) \t " *
                 "f_evo = $(cosmo.params.????_evo) \t s_b = $(cosmo.params.s_b)")
     println(io, "#")

     println(io, "# - CosmoParams about the Input Power Spectrum: ")
     my_println_dict(io, cosmo.params.IPS; pref="#\t ", N=2)
     println(io, "#")

     println(io, "# - CosmoParams about the Input Power Spectrum Tools: ")
     my_println_dict(io, cosmo.params.IPSTools; pref="#\t ", N=3)
     println(io, "#")

     println(io, "# - CosmoParams about the Integrated Window Function F: ")
     my_println_dict(io, cosmo.params.WFI; pref="#\t ", N=3)
     println(io, "#")

     println(io, "# - Computed quantities: ")
     println(io, "# \t effective redshift z_eff = $(cosmo.z_eff) ")
     println(io, "# \t comoving s_min = $(cosmo.s_min) Mpc/h_0")
     println(io, "# \t comoving s_max = $(cosmo.s_max) Mpc/h_0")
     println(io, "# \t comoving s_eff = $(cosmo.s_eff) Mpc/h_0")
     println(io, "# \t Volume of the survey V_survey = $(cosmo.volume) Mpc^3/h_0^3")
     println(io, "# \t ??_0 = $(cosmo.tools.??_0)")
     println(io, "# \t ??_1 = $(cosmo.tools.??_1)")
     println(io, "# \t ??_2 = $(cosmo.tools.??_2)")
     println(io, "# \t ??_3 = $(cosmo.tools.??_3)")
     println(io, "# \t (where ??_i = \\int_{k_{min}}^{k_{max}}\\frac{dq}{2 ??^2} q^{2-i} P(q))")
     println(io, "# ")
end

parameters_used(cosmo::Cosmology) = parameters_used(stdout, cosmo)


"""
     parameters_used(io::IO, cosmo::Cosmology)
     parameters_used(cosmo::Cosmology) = parameters_used(stdout, cosmo)

Writes in the `io` stream all the information concerning the input Cosmology `cosmo`.

See also: [`Cosmology`](@ref)
"""
parameters_used


end # module
