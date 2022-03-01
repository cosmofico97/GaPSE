var documenterSearchIndex = {"docs":
[{"location":"BackgroundData/","page":"Background Data","title":"Background Data","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"BackgroundData/#Background-functions","page":"Background Data","title":"Background functions","text":"","category":"section"},{"location":"BackgroundData/","page":"Background Data","title":"Background Data","text":"GaSPE.ℋ0 \nGaPSE.z_eff","category":"page"},{"location":"AutoCorrelations/","page":"Auto Correlations","title":"Auto Correlations","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"AutoCorrelations/#Auto-Correlation-functions","page":"Auto Correlations","title":"Auto-Correlation functions","text":"","category":"section"},{"location":"AutoCorrelations/","page":"Auto Correlations","title":"Auto Correlations","text":"GaPSE.ξ_Doppler\nGaPSE.ξ_Lensing\nGaPSE.ξ_LocalGP\nGaPSE.ξ_IntegratedGP\nGaPSE.integrand_ξ_Lensing\nGaPSE.integrand_ξ_IntegratedGP","category":"page"},{"location":"AutoCorrelations/#GaPSE.ξ_Doppler","page":"Auto Correlations","title":"GaPSE.ξ_Doppler","text":" ξ_Doppler(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64\n\nReturn the Doppler auto-correlation function, defined as follows:\n\nxi^v_parallelv_parallel (s_1 s_2 costheta) \n= D_1 D_2 f_1 f_2 mathcalH_1 mathcalH_2 mathcalR_1 mathcalR_2 \n(J_00 I^0_0(s) + J_02I^0_2(s) + J_04I^0_4(s) + J_20I^2_0(s))\n\nwhere D_1 = D(s_1), D_2 = D(s_2) and so on, mathcalH = a H,  y = costheta = hatmathbfs_1 dot hatmathbfs_2 and  the J coefficients are given by:\n\nbeginalign*\n    J_00 (s_1 s_2 y)  = frac145 (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \n    J_02 (s_1 s_2 y)  = frac263 (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \n    J_04 (s_1 s_2 y)  = frac1105 (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \n    J_20 (s_1 s_2 y)  = frac13 y s^2\nendalign*\n\nInputs\n\nP1::Point and P2::Point: Point where the CF has to be calculated; they contain all the  data of interest needed for this calculus (comoving distance, growth factor and so on)\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nSee also: Point, Cosmology\n\n\n\n\n\n","category":"function"},{"location":"AutoCorrelations/#GaPSE.ξ_Lensing","page":"Auto Correlations","title":"GaPSE.ξ_Lensing","text":" ξ_Lensing(s1, s2, y, cosmo::Cosmology;\n      en::Float64 = 1e6,\n      Δχ_min::Float64 = 1e-3,\n      N_χs::Integer = 100) :: Float64\n\nReturn the Lensing auto-correlation function  xi^kappakappa (s_1 s_2 costheta), defined as follows:\n\nxi^kappakappa (s_1 s_2 costheta) = \nint_0^s_1 mathrmd chi_1 int_0^s_2 mathrmd chi_2 \nfrac12\nfrac\n     mathcalH_0^4 Omega_ mathrmM0^2 D_1 D_2 (chi_1 - s_1)(chi_2 - s_2)\n\n     s_1 s_2 a(chi_1) a(chi_2) \n(J_00  I^0_0(chi) + J_02  I^0_2(chi) + \n     J_31  I^3_1(chi) + J_22  I^2_2(chi))\n\nwhere D_1 = D(chi_1), D_2 = D(chi_2) and so on, mathcalH = a H,  chi = sqrtchi_1^2 + chi_2^2 - 2chi_1chi_2costheta,  y = costheta = hatmathbfs_1 dot hatmathbfs_2)  and the J coefficients are given by \n\nbeginalign*\n    J_00  = - frac3 chi_1^2 chi_2^24 chi^4 (y^2 - 1) \n               (8 y (chi_1^2 + chi_2^2) - 9 chi_1 chi_2 y^2 - 7 chi_1 chi_2) \n    J_02  = - frac3 chi_1^2 chi_2^22 chi^4 (y^2 - 1)\n               (4 y (chi_1^2 + chi_2^2) - 3 chi_1 chi_2 y^2 - 5 chi_1 chi_2) \n    J_31  = 9 y chi^2 \n    J_22  = frac9 chi_1 chi_24 chi^4\n                2 (chi_1^4 + chi_2^4) (7 y^2 - 3) \n                 - 16 y chi_1 chi_2 (chi_1^2 + chi_2^2) (y^2+1) \n               + chi_1^2 chi_2^2 (11 y^4 + 14 y^2 + 23)\nendalign*\n\nThe computation is made applying trapz (see the  Trapz Julia package) to the integrand function integrand_ξ_Lensing.\n\nInputs\n\ns1 and s2: comovign distances where the function must be evaluated\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nen::Float64 = 1e6: just a float number used in order to deal better  with small numbers;\nΔχ_min::Float64 = 1e-6 : when Deltachi = sqrtchi_1^2 + chi_2^2 - 2  chi_1 chi_2 y to 0^+, some I_ell^n term diverges, but the overall parenthesis has a known limit:\n   lim_chito0^+ (J_00  I^0_0(chi) + J_02  I^0_2(chi) + \n        J_31  I^3_1(chi) + J_22  I^2_2(chi)) = \n        frac415  (5  sigma_2 + frac23  σ_0 s_1^2  chi_2^2)\nSo, when it happens that chi  Deltachi_mathrmmin, the function considers this limit as the result of the parenthesis instead of calculating it in the normal way; it prevents computational divergences.\nN_χs::Integer = 100: number of points to be used for sampling the integral along the ranges (0, s1) (for χ1) and (0, s1) (for χ2); it has been checked that with N_χs ≥ 50 the result is stable.\n\nSee also: integrand_ξ_Lensing, integrand_on_mu_Lensing integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"AutoCorrelations/#GaPSE.ξ_LocalGP","page":"Auto Correlations","title":"GaPSE.ξ_LocalGP","text":" ξ_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64\n\nReturn the local gravitational potential auto-correlation function,  defined as follows:\n\nxi^phiphi (s_1 s_2 costheta) = \n     frac9 mathcalH_0^4 Omega_M0^2 D(s_1) D(s_2)s^44 a(s_1) a(s_2)\n     (1 + mathcalR_1 + mathcalR_2 + mathcalR_1mathcalR_2)\n     tildeI^4_0(s)\n\nwhere D_1 = D(s_1), D_2 = D(s_2) and so on, mathcalH = a H,  y = costheta = hatmathbfs_1 dot hatmathbfs_2 and:\n\ntildeI^4_0 (s) = int_0^infty fracmathrmdq2pi^2 \n          q^2  P(q)  fracj_0(q s) - 1(q s)^4\n\nInputs\n\nP1::Point and P2::Point: Point where the CF has to be calculated; they contain all the  data of interest needed for this calculus (comoving distance, growth factor and so on)\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nSee also: Point, Cosmology\n\n\n\n\n\n","category":"function"},{"location":"AutoCorrelations/#GaPSE.ξ_IntegratedGP","page":"Auto Correlations","title":"GaPSE.ξ_IntegratedGP","text":" ξ_IntegratedGP(s1, s2, y, cosmo::Cosmology; \n      en::Float64 = 1e10,\n      N_χs::Integer = 100) :: Float64\n\nReturn the integrated gravitational potential auto-correlation function  xi^intphiintphi(s_1 s_2 costheta), defined as follows:\n\nxi^intphiintphi (s_1 s_2 costheta) = \n     int_0^s_1 mathrmd chi_1 int_0^s_2mathrmd chi_2 \n     J_40(s_1 s_2 y chi_1 chi_2)  tildeI^4_0(chi)\n\nwhere chi = sqrtchi_1^2 + chi_2^2 - 2  chi_1  chi_2  y, y = costheta = hatmathbfs_1 dot hatmathbfs_2 and:\n\nbeginalign*\n     J_40(s_1 s_2 y chi_1 chi_2) = \n          frac\n                    9 mathcalH_0^4 Omega_M0^2 D(chi_1) D(chi_2) chi^4\n              a(chi_1) a(chi_2) s_1 s_2 \n          (s_2 mathcalH(chi_2) mathcalR(s_2) (f(chi_2)-1) - 1) \n          (s_1 mathcalH(chi_1) mathcalR(s_1) (f(chi_1)-1) - 1) 5pt\n     tildeI^4_0 (s) = int_0^infty fracmathrmdq2pi^2 \n          q^2  P(q)  fracj_0(q s) - 1(q s)^4\nendaling*\n\nand P(q) is the input power spectrum.\n\nThe computation is made applying trapz (see the  Trapz Julia package) to the integrand function integrand_ξ_Lensing.\n\nInputs\n\ns1 and s2: comovign distances where the function must be evaluated\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nen::Float64 = 1e10: just a float number used in order to deal better  with small numbers.\nN_χs::Integer = 100: number of points to be used for sampling the integral along the ranges (0, s1) (for χ1) and (0, s1) (for χ2); it has been checked that with N_χs ≥ 50 the result is stable.\n\nSee also: integrand_ξ_IntegratedGP, integrand_on_mu_IntegratedGP integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"AutoCorrelations/#GaPSE.integrand_ξ_Lensing","page":"Auto Correlations","title":"GaPSE.integrand_ξ_Lensing","text":" integrand_ξ_Lensing(\n      IP1::Point, IP2::Point,\n      P1::Point, P2::Point,\n      y, cosmo::Cosmology;\n      Δχ_min::Float64 = 1e-4) :: Float64\n\nReturn the integrand of the Lensing auto-correlation function  xi^kappakappa (s_1 s_2 costheta), i.e. the function  f(s_1 s_2 y chi_1 chi_2) defined as follows:  \n\nf(s_1 s_2 y chi_1 chi_2) = \nfrac12\nfrac\n     mathcalH_0^4 Omega_ mathrmM0^2 D_1 D_2 (chi_1 - s_1)(chi_2 - s_2)\n\n     s_1 s_2 a(chi_1) a(chi_2) \n(J_00  I^0_0(chi) + J_02  I^0_2(chi) + \n     J_31  I^3_1(chi) + J_22  I^2_2(chi))\n\nwhere D_1 = D(chi_1), D_2 = D(chi_2) and so on, mathcalH = a H,  chi = sqrtchi_1^2 + chi_2^2 - 2chi_1chi_2costheta,  y = costheta = hatmathbfs_1 dot hatmathbfs_2)  and the J coefficients are given by \n\nbeginalign*\n    J_00  = - frac3 chi_1^2 chi_2^24 chi^4 (y^2 - 1) \n               (8 y (chi_1^2 + chi_2^2) - 9 chi_1 chi_2 y^2 - 7 chi_1 chi_2) \n    J_02  = - frac3 chi_1^2 chi_2^22 chi^4 (y^2 - 1)\n               (4 y (chi_1^2 + chi_2^2) - 3 chi_1 chi_2 y^2 - 5 chi_1 chi_2) \n    J_31  = 9 y chi^2 \n    J_22  = frac9 chi_1 chi_24 chi^4\n                2 (chi_1^4 + chi_2^4) (7 y^2 - 3) \n                 - 16 y chi_1 chi_2 (chi_1^2 + chi_2^2) (y^2+1) \n               + chi_1^2 chi_2^2 (11 y^4 + 14 y^2 + 23)\nendalign*\n\nInputs\n\nIP1::Point and IP2::Point: Point inside the integration limits, placed  at comoving distance χ1 and χ2 respectively.\nP1::Point and P2::Point: extreme Point of the integration, placed  at comoving distance s1 and s2 respectively.\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nΔχ_min::Float64 = 1e-6 : when Deltachi = sqrtchi_1^2 + chi_2^2 - 2  chi_1 chi_2 y to 0^+, some I_ell^n term diverges, but the overall parenthesis has a known limit:\n   lim_chito0^+ (J_00  I^0_0(chi) + J_02  I^0_2(chi) + \n        J_31  I^3_1(chi) + J_22  I^2_2(chi)) = \n        frac415  (5  sigma_2 + frac23  σ_0 s_1^2  chi_2^2)\nSo, when it happens that chi  Deltachi_mathrmmin, the function considers this limit as the result of the parenthesis instead of calculating it in the normal way; it prevents computational divergences.\n\nSee also: ξ_Lensing, integrand_on_mu_Lensing integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"AutoCorrelations/#GaPSE.integrand_ξ_IntegratedGP","page":"Auto Correlations","title":"GaPSE.integrand_ξ_IntegratedGP","text":" integrand_ξ_IntegratedGP(IP1::Point, IP2::Point,\n      P1::Point, P2::Point,\n      y, cosmo::Cosmology) :: Float64\n\nReturn the integrand of the integrated gravitational potential  auto-correlation function xi^intphiintphi (s_1 s_2 costheta),  i.e. the function f(s_1 s_2 y chi_1 chi_2) defined as follows:  \n\nf(s_1 s_2 y chi_1 chi_2) = J_40(s_1 s_2 y chi_1 chi_2) tildeI^4_0(chi)\n\nwhere chi = sqrtchi_1^2 + chi_2^2 - 2  chi_1  chi_2  y, y = costheta = hatmathbfs_1 dot hatmathbfs_2 and:\n\nbeginsplit\n     J_40(s_1 s_2 y chi_1 chi_2)  = \n          frac\n                    9 mathcalH_0^4 Omega_M0^2 D(chi_1) D(chi_2) chi^4\n              a(chi_1) a(chi_2) s_1 s_2 \n          (s_2 mathcalH(chi_2) mathcalR(s_2) (f(chi_2)-1) - 1) \n          (s_1 mathcalH(chi_1) mathcalR(s_1) (f(chi_1)-1) - 1)5pt\n     tildeI^4_0 (s) = int_0^infty fracmathrmdq2pi^2 \n          q^2  P(q)  fracj_0(q s) - 1(q s)^4\nendsplit\n\nInputs\n\nIP1::Point and IP2::Point: Point inside the integration limits, placed  at comoving distance χ1 and χ2 respectively.\nP1::Point and P2::Point: extreme Point of the integration, placed  at comoving distance s1 and s2 respectively.\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nSee also: ξ_IntegratedGP, integrand_on_mu_IntegratedGP integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"F_evaluation/","page":"The F window function","title":"The F window function","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"F_evaluation/#The-F-window-function","page":"The F window function","title":"The F window function","text":"","category":"section"},{"location":"F_evaluation/","page":"The F window function","title":"The F window function","text":"GaPSE.F\nGaPSE.integrand_F\nGaPSE.F_map","category":"page"},{"location":"F_evaluation/#GaPSE.F","page":"The F window function","title":"GaPSE.F","text":" F(x, μ, θ_max; tolerance=1e-8) :: Tuple{Float64, Float64}\n\nThe function F(xmu theta_mathrmmax), defined as follows:\n\nbeginsplit\nF(xmu theta_mathrmmax) =  4pi \n    int_0^theta_mathrmmax mathrmdtheta_1 int_0^pi mathrmd theta  \n     Thetaleft(frac\n        x cos theta + cos theta_1sqrtx^1+2+2xmu - \n        cos(theta_mathrmmax) \n        right) \n     Theta(mu-cos(theta+theta_1)) \n    Theta(cos(theta - theta_1)-mu) \n    fracsinthetasintheta_1\n        sqrt(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2\nendsplit\n\ntolerance is a parameter needed in case of small negative denominator: the Heaviside theta function mathematically prevent that  mathrmden=(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2 becomes negative, but computationally might happen that mathrmden results as a very small negative number (for instance -1.2368946523-18); in this case tolerance solve the problem, returning 0 if 0-mathrmden mathrmtolerance.\n\nThe double integral is performed with hcubature function from the Julia Package HCubature; rtol, atol and all the kwargs insert into F  are directly transferred to hcubature. \n\nPAY ATTENTION: do not set too small atol and rtol, or the computation can easily become overwhelming! \n\nSee also: F_map, integrand_F\n\n\n\n\n\n","category":"function"},{"location":"F_evaluation/#GaPSE.integrand_F","page":"The F window function","title":"GaPSE.integrand_F","text":"integrand_F(θ, θ_1, x, μ, θ_max; tolerance=1e-8) :: Float64\n\nReturn the integrand of the function F(xmu theta_mathrmmax), i.e the  function f(xmu theta theta_1 theta_mathrmmax):\n\nbeginsplit\nf(xmu theta theta_1 theta_mathrmmax) = \n     Thetaleft( frac\n        x cos theta + cos theta_1sqrtx^1+2+2xmu - \n        cos(theta_mathrmmax) \n        right) \n     times  Theta(mu-cos(theta+theta_1))   times \n     quad Theta(cos(theta - theta_1)-mu)  times \n    frac4pi sinthetasintheta_1\n        sqrt(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2\nendsplit\n\nbeginequation\nF(xmu theta_mathrmmax) = int_0^theta_mathrmmax \n        mathrmdtheta_1 int_0^pi mathrmd theta \n         f(xmu theta theta_1 theta_mathrmmax)\nendequation\n\ntolerance is a parameter needed in case of small negative denominator: the Heaviside theta function mathematically prevent that  mathrmden=(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2 becomes negative, but computationally might happen that mathrmden results as a very small negative number (for instance -1.2368946523-18); in this case tolerance solve the problem, returning 0 if 0-mathrmden mathrmtolerance\n\nSee also: F, F_map\n\n\n\n\n\n","category":"function"},{"location":"F_evaluation/#GaPSE.F_map","page":"The F window function","title":"GaPSE.F_map","text":" F_map(x_step::Float64 = 0.01, μ_step::Float64 = 0.01;\n      out = \"data/F_map.txt\", x1 = 0, x2 = 3, μ1 = -1, μ2 = 1, kwargs...)\n F_map(xs::Vector{Float64}, μs::Vector{Float64};\n      out = \"data/F_map.txt\", kwargs...)\n\nEvaluate the window function F(xmu theta_mathrmmax) in a rectangual grid  of mu and x values.\n\nIn the first method you specify start, stop and step for x and μ manually, while with thr second one you pass the values (through a vector )you want to calculate  the function in.\n\nSee also: F_map, integrand_F\n\n\n\n\n\n","category":"function"},{"location":"PowerSpectrum/","page":"Power Spectrum Multipoles","title":"Power Spectrum Multipoles","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"PowerSpectrum/#The-Power-Spectrum-multipole-computation","page":"Power Spectrum Multipoles","title":"The Power Spectrum multipole computation","text":"","category":"section"},{"location":"PowerSpectrum/","page":"Power Spectrum Multipoles","title":"Power Spectrum Multipoles","text":"GaPSE.PS_multipole","category":"page"},{"location":"PowerSpectrum/#GaPSE.PS_multipole","page":"Power Spectrum Multipoles","title":"GaPSE.PS_multipole","text":" PS_multipole(in::String, int_s_min = 0.0, int_s_max = 1000.0; \n      L::Integer = 0, N::Integer = 1024, pr::Bool = true)\n PS_multipole(\n      f_in::Union{Function,Dierckx.Spline1D}, int_s_min = 0.0, int_s_max = 1000.0; \n      L::Integer = 0, N::Integer = 1024, pr::Bool = true\n      ) :: Tuple{Vector{Float64}, Vector{Float64}}\n\nReturn the L-order multipole from the input function f_in, through the following Fast Fourier Transform and the effective redshift approximation:\n\nP_L(k) = frac2 L + 1A^ (-i)^L  phi(s_mathrmeff) int_0^infty \n        mathrmd s  s^2  j_L(ks)  f_mathrmin(s)  \n        quad  A^ = 8  pi^2\n\nThis expression can be easily obtained from the standard one:\n\nbeginsplit\n    P_L(k) = frac2 L + 1A (-i)^L  \n        int_0^infty mathrmd s_1  s_1^2 \n        int_0^infty mathrmd s  s^2 \n        int_-1^+1 mathrmd mu \n        j_L(ks)  xi(s_1 s mu)  phi(s_1)  phi(s_2) \n        mathcalL_L(mu) Fleft(fracss_1 mu right) \n        mathrmwith  s_2 = s_2(s_1 s μ) = sqrts_1^2 + s^2 + 2s_1smu\n         \n         quad A(s_mathrmmax s_mathrmmin theta_mathrmmax) \n        = 2  pi  V\nendsplit\n\nwith the definition\n\nf_mathrmin(s_1 s) =  int_-1^+1 mathrmd mu \n        xi(s_1 s mu)  phi(s_2) \n        mathcalL_L(mu)  Fleft(fracss_1 mu right)\n\nand the application of the effective redshift approximation.\n\nSee also: V_survey, A, A_prime\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/","page":"Input Power Spectrum Tools","title":"Input Power Spectrum Tools","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"IPSTools/#Tool-functions","page":"Input Power Spectrum Tools","title":"Tool functions","text":"","category":"section"},{"location":"IPSTools/","page":"Input Power Spectrum Tools","title":"Input Power Spectrum Tools","text":"GaPSE.spline_F\nGaPSE.I00\nGaPSE.I20\nGaPSE.I40\nGaPSE.I02\nGaPSE.I22\nGaPSE.I31\nGaPSE.I13\nGaPSE.I11\nGaPSE.σ_0\nGaPSE.σ_1\nGaPSE.σ_2\nGaPSE.σ_3 \nGaPSE.ϕ\nGaPSE.W\nGaPSE.V_survey\nGaPSE.A\nGaPSE.A_prime\nGaPSE.s\nGaPSE.μ\nGaPSE.s2\nGaPSE.y","category":"page"},{"location":"IPSTools/#GaPSE.spline_F","page":"Input Power Spectrum Tools","title":"GaPSE.spline_F","text":" spline_F(x, μ, str::WindowF)) :: Float64\n\nReturn the 2-dim spline value of F in the given (x,μ), where F is defined in the input WindowF. The spline is obtained through the interpolate function of the  GridInterpolations Julia package.\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.ϕ","page":"Input Power Spectrum Tools","title":"GaPSE.ϕ","text":" ϕ(s, s_min, s_max) :: Float64\n\nRadial part of the survey window function. Return 1.0 if is true that s_mathrmmin le s le s_mathrmmax and 0.0 otherwise.\n\nIn this software we made the assuption that the survey window function can be separated into a radial and angular part, i.e.:\n\n     phi(mathbfs) = phi(s)  W(hats)\n\nSee also: W\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.W","page":"Input Power Spectrum Tools","title":"GaPSE.W","text":" W(θ; θ_max = θ_MAX) :: Float64\n\nAngular part of the survey window function. Return 1.0 if is true that 00 leq theta le theta_mathrmmax and 0.0 otherwise. It is implicitly assumed an azimutal simmetry of the survey.\n\nIn this software we made the assuption that the survey window function can be separated into a radial and angular part, i.e.:\n\n     phi(mathbfs) = phi(s)  W(hats)\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.V_survey","page":"Input Power Spectrum Tools","title":"GaPSE.V_survey","text":" V_survey(s_min = s_MIN, s_max = s_MAX, θ_max = θ_MAX) :: Float64\n\nReturn the volume of a survey with azimutal simmetry, i.e.:\n\nbeginsplit\n    V(s_mathrmmax s_mathrmmin theta_mathrmmax) =  C_mathrmup - C_mathrmdown + TC \n    C_mathrmup = fracpi3 s_mathrmmax^3  \n        (1 - costheta_mathrmmax)^2  (2 + costheta_mathrmmax) \n    C_mathrmdown = fracpi3 s_mathrmmin^3  \n        (1 - costheta_mathrmmax)^2  (2 + costheta_mathrmmax) \n    TC = fracpi3 (s_mathrmmax^2 + s_mathrmmin^2 + \n        s_mathrmmax s_mathrmmin)   (s_mathrmmax - s_mathrmmin) \n        costheta_mathrmmax sin^2theta_mathrmmax\nendsplit\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.A","page":"Input Power Spectrum Tools","title":"GaPSE.A","text":" A(s_min, s_max, θ_max) :: Float64\n\nReturn the Power Spectrum multipole normalization coefficient A, i.e.:\n\n     A(s_mathrmmax s_mathrmmin theta_mathrmmax)= 2  pi  \n     V(s_mathrmmax s_mathrmmin theta_mathrmmax)\n\nwhere V(s_mathrmmax s_mathrmmin theta_mathrmmax) is the  survey volume.\n\nPay attention: this is NOT used for the normalization of PS, see instead A_prime\n\nSee also: V_survey\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.A_prime","page":"Input Power Spectrum Tools","title":"GaPSE.A_prime","text":" A_prime :: Float64\n\nIt's the Power Spectrum multipole normalization coefficient A^, i.e.:\n\n     A^ = frac3  A (s_mathrmmax^3 - s_mathrmmin^3) = 8 pi^2\n\nSee also: A, V_survey\n\n\n\n\n\n","category":"constant"},{"location":"IPSTools/#GaPSE.s","page":"Input Power Spectrum Tools","title":"GaPSE.s","text":" s(s1, s2, y) :: Float64\n\nReturn the value s = sqrts_1^2 + s_2^2 - 2  s_1  s_2  y\n\nSee also: μ, s2, y\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.μ","page":"Input Power Spectrum Tools","title":"GaPSE.μ","text":" μ(s1, s2, y) :: Float64\n\nReturn the value mu=hatmathbfs_1dothatmathbfs, defined as:\n\nmu = mu(s_1 s_2 y) = fracy  s_2 - s_1s(s_1 s_2 y) \nquad s(s_1 s_2 y) = sqrts_1^2 + s^2 - 2  s_1  s_2  y\n\nwith y=costheta=hatmathbfs_1dothatmathbfs and where s is  obtained from the function s\n\nSee also: s, s2, y\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.s2","page":"Input Power Spectrum Tools","title":"GaPSE.s2","text":" s2(s1, s, μ) :: Float64\n\nReturn the value s_2 = sqrts_1^2 + s^2 + 2  s_1  s  mu\n\nSee also: s, μ, y\n\n\n\n\n\n","category":"function"},{"location":"IPSTools/#GaPSE.y","page":"Input Power Spectrum Tools","title":"GaPSE.y","text":" y(s1, s, μ) :: Float64\n\nReturn the value y=costheta, defined as:\n\ny = y(s_1 s mu) = fracmu  s + s_1s2(s_1 s mu) \nquad s_2 = sqrts_1^2 + s^2 + 2  s_1  s  mu\n\nwith mu=hatmathbfs_1dothatmathbfs_2 and  where s_2 is btained from the function s2\n\nSee also: s, μ, s2\n\n\n\n\n\n","category":"function"},{"location":"","page":"Introduction","title":"Introduction","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"#GaPSE.jl-:-a-Galaxy-Power-Spectrum-Estimator","page":"Introduction","title":"GaPSE.jl : a Galaxy Power Spectrum Estimator","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This is the documentation of GaPSE.jl package, an implementation of a Galaxy Power Spectrum Estimator written in Julia.","category":"page"},{"location":"#Documentation","page":"Introduction","title":"Documentation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The documentation was built using Documenter.jl.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using Dates # hide\nprintln(\"Documentation built on $(now()) using Julia $(VERSION).\") # hide","category":"page"},{"location":"#Contents","page":"Introduction","title":"Contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"#Index","page":"Introduction","title":"Index","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"CrossCorrelations/","page":"CrossCorrelations","title":"CrossCorrelations","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"CrossCorrelations/#Cross-Correlation-functions","page":"CrossCorrelations","title":"Cross-Correlation functions","text":"","category":"section"},{"location":"CrossCorrelations/","page":"CrossCorrelations","title":"CrossCorrelations","text":"GaPSE.ξ_Lensing_Doppler\nGaPSE.ξ_Doppler_Lensing\nGaPSE.ξ_Doppler_LocalGP\nGaPSE.ξ_LocalGP_Doppler\nGaPSE.ξ_Doppler_IntegratedGP\nGaPSE.ξ_IntegratedGP_Doppler\nGaPSE.ξ_Lensing_LocalGP\nGaPSE.ξ_LocalGP_Lensing,\nGaPSE.ξ_Lensing_IntegratedGP\nGaPSE.ξ_IntegratedGP_Lensing\nGaPSE.ξ_LocalGP_IntegratedGP\nGaPSE.ξ_IntegratedGP_LocalGP","category":"page"},{"location":"CrossCorrelations/#GaPSE.ξ_Lensing_Doppler","page":"CrossCorrelations","title":"GaPSE.ξ_Lensing_Doppler","text":" ξ_Lensing_Doppler(s1, s2, y, cosmo::Cosmology;\n      en::Float64 = 1e6, N_χs::Integer = 100):: Float64\n\nReturn the Lensing-Doppler cross-correlation function  xi^v_parallelkappa (s_1 s_2 costheta), defined as follows:\n\nxi^v_parallelkappa (s_1 s_2 costheta) = \n     mathcalH_0^2 Omega_M0 D(s_2) f(s_2) mathcalH(s_2) mathcalR(s_2) \n     int_0^s_1 mathrmd chi_1 \n     frac D(chi_1) (chi_1 - s_1) a(chi_1) s_1 \n     left(\n          J_00 I^0_0(Deltachi_1) + J_02 I^0_2(Deltachi_1) \n          + J_04 I^0_4(Deltachi_1) + J_20 I^2_0(Deltachi_1)\n     right)\n\nwhere mathcalH = a H,  Deltachi_1= sqrtchi_1^2 + s_2^2 - 2 chi_1 s_2 costheta,  y = costheta = hatmathbfs_1 dot hatmathbfs_2)  and the J coefficients are given by:\n\nbeginalign*\n     J_00  = frac115(chi_1^2 y + chi_1(4 y^2 - 3) s_2 - 2 y s_2^2) \n     J_02  = frac142 Deltachi_1^2 \n          (4 chi_1^4 y + 4 chi_1^3 (2 y^2 - 3) s_2 + chi_1^2 y (11 - 23 y^2) s_2^2 + \n          chi_1 (23 y^2 - 3) s_2^3 - 8 y s_2^4) \n     J_04  = frac170 Deltachi_1^2\n          (2 chi_1^4 y + 2 chi_1^3 (2y^2 - 3) s_2 - chi_1^2 y (y^2 + 5) s_2^2 + \n          chi_1 (y^2 + 9) s_2^3 - 4 y s_2^4) \n     J_20  = y Deltachi_1^2\nendalign*\n\nThe computation is made applying trapz (see the  Trapz Julia package) to the integrand function integrand_ξ_Lensing_Doppler.\n\nInputs\n\ns1 and s2: comovign distances where the function must be evaluated\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nen::Float64 = 1e6: just a float number used in order to deal better  with small numbers;\nΔχ_min::Float64 = 1e-6 : when Deltachi = sqrtchi_1^2 + chi_2^2 - 2  chi_1 chi_2 y to 0^+, some I_ell^n term diverges, but the overall parenthesis has a known limit:\n   lim_chito0^+ (J_00  I^0_0(chi) + J_02  I^0_2(chi) + \n        J_31  I^3_1(chi) + J_22  I^2_2(chi)) = \n        frac415  (5  sigma_2 + frac23  σ_0 s_1^2  chi_2^2)\nSo, when it happens that chi  Deltachi_mathrmmin, the function considers this limit as the result of the parenthesis instead of calculating it in the normal way; it prevents computational divergences.\nN_χs::Integer = 100: number of points to be used for sampling the integral along the ranges (0, s1) (for χ1) and (0, s1) (for χ2); it has been checked that with N_χs ≥ 50 the result is stable.\n\nSee also: integrand_ξ_Lensing_Doppler, int_on_mu_Lensing_Doppler integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"CrossCorrelations/#GaPSE.ξ_Doppler_LocalGP","page":"CrossCorrelations","title":"GaPSE.ξ_Doppler_LocalGP","text":" ξ_Doppler_LocalGP(P1::Point, P2::Point, y, cosmo::Cosmology) :: Float64\n\nReturn the Doppler-LocalGP cross-correlation function, defined as follows:\n\nxi^v_parallelphi (s_1 s_2 costheta) = \n     frac32 a(s_2) mathcalH(s_1) f(s_1) D(s_1)\n     mathcalR(s_1) mathcalH_0^2 Omega_M0 D(s_2)\n     (1 + mathcalR(s_2)) (s_1 - s_2costheta) s^2 I^3_1(s)\n\nwhere mathcalH = a H, y = costheta = hatmathbfs_1 dot hatmathbfs_2 and :\n\nI^n_l(s) = int_0^infty fracmathrmdq2pi^2 q^2  P(q)  fracj_l(qs)(q s)^n\n\nInputs\n\nP1::Point and P2::Point: Point where the CF has to be calculated; they contain all the  data of interest needed for this calculus (comoving distance, growth factor and so on)\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nSee also: Point, Cosmology\n\n\n\n\n\n","category":"function"},{"location":"CrossCorrelations/#GaPSE.ξ_Doppler_IntegratedGP","page":"CrossCorrelations","title":"GaPSE.ξ_Doppler_IntegratedGP","text":" ξ_Doppler_IntegratedGP(s1, s2, y, cosmo::Cosmology;\n      en::Float64 = 1e6, N_χs::Integer = 100):: Float64\n\nReturn the Doppler-LocalGP cross-correlation function  xi^v_parallelintphi (s_1 s_2 costheta), defined as follows:\n\nxi^v_parallelintphi (s_1 s_2 costheta) = \n     3 mathcalH(s_1) f(s_1) D(s_1) mathcalH_0^2 Omega_M0 mathcalR(s_1) \n     int_0^s_2 mathrmdchi_2   J_31   I^3_1(chi)\n\nwhere mathcalH = a H,  chi = sqrts_1^2 + chi_2^2 - 2 s_1 chi_2 costheta,  y = costheta = hatmathbfs_1 dot hatmathbfs_2)  and:\n\nJ_31 = \n     fracD(chi_2) (s_1 - chi_2 costheta)a(chi_2) chi^2 \n     left(\n          - frac1s_2 + mathcalR(s_2) mathcalH(chi_2) (f(chi_2) - 1)\n     right)\n\nThe computation is made applying trapz (see the  Trapz Julia package) to the integrand function integrand_ξ_Doppler_IntegratedGP.\n\nInputs\n\ns1 and s2: comovign distances where the function must be evaluated\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nen::Float64 = 1e6: just a float number used in order to deal better  with small numbers;\nN_χs::Integer = 100: number of points to be used for sampling the integral along the ranges (0, s1) (for χ1) and (0, s1) (for χ2); it has been checked that with N_χs ≥ 50 the result is stable.\n\nSee also: integrand_ξ_Doppler_IntegratedGP, int_on_mu_Doppler_IntegratedGP integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"CrossCorrelations/#GaPSE.ξ_Lensing_LocalGP","page":"CrossCorrelations","title":"GaPSE.ξ_Lensing_LocalGP","text":" ξ_Lensing_LocalGP(s1, s2, y, cosmo::Cosmology;\n      en::Float64 = 1e6, N_χs::Integer = 100):: Float64\n\nReturn the Lensing-LocalGP cross-correlation function  xi^kappa phi (s_1 s_2 costheta), defined as follows:\n\nxi^kappa phi (s_1 s_2 costheta) = \n     frac\n          9 mathcalH_0^4 Omega_M0^2 D(s_2) (1 + mathcalR(s_2)) s_2\n     4 a(s_2) s_1 \n     int_0^s_1 mathrmdchi_1 fracD(chi_1)(s_1 - chi_1) a(chi_1)\n     left( J_31 I^3_1(Deltachi_1) +  J_22 I^2_2(Deltachi_1) right)\n\nwhere mathcalH = a H,  Deltachi_1 = sqrtchi_1^2 + s_2^2 - 2 chi_1 s_2costheta,  y = costheta = hatmathbfs_1 dot hatmathbfs_2)  and the J coefficients are given by \n\nbeginalign*\n     J_31  = -2 y Deltachi_1^2 \n     J_22  = chi_1 s_2 (1 - y^2)\nendalign*\n\nThe computation is made applying trapz (see the  Trapz Julia package) to the integrand function integrand_ξ_Lensing_LocalGP.\n\nInputs\n\ns1 and s2: comovign distances where the function must be evaluated\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nen::Float64 = 1e6: just a float number used in order to deal better  with small numbers;\nN_χs::Integer = 100: number of points to be used for sampling the integral along the ranges (0, s1) (for χ1) and (0, s1) (for χ2); it has been checked that with N_χs ≥ 50 the result is stable.\n\nSee also: integrand_ξ_Lensing_LocalGP, int_on_mu_Lensing_LocalGP integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"},{"location":"CrossCorrelations/#GaPSE.ξ_LocalGP_IntegratedGP","page":"CrossCorrelations","title":"GaPSE.ξ_LocalGP_IntegratedGP","text":" ξ_LocalGP_IntegratedGP(s1, s2, y, cosmo::Cosmology;\n      en::Float64 = 1e6, N_χs::Integer = 100):: Float64\n\nReturn the LocalGP-IntegratedGP cross-correlation function  xi^v_parallelint phi (s_1 s_2 costheta), defined as follows:\n\nxi^v_parallelint phi (s_1 s_2 costheta) = \n     frac9 mathcalH_0^4 Omega_M0^2 D(s_1) (mathcalR(s_1) +1)2 a(s_1) \n     int_0^s_2 mathrmdchi_2 fracD(chi_2) Deltachi_2^4 a(chi_2)\n     left(\n          mathcalH(chi_2)( f(chi_2) - 1) mathcalR(s_2) - frac1s_2\n     right) tildeI^4_0(Deltachi_2)\n\nwhere mathcalH = a H,  Deltachi_2 = sqrts_1^2 + chi_2^2 - 2 s_1 chi_2 costheta,  y = costheta = hatmathbfs_1 dot hatmathbfs_2).\n\nThe computation is made applying trapz (see the  Trapz Julia package) to the integrand function integrand_ξ_LocalGP_IntegratedGP.\n\nInputs\n\ns1 and s2: comovign distances where the function must be evaluated\ny: the cosine of the angle between the two points P1 and P2\ncosmo::Cosmology: cosmology to be used in this computation\n\nOptional arguments\n\nen::Float64 = 1e6: just a float number used in order to deal better  with small numbers;\nN_χs::Integer = 100: number of points to be used for sampling the integral along the ranges (0, s1) (for χ1) and (0, s1) (for χ2); it has been checked that with N_χs ≥ 50 the result is stable.\n\nSee also: integrand_ξ_LocalGP_IntegratedGP, int_on_mu_LocalGP_IntegratedGP integral_on_mu, ξ_multipole\n\n\n\n\n\n","category":"function"}]
}
