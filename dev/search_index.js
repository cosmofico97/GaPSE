var documenterSearchIndex = {"docs":
[{"location":"Auto_doppler/","page":"Doppler Auto-Correlation functions","title":"Doppler Auto-Correlation functions","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"Auto_doppler/#Doppler-Auto-Correlation-functions","page":"Doppler Auto-Correlation functions","title":"Doppler Auto-Correlation functions","text":"","category":"section"},{"location":"Auto_doppler/","page":"Doppler Auto-Correlation functions","title":"Doppler Auto-Correlation functions","text":"GaPSE.ξ_doppler\nGaSPE.integrand_on_mu_doppler\nGaPSE.integral_on_mu_doppler\nGaPSE.map_integral_on_mu_doppler\nGaPSE.print_map_int_on_mu_doppler","category":"page"},{"location":"Auto_doppler/#GaPSE.ξ_doppler","page":"Doppler Auto-Correlation functions","title":"GaPSE.ξ_doppler","text":" ξ_doppler(s1, s2, y; enhancer = 1, tol = 1)\n\nReturn the Doppler auto-correlation function, defined as follows:\n\nbeginequation\n    xi^v_parallelv_parallel (s_1 s_2 costheta) \n    = D_1 D_2 f_1 f_2 mathcalH_1 mathcalH_2 mathcalR_1 mathcalR_2 \n    (J_00 I^0_0(s) + J_02I^0_2(s) + J_04I^0_4(s) + J_20I^2_0(s))\nendequation\n\nwhere D_1 = D(s_1), D_2 = D(s_2) and so on, mathcalH = a H and  the J coefficients are given by (with y = costheta):\n\nbeginalign\n    J_00 (s_1 s_2 y)  = frac145 (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \n    J_02 (s_1 s_2 y)  = frac263 (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \n    J_04 (s_1 s_2 y)  = frac1105 (y^2 s_1 s_2 - 2y(s_1^2 + s_2^2) + 3s_1 s_2) \n    J_20 (s_1 s_2 y)  = frac13 y s^2\nendalign\n\nenhancer is just a float number used in order to deal better with small numbers; the returned value is actually enhancer * xi_doppler, where xi_doppler is the true value calculated as shown before.\n\ntol is a flaot number that set the minimum distance of interest between the two galaxies placed at mathbfs_1 and mathbfs_2; the output is set to 0.0 if it is true that:\n\nbeginsplit\n     s = s(s_1 s_2 y) = sqrts_1^2 + s_2^2 - 2  s_1  s_2  y  ge mathrmtol \n\n     mathbfwhere   y=costheta = hatmathbfs_1dothatmathbfs_2\nendsplit\n\n\n\n\n\n","category":"function"},{"location":"Auto_lensing/","page":"Lensing Auto-Correlation functions","title":"Lensing Auto-Correlation functions","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"Auto_lensing/#Lensing-Auto-Correlation-functions","page":"Lensing Auto-Correlation functions","title":"Lensing Auto-Correlation functions","text":"","category":"section"},{"location":"Auto_lensing/","page":"Lensing Auto-Correlation functions","title":"Lensing Auto-Correlation functions","text":"GaPSE.integand_ξ_lensing\nGaPSE.ξ_lensing\nGaSPE.integrand_on_mu_lensing\nGaPSE.integral_on_mu_lensing\nGaPSE.map_integral_on_mu_lensing\nGaPSE.print_map_int_on_mu_lensing","category":"page"},{"location":"Auto_lensing/#GaPSE.ξ_lensing","page":"Lensing Auto-Correlation functions","title":"GaPSE.ξ_lensing","text":" ξ_lensing(s1, s2, y; tol = 0.5, enhancer = 1, Δχ_min = 1e-6, kwargs...)\n\nReturn the Lensing Auto-correlation function  xi^kappakappa (s_1 s_2 costheta) , defined as follows:\n\nbeginequation\n    xi^kappakappa (s_1 s_2 costheta) = \n    int_0^s_1 mathrmd chi_1 int_0^s_2 mathrmd chi_2 \n     frac12\n     frac\n          mathcalH_0^4 Omega_ mathrmM0^2 D_1 D_2 (chi_1 - s_1)(chi_2 - s_2)\n     \n          s_1 s_2 a(chi_1) a(chi_2) \n     (J_00  I^0_0(chi) + J_02  I^0_2(chi) + \n          J_31  I^3_1(chi) + J_22  I^2_2(chi))\nendequation\n\nwhere D_1 = D(chi_1), D_2 = D(chi_2) and so on, mathcalH = a H,  chi = sqrtchi_1^2 + chi_2^2 - 2chi_1chi_2costheta  and the J coefficients are given by (with y = costheta)\n\nbeginalign\n    J_00  = - frac3 chi_1^2 chi_2^24 chi^4 (y^2 - 1) \n               (8 y (chi_1^2 + chi_2^2) - 9 chi_1 chi_2 y^2 - 7 chi_1 chi_2) \n    J_02  = - frac3 chi_1^2 chi_2^22 chi^4 (y^2 - 1)\n               (4 y (chi_1^2 + chi_2^2) - 3 chi_1 chi_2 y^2 - 5 chi_1 chi_2) \n    J_31  = 9 y chi^2 \n    J_22  = frac9 chi_1 chi_24 chi^4\n               2 chi_1^4 (7 y^2 - 3) - 16chi_1^3chi_2y(y^2+1) \n               + chi_1^2 chi_2^2 (11 y^4 + 14 y^2 + 23) - 16 chi_1 chi_2^3 y (y^2 + 1) \n               + 2chi_2^4(7y^2-3)\nendalign\n\nThe computation is made applying hcubature (see the  Hcubature Julia package) to the integrand function integrand_ξ_lensing.\n\nOptional arguments\n\ntol = 0.5 : \n\nSee also: integrand_ξ_lensing\n\n\n\n\n\n","category":"function"},{"location":"Background_functions/","page":"Background functions","title":"Background functions","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"Background_functions/#Background-functions","page":"Background functions","title":"Background functions","text":"","category":"section"},{"location":"Background_functions/","page":"Background functions","title":"Background functions","text":"","category":"page"},{"location":"Tool_functions/","page":"Tool functions","title":"Tool functions","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"Tool_functions/#Tool-functions","page":"Tool functions","title":"Tool functions","text":"","category":"section"},{"location":"Tool_functions/","page":"Tool functions","title":"Tool functions","text":"GaPSE.spline_F\nGaPSE.I00\nGaPSE.I20\nGaPSE.I40\nGaPSE.I02\nGaPSE.I22\nGaPSE.I31\nGaPSE.I13\nGaPSE.I11\nGaPSE.σ_0\nGaPSE.σ_1\nGaPSE.σ_2\nGaPSE.σ_3 \nGaPSE.ϕ\nGaPSE.W\nGaPSE.V_survey\nGaPSE.A\nGaPSE.A_prime","category":"page"},{"location":"Tool_functions/#GaPSE.spline_F","page":"Tool functions","title":"GaPSE.spline_F","text":" spline_F(x, μ) :: Float64\n\nReturn the 2-dim spline value of F in the given (x,μ) input. The spline is obtained through the interpolate function of the  GridInterpolations Julia package applied to the input FILE_F_MAP.\n\n\n\n\n\n","category":"function"},{"location":"Tool_functions/#GaPSE.I00","page":"Tool functions","title":"GaPSE.I00","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I20","page":"Tool functions","title":"GaPSE.I20","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I40","page":"Tool functions","title":"GaPSE.I40","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I02","page":"Tool functions","title":"GaPSE.I02","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I22","page":"Tool functions","title":"GaPSE.I22","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I31","page":"Tool functions","title":"GaPSE.I31","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I13","page":"Tool functions","title":"GaPSE.I13","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.I11","page":"Tool functions","title":"GaPSE.I11","text":" I00, I20, I40, I02, I22, I31, I13, I11 ::Float64\n\nReturn the value of the integral:\n\nI_ell^n(s)=int_0^infty fracmathrmd q2 pi^2 q^2  P(q) \n     fracj_ell(qs)(qs)^n\n\nwhere, for a generic Iab name, ell is the FIRST number (a) and  n the second (b).\n\nThese function are obtained through a Spline1D (from the  Dierckx Julia package) of the Spherical Bessel Transform function xicalc (from the  TwoFAST Julia package) applied to the  input Power Spectrum P(q).\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.σ_0","page":"Tool functions","title":"GaPSE.σ_0","text":" σ_0, σ_1, σ_2, σ_3\n\nThese are the results of the following integral:\n\n     sigma_i = int_0^infty fracmathrmd q2 pi^2 q^2-i  P(q) \n\nwhere  P(q) is the input Power Spectrum.\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.σ_1","page":"Tool functions","title":"GaPSE.σ_1","text":" σ_0, σ_1, σ_2, σ_3\n\nThese are the results of the following integral:\n\n     sigma_i = int_0^infty fracmathrmd q2 pi^2 q^2-i  P(q) \n\nwhere  P(q) is the input Power Spectrum.\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.σ_2","page":"Tool functions","title":"GaPSE.σ_2","text":" σ_0, σ_1, σ_2, σ_3\n\nThese are the results of the following integral:\n\n     sigma_i = int_0^infty fracmathrmd q2 pi^2 q^2-i  P(q) \n\nwhere  P(q) is the input Power Spectrum.\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.σ_3","page":"Tool functions","title":"GaPSE.σ_3","text":" σ_0, σ_1, σ_2, σ_3\n\nThese are the results of the following integral:\n\n     sigma_i = int_0^infty fracmathrmd q2 pi^2 q^2-i  P(q) \n\nwhere  P(q) is the input Power Spectrum.\n\n\n\n\n\n","category":"constant"},{"location":"Tool_functions/#GaPSE.ϕ","page":"Tool functions","title":"GaPSE.ϕ","text":" ϕ(s; s_min = s_min, s_max = s_max) :: Float64\n\nRadial part of the survey window function. Return 1.0 if is true that s_mathrmmin le s le s_mathrmmax and 0.0 otherwise.\n\nIn this software we made the assuption that the survey window function can be separated into a radial and angular part, i.e.:\n\n     phi(mathbfs) = phi(s)  W(hats)\n\nSee also: W\n\n\n\n\n\n","category":"function"},{"location":"Tool_functions/#GaPSE.W","page":"Tool functions","title":"GaPSE.W","text":" W(θ; θ_max = θ_MAX) :: Float64\n\nAngular part of the survey window function. Return 1.0 if is true that 00 leq theta le theta_mathrmmax and 0.0 otherwise. It is implicitly assumed an azimutal simmetry of the survey.\n\nIn this software we made the assuption that the survey window function can be separated into a radial and angular part, i.e.:\n\n     phi(mathbfs) = phi(s)  W(hats)\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"Tool_functions/#GaPSE.V_survey","page":"Tool functions","title":"GaPSE.V_survey","text":" V_survey(s_min = s_min, s_max = s_max, θ_max = θ_MAX) :: Float64\n\nReturn the volume of a survey with azimutal simmetry, i.e.:\n\nbeginsplit\n    V(s_mathrmmax s_mathrmmin theta_mathrmmax) =  C_mathrmup - C_mathrmdown + TC \n    C_mathrmup = fracpi3 s_mathrmmax^3  \n        (1 - costheta_mathrmmax)^2  (2 + costheta_mathrmmax) \n    C_mathrmdown = fracpi3 s_mathrmmin^3  \n        (1 - costheta_mathrmmax)^2  (2 + costheta_mathrmmax) \n    TC = fracpi3 (s_mathrmmax^2 + s_mathrmmin^2 + \n        s_mathrmmax s_mathrmmin)   (s_mathrmmax - s_mathrmmin) \n        costheta_mathrmmax sin^2theta_mathrmmax\nendsplit\n\n\n\n\n\n","category":"function"},{"location":"Tool_functions/#GaPSE.A","page":"Tool functions","title":"GaPSE.A","text":" A(s_min = s_min, s_max = s_max, θ_max = θ_MAX) :: Float64\n\nReturn the Power Spectrum multipole normalization coefficient A, i.e.:\n\n     A(s_mathrmmax s_mathrmmin theta_mathrmmax)= 2  pi  \n     V(s_mathrmmax s_mathrmmin theta_mathrmmax)\n\nwhere V(s_mathrmmax s_mathrmmin theta_mathrmmax) is the  survey volume.\n\nPay attention: this is NOT used for the normalization of PS, see instead A_prime\n\nSee also: V_survey\n\n\n\n\n\n","category":"function"},{"location":"Tool_functions/#GaPSE.A_prime","page":"Tool functions","title":"GaPSE.A_prime","text":" A_prime :: Float64\n\nIt's the Power Spectrum multipole normalization coefficient A^, i.e.:\n\n     A^ = frac3  A (s_mathrmmax^3 - s_mathrmmin^3) = 4 pi^2\n\nSee also: A, V_survey\n\n\n\n\n\n","category":"constant"},{"location":"F_evaluation/","page":"The F window function","title":"The F window function","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"F_evaluation/#The-F-window-function","page":"The F window function","title":"The F window function","text":"","category":"section"},{"location":"F_evaluation/","page":"The F window function","title":"The F window function","text":"GaPSE.F\nGaPSE.integrand_F\nGaPSE.F_map","category":"page"},{"location":"F_evaluation/#GaPSE.F","page":"The F window function","title":"GaPSE.F","text":" F(x, μ, θ_max; tolerance=1e-8) :: Tuple{Float64, Float64}\n\nThe function F(xmu theta_mathrmmax), defined as follows:\n\nbeginsplit\nF(xmu theta_mathrmmax) =  4pi \n    int_0^theta_mathrmmax mathrmdtheta_1 int_0^pi mathrmd theta  \n     Thetaleft(frac\n        x cos theta + cos theta_1sqrtx^1+2+2xmu - \n        cos(theta_mathrmmax) \n        right) \n     Theta(mu-cos(theta+theta_1)) \n    Theta(cos(theta - theta_1)-mu) \n    fracsinthetasintheta_1\n        sqrt(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2\nendsplit\n\ntolerance is a parameter needed in case of small negative denominator: the Heaviside theta function mathematically prevent that  mathrmden=(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2 becomes negative, but computationally might happen that mathrmden results as a very small negative number (for instance -1.2368946523-18); in this case tolerance solve the problem, returning 0 if 0-mathrmden mathrmtolerance.\n\nThe double integral is performed with hcubature function from the Julia Package HCubature; rtol, atol and all the kwargs insert into F  are directly transferred to hcubature. \n\nPAY ATTENTION: do not set too small atol and rtol, or the computation can easily become overwhelming! \n\nSee also: F_map, integrand_F\n\n\n\n\n\n","category":"function"},{"location":"F_evaluation/#GaPSE.integrand_F","page":"The F window function","title":"GaPSE.integrand_F","text":"integrand_F(θ, θ_1, x, μ, θ_max; tolerance=1e-8) :: Float64\n\nReturn the integrand of the function F(xmu theta_mathrmmax), i.e the  function f(xmu theta theta_1 theta_mathrmmax):\n\nbeginsplit\nf(xmu theta theta_1 theta_mathrmmax) = \n     Thetaleft( frac\n        x cos theta + cos theta_1sqrtx^1+2+2xmu - \n        cos(theta_mathrmmax) \n        right) \n     times  Theta(mu-cos(theta+theta_1))   times \n     quad Theta(cos(theta - theta_1)-mu)  times \n    frac4pi sinthetasintheta_1\n        sqrt(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2\nendsplit\n\nbeginequation\nF(xmu theta_mathrmmax) = int_0^theta_mathrmmax \n        mathrmdtheta_1 int_0^pi mathrmd theta \n         f(xmu theta theta_1 theta_mathrmmax)\nendequation\n\ntolerance is a parameter needed in case of small negative denominator: the Heaviside theta function mathematically prevent that  mathrmden=(sinthetasintheta_1)^2-(costhetacostheta_1-mu)^2 becomes negative, but computationally might happen that mathrmden results as a very small negative number (for instance -1.2368946523-18); in this case tolerance solve the problem, returning 0 if 0-mathrmden mathrmtolerance\n\nSee also: F, F_map\n\n\n\n\n\n","category":"function"},{"location":"F_evaluation/#GaPSE.F_map","page":"The F window function","title":"GaPSE.F_map","text":" F_map(x_step::Float64 = 0.01, μ_step::Float64 = 0.01;\n      out = \"data/F_map.txt\", x1 = 0, x2 = 3, μ1 = -1, μ2 = 1, kwargs...)\n F_map(xs::Vector{Float64}, μs::Vector{Float64};\n      out = \"data/F_map.txt\", kwargs...)\n\nEvaluate the window function F(xmu theta_mathrmmax) in a rectangual grid  of mu and x values.\n\nIn the first method you specify start, stop and step for x and μ manually, while with thr second one you pass the values (through a vector )you want to calculate  the function in.\n\nSee also: F_map, integrand_F\n\n\n\n\n\n","category":"function"},{"location":"Power_Spectrum/","page":"Power Spectrum Multipoles","title":"Power Spectrum Multipoles","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"Power_Spectrum/#The-Power-Spectrum-multipole-computation","page":"Power Spectrum Multipoles","title":"The Power Spectrum multipole computation","text":"","category":"section"},{"location":"Power_Spectrum/","page":"Power Spectrum Multipoles","title":"Power Spectrum Multipoles","text":"GaPSE.PS_multipole","category":"page"},{"location":"Power_Spectrum/#GaPSE.PS_multipole","page":"Power Spectrum Multipoles","title":"GaPSE.PS_multipole","text":" PS_multipole(in::String; int_s_min = 1e-2, int_s_max = 2 * s_max, N = 128,\n      L::Integer = 0, pr::Bool = true)\n PS_multipole(f_in::Union{Function,Dierckx.Spline1D};\n      L::Integer = 0,  N = 128, \n      int_s_min = 1e-2, int_s_max = 2 * s_max,\n      pr::Bool = true) :: Tuple{Vector{Float64}, Vector{Float64}}\n\nReturn the L-order multipole from the input function f_in, through the following Fast Fourier Transform and the effective redshift approximation:\n\nP_L(k) = frac2 L + 1A^ (-i)^L  phi(s_mathrmeff) int_0^infty \n        mathrmd s  s^2  j_L(ks)  f_mathrmin(s)  \n        quad  A^ = 4  pi^2\n\nThis expression can be easily obtained from the standard one:\n\nbeginsplit\n    P_L(k) = frac2 L + 1A (-i)^L  \n        int_0^infty mathrmd s_1  s_1^2 \n        int_0^infty mathrmd s  s^2 \n        int_-1^+1 mathrmd mu \n        j_L(ks)  xi(s_1 s mu)  phi(s_1)  phi(s_2) \n        mathcalL_L(mu) Fleft(fracss_1 mu right) \n        mathrmwith  s_2 = s_2(s_1 s μ) = sqrts_1^2 + s^2 + 2s_1smu\n         \n         quad A(s_mathrmmax s_mathrmmin theta_mathrmmax) \n        = 2  pi  V\nendsplit\n\nwith the definition\n\nf_mathrmin(s_1 s) =  int_-1^+1 mathrmd mu \n        xi(s_1 s mu)  phi(s_2) \n        mathcalL_L(mu)  Fleft(fracss_1 mu right)\n\nand the application of the effective redshift approximation.\n\nSee also: V_survey, A, A_prime\n\n\n\n\n\n","category":"function"},{"location":"","page":"Introduction","title":"Introduction","text":"DocTestSetup = quote\n    using GaPSE\nend","category":"page"},{"location":"#GaPSE.jl-:-a-Galaxy-Power-Spectrum-Estimator","page":"Introduction","title":"GaPSE.jl : a Galaxy Power Spectrum Estimator","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This is the documentation of GaPSE.jl package, an implementation of a Galaxy Power Spectrum Estimator written in Julia.","category":"page"},{"location":"#Documentation","page":"Introduction","title":"Documentation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The documentation was built using Documenter.jl.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using Dates # hide\nprintln(\"Documentation built on $(now()) using Julia $(VERSION).\") # hide","category":"page"},{"location":"#Contents","page":"Introduction","title":"Contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"#Index","page":"Introduction","title":"Index","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"}]
}
