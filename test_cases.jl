function Agr5_23()
    n = 2^14
    T0 = 28.0e-15
    P0 = 3.01e8
    C0 = 0.
    T_window = 50

    alpha = 0.
    b2 = -11.830
    b3 = 8.1038e-2
    b4 = -9.5205e-5
    b5 = 2.0737e-7
    b6 = -5.3943e-10
    b7 = 1.3486e-12
    b8 = -2.5495e-15
    b9 = 3.0524e-18
    b10 = -1.7140e-21
    beta_pskm = [b2, b3, b4, b5, b6, b7, b8, b9, b10]
    beta = ps_km2s_m(beta_pskm)
    # beta = (-5.92e-20, 2.98e-34)
    gamma = 1.
    wl = 2.64e-6
    steep = 1.im / (2pi * 3e8 / wl)
    t_raman = 2.80e-15
    L = 5.3e-8
    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
end

function soliton(N=1)
    n = 2^14
    T0 = 56.7e-15
    C0 = 0.
    T_window = 20

    alpha = 0.
    beta2 = -1.e-26
    gamma = 1e-2
    steep = 0.
    t_raman = 0. # 2.80e-15

    P0 = N^2 * abs(beta2) / (gamma * T0^2)
    L = pi / 2 * T0^2 / abs(beta2)
    run_scalar(n, T_window, alpha, [beta2], gamma, t_raman, steep, L, T0, P0, C0)
end

function Heidt()
    n = 2^14
    T0 = 28.4e-15
    P0 = 10.e3
    C0 = 0. 
    T_window = 8.e-12/T0

    alpha = 0.
    b2 = -11.830
    b3 = 8.1038e-2
    b4 = -9.5205e-5
    b5 = 2.0737e-7
    b6 = -5.3943e-10
    b7 = 1.3486e-12
    b8 = -2.5495e-15
    b9 = 3.0524e-18
    b10 = -1.7140e-21
    beta = ps_km2s_m([b2, b3, b4, b5, b6, b7, b8, b9, b10])
    gamma = 0.11
    steep = 1.im * 0.56e-15    
    t_raman = 2.80e-15
    
    L = 0.15
    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
end

function Agr13_13()
    n = 2^14
    T0 = 85.0e-15
    P0 = 10.e3
    C0 = 0.
    T_window = 5.e-12 / T0

    alpha = 0.
    b2 = -12.76
    b3 = 8.119e-2
    b4 = -1.321e-4
    b5 = 3.032e-7
    b6 = -4.196e-10
    b7 = 2.57e-13
    beta = ps_km2s_m([b2, b3, b4, b5, b6, b7])
    gamma = 0.10
    steep = 1.im * 0.56e-15    
    t_raman = 2.80e-15

    L = 0.10
    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
end

function Agr3_6(case=0)
    n = 2^12
    T0 = 85.0e-15
    P0 = 10.e3
    C0 = 0.  
    T_window = 20

    alpha = 0.
    b3 = 8.119e-2
    b2 = case == 0 ? 0. : b3 / (T0 / 1.e-12)
    betha = ps_km2s_m([b2, b3])
    gamma = 0.00
    steep = 0 # 1.im * 0.56e-15    
    t_raman = 0 # 2.80e-15

    L = 5 * T0^3 / abs(betha[2])
    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0, 1)
end

function Agr3_11(case=0)
    n = 2^12
    T0 = 0.5e-12 / sqrt(4*log(2)) # T0 != FWHM, hence the difference from the text
    P0 = 10.e3
    C0 = 0.  
    T_window = 20

    alpha = 0.
    b2 = 0.
    b3 = case == 0 ? 0.124 : -0.076
    betha = ps_km2s_m([b2, b3])
    gamma = 0.00
    steep = 0 # 1.im * 0.56e-15    
    t_raman = 0 # 2.80e-15

    L = 2.5e3
    run_scalar(n, T_window, alpha, betha, gamma, t_raman, steep, L, T0, P0, C0, 1)
end

function Agr4_15(case=0)
    n = 2^12
    T0 = 85.0e-15
    P0 = 10.e3
    C0 = 0.  
    T_window = 20

    alpha = 0.
    b3 = 8.119e-2
    b2 = case == 0 ? 0. : b3 / (T0 / 1.e-12)
    beta = ps_km2s_m([b2, b3])
    gamma = abs(betha[2]) / (P0 * T0^3)
    steep = 0 # 1.im * 0.56e-15    
    t_raman = 0 # 2.80e-15

    L = 5 * T0^3 / abs(beta[2])
    run_scalar(n, T_window, alpha, betha, gamma, t_raman, steep, L, T0, P0, C0, 1)
end

function Agr4_16(N=10)
    n = 2^12
    T0 = 85.0e-15
    P0 = 10.e3
    C0 = 0.  
    T_window = 15

    alpha = 0.
    b3 = 8.119e-2
    b2 = 0
    beta = ps_km2s_m([b2, b3])
    gamma = N^2 * abs(beta[2]) / (P0 * T0^3)
    steep = 0 # 1.im * 0.56e-15    
    t_raman = 0 # 2.80e-15

    L = 0.2 * T0^3 / abs(beta[2])
    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0, 1)
end

function Agr4_19(k=20)
    n = 2^12
    T0 = 1.
    P0 = 1.
    C0 = 0.  
    T_window = 15

    alpha = 0.
    beta = [0., 0.]
    gamma = 1.
    steep = 0.01im # 1.im * 0.56e-15    
    t_raman = 0 # 2.80e-15

    L = k / (gamma * P0)
    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0, 1)
end

function Agr4_20(case=1.)
    n = 2^12
    T0 = 1.
    P0 = 1.
    C0 = 0.  
    T_window = 50

    alpha = 0.
    gamma = 1.
    sign = case == 1 ? 1 : -1
    beta = [sign * 1. / 2^2, 0.]
    
    steep = 1.e-30im # 1.im * 0.56e-15    
    t_raman = 0.03 # 2.80e-15
    
    L = 8 * T0^2/ abs(beta[1])

    run_scalar(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0, 1)
end