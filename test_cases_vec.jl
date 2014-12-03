function soliton_vec(N=1)
    n = 2^14
    T0 = 56.7e-14
    C0 = 0.
    theta = pi / 4
    T_window = 40

    alpha = 0.
    betha2 = -1.e-26
    dbetha = 0.2 * 2 * abs(betha2) / T0
    gamma = 1e-2
    steep = 0.

    P0 = N^2 * abs(betha2) / (gamma * T0^2)
    L = 10 * T0^2 / abs(betha2)
    # L = pi / 2 * T0^2 / abs(betha2)
    run_vec(n, T_window, alpha, [betha2], dbetha, gamma, L, T0, P0, C0, theta)
end

function Agr3_6_vec(case=0, theta=0.)
    n = 2^12
    T0 = 85.0e-15
    P0 = 10.e3
    C0 = 0.  
    T_window = 20

    alpha = 0.
    b3 = 8.119e-2
    b2 = case == 0 ? 0. : b3 / (T0 / 1.e-12)
    betha = ps_km2s_m([b2, b3])
    gamma = 0.

    L = 5 * T0^3 / abs(betha[2])
    run_vec(n, T_window, alpha, betha, 0., gamma, L, T0, P0, C0, theta, 1)
end