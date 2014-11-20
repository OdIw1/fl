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