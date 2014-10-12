function soliton_vec(N=1)
    n = 2^14
    T0 = 56.7e-15
    C0 = 0.
    T_window = 20

    alpha = 0.
    beta2 = -1.e-26
    gamma = 1e-2
    steep = 0.

    P0 = N^2 * abs(beta2) / (gamma * T0^2)
    L = pi / 2 * T0^2 / abs(beta2)
    run_vec(n, T_window, alpha, beta2, gamma, L, T0, P0, C0)
end