
    P1 = W4(a1)
    P2 = W4(a2)
    P3 = W2(a3)
    P23 = P2 * P3

    u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
        rk4ip_vec!(uX, uY, L, 1.e-10L, t, w, alpha, betha, dbetha, gamma,
                   gain, gain_bandwidth, E_sat, fft_plan!, ifft_plan!, 0, 0)

    spectrum!(uX, UX, ifft_plan!, T);                spectrum!(uY, UY, ifft_plan!, T)

