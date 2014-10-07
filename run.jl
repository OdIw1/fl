function run5_23()
    n = 2^14
    T_window = 50
    alpha = 0.
    # b2 = -11.830
    # b3 = 8.1038e-2
    # b4 = -9.5205e-5
    # b5 = 2.0737e-7
    # b6 = -5.3943e-10
    # b7 = 1.3486e-12
    # b8 = -2.5495e-15
    # b9 = 3.0524e-18
    # b10 = -1.7140e-21
    # beta_pskm = (b2, b3, b4, b5, b6, b7, b8, b9, b10)
    # beta = beta_pskm_to_sm(beta_pskm...)
    beta = (-5.92e-20, 2.98e-34)
    gamma = 1.
    T0 = 28.0e-15
    P0 = 3.01e8
    C0 = 0.
    wl = 2.64e-6
    steep = 1.im / (2pi * 3e8 / wl)
    t_raman = 2.80e-15
    L = 5.3e-8
    run(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
end

function run_soliton2nd()
    n = 2^14
    T_window = 100
    alpha = 0.
    beta = 1.e-26
    gamma = 1e-2
    T0 = 56.7e-15
    P0 = 1.24e3
    C0 = 0.
    steep = 0.
    t_raman = 0. # 2.80e-15
    L = 10.
    run(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
end

function run_Heidt()
    n = 2^13
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
    beta = beta_pskm_to_sm(b2, b3)
    gamma = 0.11

    T0 = 34.0e-15
    P0 = 1.e10
    C0 = 0.
    
    steep = 1.im * 0.56e-15    
    t_raman = 2.80e-15
    L = 15e-2
    run(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
end

function run(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0)
    ln, ld, soliton_order = pulse_propagation_params(T0, P0, gamma, beta...)
    @show ln, ld
    @show soliton_order

    FFTW.set_num_threads(1)
    BLAS.blas_set_num_threads(1) # this function is erroneously reported in REPL

    T = T_window * T0
    t = t_grid(n, T)
    w = w_grid(n, T)
    u0 = secant_pulse(T0, P0, C0, t)
    E0 = sum(abs2(u0)) * (t[end] - t[end-1])

    u = copy(u0)
    u1 = similar(u0)
    U0 = similar(u0)
    U1 = similar(u0)
    fft_plan! = plan_fft!(u1, (1,), FFTW.MEASURE)
    ifft_plan! = plan_ifft!(u1, (1,), FFTW.MEASURE)

    spectrum!(u0, U0, ifft_plan!, T)
    
    # add directory creation    
    outdir = "out"
    fwrite(joinpath(outdir, "u0.tsv"), abs2(u0))
    fwrite(joinpath(outdir, "U0.tsv"), abs2(U0))

    (u1, u_plot, U_plot, n_steps, n_steps_rejected, steps) = 
        rk4ip(u, L, 1.e-10L, t, w, alpha, beta, gamma, steep, t_raman,
              fft_plan!, ifft_plan!)

    fwrite(joinpath(outdir, "u_log_plot.tsv"), clamp_log_plot(u_plot))
    fwrite(joinpath(outdir, "U_log_plot.tsv"), clamp_log_plot(U_plot))
    fwrite(joinpath(outdir, "u_plot.tsv"), clamp_plot(u_plot))
    fwrite(joinpath(outdir, "U_plot.tsv"), clamp_plot(U_plot))  

    spectrum!(u1, U1, ifft_plan!, T)
    fwrite(joinpath(outdir, "u1.tsv"), abs2(u1))
    fwrite(joinpath(outdir, "U1.tsv"), abs2(U1))

    fwrite(joinpath(outdir, "steps.tsv"), steps)

    Ef = sum(abs2(u)) * (t[end] - t[end-1])
    @show (E0, Ef)
    @show (n_steps, n_steps_rejected)
end