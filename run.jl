function run(n, T_window, alpha, beta, gamma, t_raman, steep, L, T0, P0, C0, shape=0)
    ln, ld, soliton_order = pulse_propagation_params(T0, P0, gamma, beta...)
    @show ln, ld
    @show soliton_order

    FFTW.set_num_threads(1)
    BLAS.blas_set_num_threads(1) # this function is erroneously reported in REPL

    T = T_window * T0
    t = t_grid(n, T)
    w = w_grid(n, T)
    u0 = shape == 0. ? secant_pulse(T0, P0, C0, t) : gaussian_pulse(shape, T0, P0, C0, t)
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
    @show (E0, Ef, E0/Ef - 1)
    @show (n_steps, n_steps_rejected)
end