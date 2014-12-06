function run_scalar(n, T_window, alpha, betha, gamma, t_raman, steep, L, T0, P0, C0, shape=0)
    ld, ln, soliton_order = pulse_propagation_params(T0, P0, betha, gamma)
    @show L
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
    v0 = similar(u0)
    v1 = similar(u0)
    fft_plan! = plan_fft!(u1, (1,), FFTW.MEASURE)
    ifft_plan! = plan_ifft!(u1, (1,), FFTW.MEASURE)

    spectrum!(u0, v0, ifft_plan!, T)
    
    # add directory creation    
    outdir = "/mnt/hgfs/VM_shared/out/"
    fwrite(joinpath(outdir, "u0.tsv"), abs2(u0))
    fwrite(joinpath(outdir, "v0.tsv"), abs2(v0))

    (u1, u_plot, v_plot, n_steps, n_steps_rejected, steps) = 
        rk4ip(u, L, 1.e-10L, t, w, alpha, betha, gamma, steep, t_raman,
              fft_plan!, ifft_plan!)

    fwrite(joinpath(outdir, "u_log_plot.tsv"), clamp_log_plot(u_plot))
    fwrite(joinpath(outdir, "v_log_plot.tsv"), clamp_log_plot(v_plot))
    fwrite(joinpath(outdir, "u_plot.tsv"), clamp_plot(u_plot))
    fwrite(joinpath(outdir, "v_plot.tsv"), clamp_plot(v_plot))  

    spectrum!(u1, v1, ifft_plan!, T)
    fwrite(joinpath(outdir, "u1.tsv"), abs2(u1))
    fwrite(joinpath(outdir, "v1.tsv"), abs2(v1))

    fwrite(joinpath(outdir, "steps.tsv"), steps)

    Ef = sum(abs2(u)) * (t[end] - t[end-1])
    @show (E0, Ef, E0/Ef - 1)
    @show (n_steps, n_steps_rejected)
end


function run_vec(n, T_window, alpha, betha, dbetha, gamma, L, T0, P0, C0, theta, shape=0)
    ld, ln, soliton_order = pulse_propagation_params(T0, P0, betha, gamma)
    @show L
    @show ln, ld
    @show soliton_order

    FFTW.set_num_threads(1)
    BLAS.blas_set_num_threads(1) # this function is erroneously reported in REPL

    T = T_window * T0
    t = t_grid(n, T)
    w = w_grid(n, T)

    u0X, u0Y = pulse_vec(shape, T0, P0, C0, theta, t)

    gain = 0.
    gain_bandwidth = 1.e40
    E_sat = 1.e40

    u1X = copy(u0X);                                 u1Y = copy(u0Y)
    v0X = copy(u0X);                                 v0Y = copy(u0Y)
    v1X = copy(u0X);                                 v1Y = copy(u0Y)
    fft_plan! = plan_fft!(copy(u0X), (1,), FFTW.MEASURE)
    ifft_plan! = plan_ifft!(copy(u0X), (1,), FFTW.MEASURE)
    

    u_plotX, u_plotY, v_plotX, v_plotY, n_steps, n_steps_rejected, steps = 
        rk4ip_vec!(u1X, u1Y, L, 1.e-10L, t, w, alpha, betha, dbetha, gamma,
                   gain, gain_bandwidth, E_sat, fft_plan!, ifft_plan!)

    # add directory creation    
    outdir = "/mnt/hgfs/VM_shared/out/"

    spectrum!(u0X, v0X, ifft_plan!, T);             spectrum!(u0Y, v0Y, ifft_plan!, T)
    spectrum!(u1X, v1X, ifft_plan!, T);             spectrum!(u1Y, v1Y, ifft_plan!, T)
    
    @outfv outdir abs2 u0 v0 u1 v1
    # @outfv outdir clamp_plot u_plot U_plot
    @outfv outdir clamp_plot u_plot v_plot
    @out outdir steps

    E0 = (sum(abs2(u0X)) + sum(abs2(u0Y))) * (t[end] - t[end-1])
    Ef = (sum(abs2(u1X)) + sum(abs2(u1Y))) * (t[end] - t[end-1])
    @show (E0, Ef, E0/Ef - 1)
    @show (n_steps, n_steps_rejected)
end