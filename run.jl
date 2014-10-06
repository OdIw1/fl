function x()

    n = 2^14
    alpha = 0.
    beta = (-5.92e-20, 2.98e-34)
    gamma = 1.
    T0 = 28.0e-15
    P0 = 3.01e8
    C0 = 0.
    wl = 2.64e-6
    steep = 0.
    t_raman = 2.80e-15
    #t_raman = 0.
    T_window = 50
    L = 5.3e-8
    T = T_window * T0

    @show L_N, L_D = lnd(T0, P0, gamma, beta)

    FFTW.set_num_threads(1)
    BLAS.blas_set_num_threads(1) # this function is erroneously reported in REPL

    t = t_grid(n, T)
    w = w_grid(n, T)
    u0 = secant_pulse(T0, P0, C0, t)
    E0 = sum(abs2(u0)) * (t[end] - t[end-1])

    u = copy(u0)
    uf = similar(u0)
    U0 = similar(u0)
    Uf = similar(u0)
    fft_plan! = plan_fft!(uf, (1,), FFTW.MEASURE)
    ifft_plan! = plan_ifft!(uf, (1,), FFTW.MEASURE)

    spectrum(u0, U0, ifft_plan!, T)
    fwrite("u0.tsv", u0)
    fwrite("U0.tsv", U0)

    (u, n_steps, n_steps_rejected, steps, u_plot) = 
        rk4ip(u, L, 1.e-4L, t, w, alpha, beta, gamma, steep, t_raman,
              fft_plan!, ifft_plan!)
    fwrite("u_plot.tsv", u_plot)
    fwrite("steps.tsv", steps)

    # ssfm(u, L, L/2^14, t, w, alpha, beta, gamma, steep, t_raman, fft_plan!, ifft_plan!)

    spectrum(u, Uf, ifft_plan!, T)
    fwrite("uf.tsv", u)
    fwrite("Uf.tsv", Uf)

    Ef = sum(abs2(u)) * (t[end] - t[end-1])
    @show (E0, Ef)
    @show (n_steps, n_steps_rejected)
end

function fwrite(fname, data)
    touch(fname)
    f = open(fname, "w+")
    writedlm(f, data)
    close(f)
end