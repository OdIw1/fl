propagate_through!(p::Pulse, M::JonesMatrix) = apply_Jones_matrix!(M, p)

function propagate_through!(p::Pulse, fiber::Fiber)
    u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
        rk4ip_vec!(p, fiber)    
end

function propagate_through!(p::Pulse, fout::FileOutput)
    UX = similar(p.uX)
    UY = similar(p.uY)
    T = (p.t[end] - p.t[1]) / 2   
    spectrum!(p.uX, UX, p.ifft_plan!, T)
    spectrum!(p.uY, UY, p.ifft_plan!, T)

    i = fout.iteration
    postfix = fout.postfix
    fwrite(outdir, "uX" * postfix, i, p.uX)
    fwrite(outdir, "uY" * postfix, i, p.uY)
    fwrite(outdir, "UX" * postfix, i, UX)
    fwrite(outdir, "UY" * postfix, i, UY)

    fout.iteration += 1
end

function run_laser!(p::Pulse, laser::LaserScheme, n_iter=2)
    for i = 1:n_iter, j = 1:length(laser)
        propagate_through!(p, laser[j])
    end
end

function run_laser(n, T_window, scheme, L, T0, P0, C0, theta, shape=0)
    P1 = W4(a1)
    P2 = W4(a2)
    P3 = W2(a3)
    P23p = P2 * P3 * [1 0; 0 0]

    uX, uY = pulse_vec(shape, T0, P0, C0, theta, t)

    outdir = mkpath_today()
    
    n_iter = 10

    for i in 1:n_iter
        apply_Jones_matrix!(P1, uX, uY)

        u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
            rk4ip_vec!(uX, uY, L, 1.e-10L, t, w, alpha, betha, dbetha, gamma,
                       gain, gain_bandwidth, E_sat, fft_plan!, ifft_plan!, 0, 0)

        apply_Jones_matrix!(P23p, uX, uY)

        spectrum!(uX, UX, ifft_plan!, T)
        spectrum!(uY, UY, ifft_plan!, T)

        fwrite(outdir, "uX", i, uX)
        fwrite(outdir, "uY", i, uY)
        fwrite(outdir, "UX", i, UX)
        fwrite(outdir, "UY", i, UY)       
    end
end

