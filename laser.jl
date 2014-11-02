immutable type Fiber{T<:Real}
    alpha::T
    betha::Vector{T} # maybe add separate bethas for both polarizations
    dbetha::T
    gamma::T
    gain::T
    gain_bandwidth::T
    saturation_energy::T
end

# no gain case
Fiber(alpha, betha, dbetha, gamma) = Fiber(alpha, betha, dbetha, gamma, 0., 1.e40, 1.e40)

# no birefrigence case
Fiber(alpha, betha, gamma) = Fiber(alpha, betha, 0., gamma)

function laser(n, T_window, alpha, betha, dbetha, gamma, L, T0, P0, C0, theta, shape=0)
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

