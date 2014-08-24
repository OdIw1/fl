module NLSE
export create_t_grid, create_w_grid

function t_grid(n, T)
# In MATLAB version 'T' was a scaling factor: time window was (-T0 T, T0 T),
# where T0 is pulse FWHM;  in Julia time window is (-T, T)
# also 'n' is no longer binary log of the number of grid points, it is the number itself
    dt = 2 * T / n
    dt * [-n/2 : (n/2 - 1)]
end

function w_grid(n, T)
    dw = pi / T
    dw * [0 : (n/2 -1), -n/2 : -1]
end

function secant_pulse(T0, P0, C0, t, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * sech(t_scaled) .* exp(-0.5im * C0 * (t_scaled).^2)
end

function gaussian_pulse(m_order, T0, P0, C0, t, t_offset)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * exp(-0.5(1 + 1im * C0) * t_scaled .^ (2m_order))
end

function integrate_RK4IP(u0, L, h0, alpha, beta, gamma, steep, t_raman)
    z = 0
    n_steps = 0
    n_steps_rejected = 0
    steps = Float64[]
    error_prev = 1
    # PI control terms
    ae = 0.7
    be = 0.4

    while z < L
    end
end

function dispersion_exponent(w, apha, beta :: (Float64, Float64))
    beta2 = beta[1]
    beta3 = beta[2]
    # FIXME: check all beta signs, coz i'm pretty sure some are wrong!!!
    -alpha / 2 + 0.5im * beta2 * w.^2 - 1im * beta3 / 6 * w.^3 
end

function step_RK4IP(u, h, disp, nonlinear_op, fft_plan, ifft_plan)
    # FIXME: maybe later vectorize all this stuff
    u1 = fft_plan(disp .* ifft_plan(u))
    k1 = fft_plan(disp .* ifft_plan(nonlinear_op(u, h)))
    k2 = nonlinear_op(u1 + k1 / 2, h)
    k3 = nonlinear_op(u1 + k2 / 2, h)
    k4 = nonlinear_op(fft_plan(disp .* ifft_plan(u1 + k3)), h)
    fft_plan(disp .* ifft_plan(u1 + (k1 + 2.(k2 +k3)) / 6)) + k4 / 6
end

function cplx_array_max(a, b)
    len = length(a)
    u = zeros(len)
    for i in 1:len
        u[i] = abs(a[i]) > abs(b[i]) ? a[i] : b[i]
    end
    return u
end

function integration_error(u1, u2, atol, rtol)
    error_scale = atol + rtol * cplx_array_max(u1, u2)
    maximum(abs((u1 -u2) ./ error_scale))
end



