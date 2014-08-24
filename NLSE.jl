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
    sqrt(P0) * exp(-0.5 * (1 + 1im* C0) * t_scaled .^ (2 * m_order))
end

function integrate_RK4IP(alpha, beta, gamma, steep, t_raman, L, h0, u0)

end