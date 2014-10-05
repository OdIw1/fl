function lnd(T0, P0, gamma, beta::Number)
# characteristic nonlinear and dispersive lenghts
# beta can be a numer or an iterable
    lnd(T0, P0, gamma, [beta])
end

function lnd(T0, P0, gamma, beta)
    L_N = 1. / (gamma * P0)
    L_D = Float64[]
    for (i, b) in enumerate(beta)
        push!(L_D, T0^(i+1) / abs(b))
    end
    L_N, L_D
end

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

function secant_pulse(T0, P0, C0, t::Vector, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * sech(t_scaled) .* exp(-0.5im * C0 * (t_scaled).^2)
end

function gaussian_pulse(m_order, T0, P0, C0, t::Vector, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * exp(-0.5(1 + 1im * C0) * t_scaled .^ (2m_order))
end

function spectrum(u, U, ifft_plan!, T)
    n = length(u)
    BLAS.blascopy!(n, u, 1, U, 1)
    ifft_plan!(U)
    fftshift!(U)
    BLAS.scal!(n, T * sqrt(2/pi), U, 1)
end

function fftshift!(u)
    n = length(u)
    if mod(n, 2) != 0 
        throw(ArgumentError("fftshift! is defined only for even-length arrays"))
    end

    shift = div(n, 2)
    for i = 1:shift
        u[i], u[shift+i] = u[shift+i], u[i]
    end
end

