module NLSE
export t_grid, w_grid, secant_pulse, gaussian_pulse, integrate_RK4IP

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

function gaussian_pulse(m_order, T0, P0, C0, t::Vector, t_offset)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * exp(-0.5(1 + 1im * C0) * t_scaled .^ (2m_order))
end

function integrate_RK4IP{T}(u0::Vector{Complex{T}},
                            t_grid::Vector{T},
                            w_grid::Vector{T},
                            fft_plan::Function,
                            ifft_plan::Function,
                            L, h0, alpha, beta, gamma, steep, t_raman)
    u = u0
    h = h0
    z = 0
    n_steps = 0
    n_steps_rejected = 0
    steps = (Float64)[]
    err_prev = 1

    dt = t_grid[end] - t_grid[end - 1]
    nonlinear_op = (u_, h_) -> nonlinear_operator(u_, h_, dt, gamma, steep, t_raman)

    tic()

    d_exp = dispersion_exponent(w_grid, alpha, beta)
    disp_full = exp(h * d_exp)
    disp_half = exp(h/2 * d_exp)

    while z < L
        # println("Step $(length(steps)+1), z = $z, h = $h")
        # full step
        u_full = step_RK4IP(u, h, disp_full, nonlinear_op, fft_plan, ifft_plan)
        # 2 half-steps
        u_half = step_RK4IP(u, h/2, disp_half, nonlinear_op, fft_plan, ifft_plan)
        u_half = step_RK4IP(u_half, h/2, disp_half, nonlinear_op, fft_plan, ifft_plan)

        err = integration_error(u_full, u_half)
        if err > 1
            n_steps_rejected += 1
            h *= scale_step_fail(err, err_prev)
            disp_full = exp(h * d_exp)
            disp_half = exp(h/2 * d_exp)
        else
            z += h
            n_steps += 1
            append!(steps, [h])             
            h *= scale_step_ok(err, err_prev)
            h = min(L - z, h)
            disp_full = exp(h * d_exp)
            disp_half = exp(h/2 * d_exp)
            err_prev = err
            u = u_half
        end
    end

    println("Time: $(toq())")

    return (u, n_steps, n_steps_rejected, steps)
end

function dispersion_exponent(w, alpha, beta :: (Float64, Float64))
    beta2 = beta[1]
    beta3 = beta[2]
    # FIXME: check all beta signs, coz i'm pretty sure some are wrong!!!
    -alpha / 2 + 0.5im * beta2 * w.^2 - 1im * beta3 / 6 * w.^3 
end

function step_RK4IP(u, h, disp, nonlinear_op::Function, fft_plan, ifft_plan)
    # FIXME: maybe later vectorize all this stuff
    u1 = fft_plan(disp .* ifft_plan(u))
    k1 = fft_plan(disp .* ifft_plan(nonlinear_op(u, h)))
    k2 = nonlinear_op(u1 + k1 / 2, h)
    k3 = nonlinear_op(u1 + k2 / 2, h)
    k4 = nonlinear_op(fft_plan(disp .* ifft_plan(u1 + k3)), h)
    fft_plan(disp .* ifft_plan(u1 + (k1 + 2.(k2 +k3)) / 6)) + k4 / 6
end

function nonlinear_operator(u, h, dt, gamma, steep, t_raman)
    u2 = abs2(u)
    if steep == 0 && t_raman ==0
        return 1im * h * gamma * u .* u2
    elseif steep == 0
        return 1im * h * gamma * u .* (u2 - t_raman * df(u2, dt))
    else
        return 1im * h * gamma * u .* (u2 + (1im * steep - t_raman) * df(u2, dt) 
            + 1im * steep * df(u, dt) .* conj(u))
    end
end

function df(x, dx)
    len = length(x)
    tmp = zeros(len)
    tmp[1] = (x[2] - x[end]) / (2dx)
    tmp[end] = (x[1] - x[end-1]) / (2dx)
    for i in 2:(len-1)
        tmp[i] = (x[i+1] - x[i-1]) / (2dx)
    end
    return tmp
end

function cplx_array_max{T <: Real}(a::Vector{Complex{T}}, b::Vector{Complex{T}})
    return [abs(a[i]) > abs(b[i]) ? a[i] : b[i] for i in 1:length(a)]  
end

function integration_error(u1::Vector, u2::Vector, atol=1.e-8, rtol=1.e-8)
    error_scale = atol + rtol * cplx_array_max(u1, u2)
    maximum(abs((u1 - u2) ./ error_scale))
end

function PI_control_factor(err, err_prev, ae=0.7, be=0.4)
    err^(-ae/5) * err_prev^(be/5)
end

function scale_step_fail(err, err_prev, ae=0.7, be=0.4)
    0.8max(1/5, PI_control_factor(err, err_prev, ae, be))
end

function scale_step_ok(err, err_prev, ae=0.7, be=0.4)
    0.8min(10, PI_control_factor(err, err_prev, ae, be))
end

end


