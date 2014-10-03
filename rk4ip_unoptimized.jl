# Runge-Kutta 4th order in the Interaction Picture methods

function integrate_RK4IP(u, t_grid, w_grid, 
                         fft_plan, ifft_plan, fft_plan!, ifft_plan!,
                         L, h, alpha, beta, gamma, steep, t_raman,
                         nt_plot = 2^7, nz_plot = 2^7)
    z = 0.
    n_steps = n_steps_rejected = 0
    steps = Float64[]
    err_prev = 1.
    dt = t_grid[end] - t_grid[end - 1]

    nonlinear_op = let dt = dt, gamma = gamma, steep = steep, t_raman = t_raman
        (u_, h_) -> nonlinear_operator(u_, h_, dt, gamma, steep, t_raman)
    end
    d_exp = dispersion_exponent(w_grid, alpha, beta)
    disp_full = exp(h * d_exp)
    disp_half = exp(h/2. * d_exp)

    do_plot = (nt_plot != 0 && nz_plot !=0)

    dz_plot = L / (nt_plot-1)
    t_plot_ind = round(linspace(1, length(t_grid), nt_plot))
    u_plot = zeros(Complex{Float64}, nz_plot, nt_plot)
    i_plot = 1
    if do_plot u_plot[i_plot,:] = u[t_plot_ind] end

    @profile @time while z < L
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
            push!(steps, h)             
            h *= scale_step_ok(err, err_prev)
            h = min(L - z, h)
            disp_full = exp(h * d_exp)
            disp_half = exp(h/2 * d_exp)
            err_prev = err
            u = u_half

            if do_plot && z >= i_plot * dz_plot
                println("step: $n_steps, z: $z")
                i_plot += 1
                u_plot[i_plot, :] = u[t_plot_ind]
            end
        end
    end

    return (u, n_steps, n_steps_rejected, steps, u_plot)
end

function dispersion_exponent(w, alpha, beta)
    # pay attention to the order of Fourier transforms, that determine
    # the sign of differentiation operator
    -alpha/2 + 1im/2 * beta[1] * w.^2 + 1im/6 * beta[2] * w.^3 
end

function step_RK4IP(u, h, disp, nonlinear_op, fft_plan, ifft_plan)
    # TODO
    # switch to inplace transforms, which are at 1.5 times faster
    # reduce GC pressure as much as possible
    # use BLAS
    u1 = fft_plan(disp .* ifft_plan(u))
    k1 = fft_plan(disp .* ifft_plan(nonlinear_op(u, h)))
    k2 = nonlinear_op(u1 + k1/2, h)
    k3 = nonlinear_op(u1 + k2/2, h)
    k4 = nonlinear_op(fft_plan(disp .* ifft_plan(u1 + k3)), h)
    fft_plan(disp .* ifft_plan(u1 + (k1 + 2.(k2 + k3))/6)) + k4/6
end

function nonlinear_operator(u, h, dt, gamma, steep, t_raman)
    u2 = abs2(u)
    (steep == 0 && t_raman ==0) && return 1im * h * gamma * u .* u2
    steep == 0 && return 1im * h * gamma * u .* (u2 - t_raman * df(u2, dt))
    return 1im * h * gamma * u .* (u2 + (1im * steep - t_raman) * df(u2, dt) 
                                    + 1im * steep * df(u, dt) .* conj(u))
end

function df(x, dx)
    len = length(x)
    tmp = zeros(x)
    tmp[1] = (x[2] - x[end]) / (2dx)
    tmp[end] = (x[1] - x[end-1]) / (2dx)
    for i in 2:(len-1)
        @inbounds tmp[i] = (x[i+1] - x[i-1]) / (2dx)
    end
    return tmp
end

# this function has led to a lot of pain, dunno why
# function cplx_array_max{T}(a::Array{Complex{T},1}, b::Array{Complex{T},1})
#     # FIXME: maybe this is not optimal at all
#     m = [(abs2(a[i]) > abs2(b[i]) ? a[i] : b[i]) for i in 1:length(a)]
#     # println("cpl_array_max: $(typeof(m)), $(length(m)), $(m[end])")
#     # m
# end

function integration_error(u1, u2, atol=1.e-6, rtol=1.e-6)
    m = [(abs(u1[i]) > abs(u2[i]) ? abs(u1[i]) : abs(u2[i])) for i in 1:length(u1)]
    error_scale = atol + rtol * m
    maximum(abs(u1 - u2) ./ error_scale)
    # maximum(abs((abs(u1 - u2) ./ error_scale)));
end

function PI_control_factor(err, err_prev, ae=0.7, be=0.4)
    err^(-ae/5) * err_prev^(be/5)
end

function scale_step_fail(err, err_prev, ae=0.7, be=0.4)
    0.8max(1/5., PI_control_factor(err, err_prev, ae, be))
end

function scale_step_ok(err, err_prev, ae=0.7, be=0.4)
    0.8min(10, PI_control_factor(err, err_prev, ae, be))
end