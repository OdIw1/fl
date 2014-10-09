# Runge-Kutta 4th order in the Interaction Picture methods

@eval function rk4ip(u, L, h, t_grid, w_grid, 
                     alpha, beta, gamma, steep, t_raman,
                     fft_plan!, ifft_plan!, nt_plot=2^9, nz_plot=2^9)
    z = 0.
    n = length(u)
    T = (t_grid[end] - t_grid[1]) / 2
    dt = (t_grid[end] - t_grid[1]) / (length(t_grid)-1)

    n_steps = n_steps_rejected = 0
    steps = Float64[]
    err_prev = 1.

    ue_ = similar(u, Float64)
    $([:($a = similar(u)) for a in 
        [:U, :_u1, :_k1, :_k2, :_k3, :_k4, :_uabs2, :_du, :ue_cplx_,
         :u_full, :u_half, :u_half2]]...)

    N! = let _uabs2 = _uabs2, _du = _du, 
             dt = dt, gamma = gamma, steep = steep, t_raman = t_raman
        if t_raman == 0. && steep == 0.
            (u_, h_) -> N_simple!(u_, h_, dt, gamma, _uabs2)
        elseif steep == 0.
            (u_, h_) -> N_raman!(u_, h_, dt, gamma, t_raman, _uabs2, _du)
        else
            (u_, h_) -> N_raman_steep!(u_, h_, dt, gamma, t_raman, steep, _uabs2, _du)
        end
    end

    d_exp = dispersion_exponent(w_grid, alpha, beta)
    disp_full = exp(h/2 * d_exp)
    disp_half = exp(h/4. * d_exp)

    # prepare plotting
    do_plot = (nt_plot != 0 && nz_plot !=0)
    dz_plot = L / (nt_plot-1)
    t_plot_ind = round(linspace(1, n, nt_plot))
    u_plot = zeros(Complex{Float64}, nz_plot, nt_plot)
    U_plot = zeros(u_plot)
    i_plot = 1
    if do_plot 
        u_plot[i_plot,:] = u[t_plot_ind]
        spectrum!(u, U, ifft_plan!, T)
        U_plot[i_plot,:] = U[t_plot_ind]
    end 

    @time @profile while z < L
        # full step
        rk4ip_step!(u, u_full, h, disp_full, N!, 
                    fft_plan!, ifft_plan!,
                    _u1, _k1, _k2, _k3, _k4)
        # 2 half-steps
        rk4ip_step!(u, u_half, h/2, disp_half, N!,
                    fft_plan!, ifft_plan!,
                    _u1, _k1, _k2, _k3, _k4)
        rk4ip_step!(u_half, u_half2, h/2, disp_half, N!,
                    fft_plan!, ifft_plan!,
                    _u1, _k1, _k2, _k3, _k4)
        
        err = integration_error_global(u_full, u_half2, ue_cplx_)
        if err > 1
            n_steps_rejected += 1
            h *= scale_step_fail(err, err_prev)
            hd2 = h/2
            hd4 = h/4
            @devec disp_full[:] = exp(hd2 .* d_exp)
            @devec disp_half[:] = exp(hd4 .* d_exp)
        else
            z += h
            n_steps += 1
            push!(steps, h)
            mod(n_steps, 100) == 0 && @show (n_steps, z, h)             
            h *= scale_step_ok(err, err_prev)
            h = min(L - z, h)
            hd2 = h/2
            hd4 = h/4
            @devec disp_full[:] = exp(hd2 .* d_exp)
            @devec disp_half[:] = exp(hd4 .* d_exp)

            err_prev = err
            BLAS.blascopy!(length(u), u_half2, 1, u, 1)

            if do_plot && z >= i_plot * dz_plot
                # println("step: $n_steps, z: $z")
                i_plot += 1
                u_plot[i_plot, :] = u[t_plot_ind]
                spectrum!(u, U, ifft_plan!, T)
                U_plot[i_plot,:] = U[t_plot_ind]                
            end
        end
    end

    return (u, u_plot, U_plot, n_steps, n_steps_rejected, steps)
end

function dispersion_exponent(w, alpha, beta)
    # pay attention to the order of Fourier transforms that determine
    # the sign of differentiation operator
    # -alpha/2 + 1im/2 * beta[1] * w.^2 + 1im/6 * beta[2] * w.^3 + ...
    res = - alpha/2 + zeros(Complex{Float64}, length(w))
    for k in 1:length(beta)
        res += (1.im / factorial(k+1) * beta[k]) * w.^(k+1) # * (-1)^(k+1)
    end
    return res
end

function rk4ip_step!(u, uf, h, disp, N!, fft_plan!, ifft_plan!,
                     u1, k1, k2, k3, k4)
    n = length(u)
    BLAS.blascopy!(n, u, 1, u1, 1)
    BLAS.blascopy!(n, u, 1, k1, 1)

    # u1 = FFT(D * IFFT(u))
    ifft_plan!(u1)
    @devec u1[:] = disp .* u1           
    fft_plan!(u1)
        
    BLAS.blascopy!(n, u1, 1, k2, 1)
    BLAS.blascopy!(n, u1, 1, k3, 1)
    BLAS.blascopy!(n, u1, 1, k4, 1)
    BLAS.blascopy!(n, u1, 1, uf, 1)
    
    # k1 = FFT(D * IFFT(N(u)))
    N!(k1, h)
    ifft_plan!(k1)
    @devec k1[:] = disp .* k1
    fft_plan!(k1)
        
    # k2 = N(u1 + k1/2)
    BLAS.axpy!(n, 0.5 + 0.im, k1, 1, k2, 1)
    N!(k2, h)
    
    # k3 = N(u1 + k2/2)
    BLAS.axpy!(n, 0.5 + 0.im, k2, 1, k3, 1)
    N!(k3, h)
    
    # k4 = N(FFT(D * IFFT(u1 + k3)))
    BLAS.axpy!(n, 1. + 0.im, k3, 1, k4, 1)
    ifft_plan!(k4)
    @devec k4[:] = disp .* k4           
    fft_plan!(k4)
    N!(k4, h)
    
    # res = FFT(D * IFFW(u1 + 1/6 k1 + 1/3 (k2 + k3)) ) + 1/6 k4
    BLAS.axpy!(n, 1/6 + 0.im, k1, 1, uf, 1)
    BLAS.axpy!(n, 1/3 + 0.im, k2, 1, uf, 1)
    BLAS.axpy!(n, 1/3 + 0.im, k3, 1, uf, 1)
    ifft_plan!(uf)
    @devec uf[:] = disp .* uf
    fft_plan!(uf)           
    BLAS.axpy!(n, 1/6 + 0.im, k4, 1, uf, 1)
end

function N_simple!(u, h, dt, gamma, _uabs2)
    n = length(u)
    map!(Abs2Fun(), _uabs2, u)
    k = 1im * h * gamma
    @devec u[:] = k .* u .* _uabs2
end    

function N_raman!(u, h, dt, gamma, t_raman, _uabs2, _du)
    n = length(u)

    # _uabs2 = |u|^2 - t_raman * d(|u|^2)/dt
    map!(Abs2Fun(), _uabs2, u)
    df!(_uabs2, _du, dt)
    BLAS.axpy!(n, -t_raman + 0.im, _du, 1, _uabs2, 1)

    k = 1im * h * gamma
    @devec u[:] = k .* u .* _uabs2
end

function N_raman_steep!(u, h, dt, gamma, t_raman, steep, _uabs2, _du)
    # steep = 1.im / w0
    n = length(u)

    # _uabs2 = |u|^2 + (steep - t_raman) * 2Re[conj(u) * d(u)/dt] + steep[conj(u) * d(u)/dt]
    map!(Abs2Fun(), _uabs2, u)
    # _du = d(u)/dt * conj(u)
    df!(u, _du, dt)
    @simd for i = 1:n
        @inbounds _du[i] = _du[i] * conj(u[i])
    end
    
    BLAS.axpy!(n, steep, _du, 1, _uabs2, 1)
    @simd for i = 1:n
        @inbounds _uabs2[i] += (steep - t_raman) * 2.real(_du[i])
    end

    k = 1im * h * gamma
    @devec u[:] = k .* u .* _uabs2
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