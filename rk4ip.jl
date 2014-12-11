# Runge-Kutta 4th order in the Interaction Picture methods
RK4IP_SCAL_DEBUG = false

rk4ip_scal!(p::Pulse, f::Fiber) = 
    rk4ip_scal!(p.uX, p.t, p.w, f.L, 1e-6f.L, f.max_steps, f.adaptive_step,
                f.alpha, f.betha, f.gamma, 0., 0., 
                f.gain, f.gain_bandwidth, f.saturation_energy,
                p.fft_plan!, p.ifft_plan!, 0, 0)

@eval function rk4ip_scal!(u, t, w, L, h0, max_steps, adaptive_step,
                          alpha, betha, gamma, steep, t_raman,
                          gain, gain_bw, saturation_energy,
                          fft_plan!, ifft_plan!, nt_plot=2^8, nz_plot=2^8)
    z = 0.
    n = length(u)
    T = calc_T(t)
    dt = calc_dt(t)

    n_steps = n_steps_rejected = 0
    steps = Float64[]
    err = err_prev = 1.

    $([:($a = similar(u)) for a in 
        [:U, :_u1, :_k1, :_k2, :_k3, :_k4, :_uabs2, :_du, :_ue_cplx,
         :u_full, :u_half, :u_half2]]...)

    N! = let _uabs2 = _uabs2, _du = _du, dt = dt, gamma = gamma
            (u_, h_) -> N_rk4ip_simple_scal!(u_, h_, dt, gamma, _uabs2)
    end

    hmin = L / max_steps
    h = (adaptive_step == ADAPTIVE_STEP) ? max(h0, hmin): hmin

    d_no_gain = D_exp_no_gain(w, alpha, betha)
    g_spec = gain_spectral_factor(w, gain, gain_bw)
    d_exp = similar(u)
    dispersion_exp_scal!(d_exp, d_no_gain, g_spec, u, dt, saturation_energy)

    disp_full = exp(h/2 * d_exp)
    disp_half = exp(h/4 * d_exp)

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

    @time while z < L
        # full step
        rk4ip_step!(u, u_full, h, disp_full, N!, 
                    fft_plan!, ifft_plan!,
                    _u1, _k1, _k2, _k3, _k4)

        # 2 half-steps
        if (adaptive_step == ADAPTIVE_STEP)
            rk4ip_step!(u, u_half, h/2, disp_half, N!,
                        fft_plan!, ifft_plan!,
                        _u1, _k1, _k2, _k3, _k4)
            rk4ip_step!(u_half, u_half2, h/2, disp_half, N!,
                        fft_plan!, ifft_plan!,
                        _u1, _k1, _k2, _k3, _k4)
            
            err = integration_error_global(u_full, u_half2, _ue_cplx)
        end

        if (adaptive_step == ADAPTIVE_STEP) & (err > 1) & (h > hmin)
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
            RK4IP_VEC_DEBUG && mod(n_steps, 100) == 0 && @show (n_steps, z, h)

            if (adaptive_step == ADAPTIVE_STEP)
                err_prev = err
                BLAS.blascopy!(n, u_half2, 1, u, 1)

                h = max(h* scale_step_ok(err, err_prev), hmin)
                h = min(L - z, h)
            else
                BLAS.blascopy!(n, u_full, 1, u, 1)
            end

            dispersion_exp_scal!(d_exp, d_no_gain, g_spec, u, dt, saturation_energy)
            hd2 = h/2
            hd4 = h/4
            @devec disp_full[:] = exp(hd2 .* d_exp)
            @devec disp_half[:] = exp(hd4 .* d_exp)

            if do_plot && z >= i_plot * dz_plot
                i_plot += 1
                u_plot[i_plot, :] = u[t_plot_ind]
                spectrum!(u, U, ifft_plan!, T)
                U_plot[i_plot,:] = U[t_plot_ind]                
            end
        end
    end

    return u_plot, U_plot, n_steps, n_steps_rejected, steps
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

function N_rk4ip_simple_scal!(u, h, dt, gamma, _uabs2)
    n = length(u)
    map!(Abs2Fun(), _uabs2, u)
    k = 1im * h * gamma
    @devec u[:] = k .* u .* _uabs2
end    

function N_rk4ip_raman_scal!(u, h, dt, gamma, t_raman, _uabs2, _du)
    n = length(u)

    # _uabs2 = |u|^2 - t_raman * d(|u|^2)/dt
    map!(Abs2Fun(), _uabs2, u)
    df!(_uabs2, _du, dt)
    BLAS.axpy!(n, -t_raman + 0.im, _du, 1, _uabs2, 1)

    k = 1im * h * gamma
    @devec u[:] = k .* u .* _uabs2
end

function N_rk4ip_raman_steep_scal!(u, h, dt, gamma, t_raman, steep, _uabs2, _du)
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