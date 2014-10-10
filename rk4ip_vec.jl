# Runge-Kutta 4th order in the Interaction Picture methods

@eval function rk4ip_vec(uX, uY, L, h, t_grid, w_grid, 
                     alpha, beta, gamma, steep, t_raman,
                     fft_plan!, ifft_plan!, nt_plot=2^9, nz_plot=2^9)
    z = 0.
    n = length(uX)
    T = (t_grid[end] - t_grid[1]) / 2
    dt = (t_grid[end] - t_grid[1]) / (length(t_grid)-1)

    n_steps = n_steps_rejected = 0
    steps = Float64[]
    err_prev = 1.

    _ueX = similar(uX, Float64);                     _ueY = similar(uX, Float64)               
    $([:($a = similar(uX)) for a in 
        [:UX, :_u1X, :_k1X, :_k2X, :_k3X, :_k4X, :_uabs2X, :_duX, :_ue_cplxX, :u_fullX, :u_halfX, :u_half2X,
         :UY, :_u1Y, :_k1Y, :_k2Y, :_k3Y, :_k4Y, :_uabs2Y, :_duY, :_ue_cplxY, :u_fullY, :u_halfY, :u_half2Y]]...)

    N! = let _uabs2A = _uabs2X, _uabs2B = _uabs2Y 
             dt = dt, gamma = gamma, steep = steep, t_raman = t_raman
        (uA_, uB_, h_) -> N_simple!(uA_, uB_, h_, dt, gamma, _uabs2A, _uabs2B)
    end

    d_exp = dispersion_exponent(w_grid, alpha, gain_bandwith, beta)
    disp_full = exp(h/2 * d_exp)
    disp_half = exp(h/4. * d_exp)

    # prepare plotting
    do_plot = (nt_plot != 0 && nz_plot !=0)
    dz_plot = L / (nt_plot-1)
    t_plot_ind = round(linspace(1, n, nt_plot))
    $([:($a = zeros(Complex{Float64}, nz_plot, nt_plot)) for a in [:u_plotX, :u_plotY, :U_plotX, :U_plotY]]...)
    
    handle_plot_data! = let t_plot_ind = t_plot_ind, uX = uX, uY = uY, UX = UX, UY = UY,
                            u_plotX = u_plotX, u_plotY = u_plotY, U_plotX = U_plotX, U_plotY = U_plotY
        (i_plot_) -> handle_plot_data!(i_plot_, t_plot_ind, uX, uY, UX, UY, u_plotX, u_plotY, U_plotX, U_plotY)
    end
    
    i_plot = 1
    do_plot && handle_plot_data!(i_plot)

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
                handle_plot_data!(i_plot)            
            end
        end
    end

    return (u, u_plot, U_plot, n_steps, n_steps_rejected, steps)
end

function handle_plot_data!(i_plot, t_plot_ind, uX, uY, UX, UY,
                           u_plotX, u_plotY, U_plotX, U_plotY)
    u_plotX[i_plot,:] = uX[t_plot_ind];                 u_plotY[i_plot,:] = uY[t_plot_ind]             
    spectrum!(uX, UX, ifft_plan!, T);                   spectrum!(uY, UY, ifft_plan!, T)
    U_plotX[i_plot,:] = UX[t_plot_ind];                 U_plotY[i_plot,:] = UY[t_plot_ind]
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

function gain_spectral_factor(w, gain_bandwidth)
    1. / (1 + w.^2 / gain_bandwidth^2)
end

function gain_saturated(uX, uY, dt, gain, saturation_energy)
    n = length(uX)
    EX = sqr(BLAS.nrm2(n, uX, 1)) * dt;             EY = sqr(BLAS.nrm2(n, uY, 1)) * dt
    energy = EX + EY
    gain_saturated = 0.5gain / (1. + energy / saturation_energy)
end

function rk4ip_step!(uX, ufX, uY, ufY, h, disp, N!, fft_plan!, ifft_plan!,
                     u1X, k1X, k2X, k3X, k4X, u1Y, k1Y, k2Y, k3Y, k4Y)
    n = length(u)
    BLAS.blascopy!(n, uX, 1, u1X, 1);               BLAS.blascopy!(n, uY, 1, u1Y, 1)
    BLAS.blascopy!(n, uX, 1, k1X, 1);               BLAS.blascopy!(n, uY, 1, k1Y, 1)

                            # u1 = FFT(D * IFFT(u))
    ifft_plan!(u1X);                                ifft_plan!(u1Y)
    @devec u1X[:] = disp .* u1X;                    @devec u1Y[:] = disp .* u1Y       
    fft_plan!(u1X);                                 fft_plan!(u1Y)   
        
    BLAS.blascopy!(n, u1X, 1, k2X, 1);              BLAS.blascopy!(n, u1Y, 1, k2Y, 1)
    BLAS.blascopy!(n, u1X, 1, k3X, 1);              BLAS.blascopy!(n, u1Y, 1, k3Y, 1)
    BLAS.blascopy!(n, u1X, 1, k4X, 1);              BLAS.blascopy!(n, u1Y, 1, k4Y, 1)
    BLAS.blascopy!(n, u1X, 1, ufX, 1);              BLAS.blascopy!(n, u1Y, 1, ufY, 1)
    
                            # k1 = FFT(D * IFFT(N(u)))
                            N!(k1X, k1Y, h);
    ifft_plan!(k1X);                                ifft_plan!(k1Y)
    @devec k1X[:] = disp .* k1X;                    @devec k1Y[:] = disp .* k1Y
    fft_plan!(k1X);                                 fft_plan!(k1Y)
        
                            # k2 = N(u1 + k1/2)
    BLAS.axpy!(n, 0.5 + 0.im, k1X, 1, k2X, 1);      BLAS.axpy!(n, 0.5 + 0.im, k1Y, 1, k2Y, 1)
                            N!(k2X, k2Y, h);
    
                            # k3 = N(u1 + k2/2)
    BLAS.axpy!(n, 0.5 + 0.im, k2X, 1, k3X, 1);      BLAS.axpy!(n, 0.5 + 0.im, k2Y, 1, k3Y, 1)
                            N!(k3X, k3Y, h);
    
                            # k4 = N(FFT(D * IFFT(u1 + k3)))
    BLAS.axpy!(n, 1. + 0.im, k3X, 1, k4X, 1);       BLAS.axpy!(n, 1. + 0.im, k3Y, 1, k4Y, 1)
    ifft_plan!(k4X);                                ifft_plan!(k4Y)
    @devec k4X[:] = disp .* k4X;                    @devec k4Y[:] = disp .* k4Y           
    fft_plan!(k4X);                                 fft_plan!(k4Y)
                            N!(k4X, k4Y, h);
    
                            # res = FFT(D * IFFW(u1 + 1/6 k1 + 1/3 (k2 + k3)) ) + 1/6 k4
    BLAS.axpy!(n, 1/6 + 0.im, k1X, 1, ufX, 1);      BLAS.axpy!(n, 1/6 + 0.im, k1Y, 1, ufY, 1)
    BLAS.axpy!(n, 1/3 + 0.im, k2X, 1, ufX, 1);      BLAS.axpy!(n, 1/3 + 0.im, k2Y, 1, ufY, 1)
    BLAS.axpy!(n, 1/3 + 0.im, k3X, 1, ufX, 1);      BLAS.axpy!(n, 1/3 + 0.im, k3Y, 1, ufY, 1)
    ifft_plan!(ufX);                                ifft_plan!(ufY)
    @devec ufX[:] = disp .* ufX;                    @devec ufY[:] = disp .* ufY
    fft_plan!(ufX);                                 fft_plan!(ufY)           
    BLAS.axpy!(n, 1/6 + 0.im, k4X, 1, ufX, 1);      BLAS.axpy!(n, 1/6 + 0.im, k4Y, 1, ufY, 1)
end

function N_simple!(uA, uB, h, dt, gamma, _uabs2A, _uabs2B)
    n = length(u)
    # maybe include this into loop too
    map!(Abs2Fun(), _uabs2A, uA);                   map!(Abs2Fun(), _uabs2B, uB)            

    k = 1im * h * gamma
    @simd for i = 1:n
        @inbounds begin
            uA[i] = k * ((_uabs2A[i] + (2/3)*_uabs2B[i]) * uA[i] + (1/3)*_u2B[i]*_u2B[i] * conj(uA[i]))
            uB[i] = k * ((_uabs2B[i] + (2/3)*_uabs2A[i]) * uB[i] + (1/3)*_u2A[i]*_u2A[i] * conj(uB[i]))
        end
    end
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