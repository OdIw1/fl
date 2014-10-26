# Runge-Kutta 4th order in the Interaction Picture methods

DEBUG = true

# maybe i should introduce separate bethas for X and Y axes
@eval function rk4ip_vec!(uX, uY, L, h, t, w, 
                     alpha, betha, dbetha, gamma, g, g_bandwidth, saturation_energy,
                     fft_plan!, ifft_plan!, nt_plot=2^8, nz_plot=2^8)
    z = 0.
    n = length(uX)
    T = (t[end] - t[1]) / 2
    dt = (t[end] - t[1]) / (n - 1)

    n_steps = n_steps_rejected = 0
    steps = Float64[]
    err_prev = 1.

    $([:($a = similar(uX)) for a in 
        postfix_vars(["X", "Y"], [:U, :_u1, :_k1, :_k2, :_k3, :_k4, :_uabs2,
                     :_du, :_ue_cplx, :u_full, :u_half, :u_half2])]...)

    N! = let dt = dt, dbetha = dbetha, gamma = gamma, _uabs2X = _uabs2X, _uabs2Y = _uabs2Y
        (uA_, uB_, h_, z_) -> N_simple!(uA_, uB_, h_, z_, dt, dbetha, gamma, _uabs2X, _uabs2Y)
    end

    d_no_gain = dispersion_without_gain(w, alpha, betha)
    g_spec = gain_spectral_factor(w, g_bandwidth)
    d_exp = similar(uX)
    dispersion_exp!(d_exp, d_no_gain, g_spec, uX, uY, dt, g, saturation_energy)
    # d_exp = dispersion_exponent(w, alpha, betha)

    disp_full = exp(h/2 * d_exp)
    disp_half = exp(h/4 * d_exp)

    # prepare plotting
    do_plot = (nt_plot != 0 && nz_plot !=0)
    dz_plot = L / (nt_plot-1)
    t_plot_ind = round(linspace(1, n, nt_plot))
    $([:($a = zeros(Complex{Float64}, nz_plot, nt_plot)) for a in [:u_plotX, :u_plotY, :U_plotX, :U_plotY]]...)
    
    handle_plot! = let t_plot_ind = t_plot_ind, uX = uX, uY = uY, UX = UX, UY = UY,
                       u_plotX = u_plotX, u_plotY = u_plotY, U_plotX = U_plotX, U_plotY = U_plotY,
                       ifft_plan! = ifft_plan!, T = T
        (i_plot_) -> handle_plot_data!(i_plot_, t_plot_ind, uX, uY, UX, UY,
                                        u_plotX, u_plotY, U_plotX, U_plotY, ifft_plan!, T)
    end
    
    i_plot = 1
    do_plot && handle_plot!(i_plot)

    @time @profile while z < L
        # full step
        rk4ip_step!(uX, uY, u_fullX, u_fullY, h, disp_full, N!, z,
                    fft_plan!, ifft_plan!,
                    _u1X, _u1Y, _k1X, _k1Y, _k2X, _k2Y, _k3X, _k3Y, _k4X, _k4Y)
        # 2 half-steps
        rk4ip_step!(uX, uY, u_halfX, u_halfY, h/2, disp_half, N!, z,
                    fft_plan!, ifft_plan!,
                    _u1X, _u1Y, _k1X, _k1Y, _k2X, _k2Y, _k3X, _k3Y, _k4X, _k4Y)
        rk4ip_step!(u_halfX, u_halfY, u_half2X, u_half2Y, h/2, disp_half, N!, z + h/2,
                    fft_plan!, ifft_plan!,
                    _u1X, _u1Y, _k1X, _k1Y, _k2X, _k2Y, _k3X, _k3Y, _k4X, _k4Y)
        
        errX = integration_error_global(u_fullX, u_half2X, _ue_cplxX)
        errY = integration_error_global(u_fullY, u_half2Y, _ue_cplxY)
        err = sqrt((sqr(errX) + sqr(errY)) / 2)

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
            DEBUG && mod(n_steps, 100) == 0 && @show (n_steps, z, h)    

            h *= scale_step_ok(err, err_prev)
            h = min(L - z, h)
            hd2 = h/2
            hd4 = h/4

            err_prev = err
            BLAS.blascopy!(n, u_half2X, 1, uX, 1);      BLAS.blascopy!(n, u_half2Y, 1, uY, 1);

            dispersion_exp!(d_exp, d_no_gain, g_spec, uX, uY, dt, g, saturation_energy)
            @devec disp_full[:] = exp(hd2 .* d_exp)
            @devec disp_half[:] = exp(hd4 .* d_exp)

            if do_plot && z >= i_plot * dz_plot
                # println("step: $n_steps, z: $z")
                i_plot += 1
                handle_plot!(i_plot)            
            end
        end
    end

    return (u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps)
end

function handle_plot_data!(i_plot, t_plot_ind, uX, uY, UX, UY,
                           u_plotX, u_plotY, U_plotX, U_plotY, ifft_plan!, T)
    u_plotX[i_plot,:] = uX[t_plot_ind];                 u_plotY[i_plot,:] = uY[t_plot_ind]             
    spectrum!(uX, UX, ifft_plan!, T);                   spectrum!(uY, UY, ifft_plan!, T)
    U_plotX[i_plot,:] = UX[t_plot_ind];                 U_plotY[i_plot,:] = UY[t_plot_ind]
end

function dispersion_without_gain(w, alpha, betha)
    # pay attention to the order of Fourier transforms that determine
    # the sign of differentiation operator
    res = zeros(Complex{Float64}, length(w))
    res -= alpha/2
    for k in 1:length(betha)
        res += (1.im / factorial(k+1) * betha[k]) * w.^(k+1) # * (-1)^(k+1)
    end
    return res
end

function gain_spectral_factor(w, g_bandwidth)
    1 ./ (1 + w.^2 / g_bandwidth^2)
end

function gain_saturated(uX, uY, dt, g, saturation_energy)
    n = length(uX)
    EX = sqr(BLAS.nrm2(n, uX, 1)) * dt;             EY = sqr(BLAS.nrm2(n, uY, 1)) * dt
    energy = EX + EY
    0.5g / (1. + energy / saturation_energy)
end

function dispersion_exp!(d_exp, d_no_gain, g_spec, uX, uY, dt, g, saturation_energy)
    g_sat = gain_saturated(uX, uY, dt, g, saturation_energy)
    @devec d_exp[:] = d_no_gain + g_sat .* g_spec
end

function rk4ip_step!(uX, uY, ufX, ufY, h, disp, N!, z, fft_plan!, ifft_plan!,
                     u1X, u1Y, k1X, k1Y, k2X, k2Y, k3X, k3Y, k4X, k4Y)
    n = length(uX)
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
                            N!(k1X, k1Y, h, z);
    ifft_plan!(k1X);                                ifft_plan!(k1Y)
    @devec k1X[:] = disp .* k1X;                    @devec k1Y[:] = disp .* k1Y
    fft_plan!(k1X);                                 fft_plan!(k1Y)
        
                            # k2 = N(u1 + k1/2)
    BLAS.axpy!(n, 0.5 + 0.im, k1X, 1, k2X, 1);      BLAS.axpy!(n, 0.5 + 0.im, k1Y, 1, k2Y, 1)
                            N!(k2X, k2Y, h, z);
    
                            # k3 = N(u1 + k2/2)
    BLAS.axpy!(n, 0.5 + 0.im, k2X, 1, k3X, 1);      BLAS.axpy!(n, 0.5 + 0.im, k2Y, 1, k3Y, 1)
                            N!(k3X, k3Y, h, z);
    
                            # k4 = N(FFT(D * IFFT(u1 + k3)))
    BLAS.axpy!(n, 1. + 0.im, k3X, 1, k4X, 1);       BLAS.axpy!(n, 1. + 0.im, k3Y, 1, k4Y, 1)
    ifft_plan!(k4X);                                ifft_plan!(k4Y)
    @devec k4X[:] = disp .* k4X;                    @devec k4Y[:] = disp .* k4Y           
    fft_plan!(k4X);                                 fft_plan!(k4Y)
                            N!(k4X, k4Y, h, z);
    
                            # res = FFT(D * IFFW(u1 + 1/6 k1 + 1/3 (k2 + k3)) ) + 1/6 k4
    BLAS.axpy!(n, 1/6 + 0.im, k1X, 1, ufX, 1);      BLAS.axpy!(n, 1/6 + 0.im, k1Y, 1, ufY, 1)
    BLAS.axpy!(n, 1/3 + 0.im, k2X, 1, ufX, 1);      BLAS.axpy!(n, 1/3 + 0.im, k2Y, 1, ufY, 1)
    BLAS.axpy!(n, 1/3 + 0.im, k3X, 1, ufX, 1);      BLAS.axpy!(n, 1/3 + 0.im, k3Y, 1, ufY, 1)
    ifft_plan!(ufX);                                ifft_plan!(ufY)
    @devec ufX[:] = disp .* ufX;                    @devec ufY[:] = disp .* ufY
    fft_plan!(ufX);                                 fft_plan!(ufY)           
    BLAS.axpy!(n, 1/6 + 0.im, k4X, 1, ufX, 1);      BLAS.axpy!(n, 1/6 + 0.im, k4Y, 1, ufY, 1)
end

function N_simple!(uX, uY, h, z, dt, dbetha, gamma, _uabs2X, _uabs2Y)
    n = length(uX)
    k = 1im * h * gamma
    B = 2/3
    C = 1/3
    phaseX = C * exp(-2im * dbetha * z);          phaseY = C * exp(2im * dbetha * z)    
    # maybe include this into loop too
    map!(Abs2Fun(), _uabs2X, uX);                   map!(Abs2Fun(), _uabs2Y, uY)                    
    for i = 1:n
        @inbounds _uX = uX[i];                      _uY = uY[i]
        @inbounds uX[i] = k * ((_uabs2X[i] + B*_uabs2Y[i])*_uX + phaseX*(_uY*_uY)*conj(_uX))
        @inbounds uY[i] = k * ((_uabs2Y[i] + B*_uabs2X[i])*_uY + phaseY*(_uX*_uX)*conj(_uY))
    end
end    