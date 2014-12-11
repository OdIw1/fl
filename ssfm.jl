# split-step Fourier method

ssfm!(p::Pulse, f::Fiber) = 
    ssfm!(p.uX, p.t, p.w, f.L, 1e-6f.L, f.max_steps, f.adaptive_step,
          f.alpha, f.betha, f.gamma, 0., 0., 
          f.gain, f.gain_bandwidth, f.saturation_energy,
          p.fft_plan!, p.ifft_plan!)

function ssfm!(u, t, w, L, h0, max_steps, adaptive_step,
               alpha, betha, gamma, steep, t_raman,
               gain, gain_bw, saturation_energy,
               fft_plan!, ifft_plan!)
    h = L / max_steps
    dt = calc_dt(t)

    _uabs2 = similar(u)
    _d_exp = similar(u)
    # _duabs2 = similar(u)
    # _du = similar(u)

    d_no_gain = D_exp_no_gain(w, alpha, betha)
    g_spec_factor = gain_spectral_factor(w, gain, gain_bw)

    # scheme 1/2D => [N, D] loop, -1/2D
    D_ssfm!(u, 0.5h, _d_exp, d_no_gain, g_spec_factor, dt, saturation_energy,
       fft_plan!, ifft_plan!)
    @time for i in 1:max_steps
        # N_ssfm_raman!(u, _uabs2, _duabs2, _du, h, dt, gamma, steep, t_raman)
        N_ssfm_simple!(u, _uabs2, h, gamma)
        D_ssfm!(u, h, _d_exp, d_no_gain, g_spec_factor, dt, saturation_energy,
           fft_plan!, ifft_plan!)
    end
    D_ssfm!(u, -0.5h, _d_exp, d_no_gain, g_spec_factor, dt, saturation_energy,
       fft_plan!, ifft_plan!)
end

function D_ssfm!(u, h, _d_exp, d_no_gain, g_spec_factor, dt, saturation_energy,
            fft_plan!, ifft_plan!)
    g_sat = gain_saturation_scal(u, dt, saturation_energy)
    @devec _d_exp[:] = exp(h .* (d_no_gain + g_sat .* g_spec_factor))

    ifft_plan!(u)
    @devec u[:] = _d_exp .* u
    fft_plan!(u)
end

function N_ssfm_simple!(u, _uabs2, h, gamma)
    n = length(u)
    map!(Abs2Fun(), _uabs2, u)
    k = 1.im * h * gamma
    BLAS.scal!(n, k, _uabs2, 1)
    @devec u[:] = exp(_uabs2) .* u
end

function N_ssfm_raman!(u, _uabs2, _duabs2, _du, h, dt, gamma, steep, t_raman)
    n = length(u)
    map!(Abs2Fun(), _uabs2, u)
    BLAS.blascopy!(length(u), _uabs2, 1, _duabs2, 1)
    df!(_duabs2, _du, dt)
    BLAS.axpy!(n, -t_raman + 0.im, _du, 1, _uabs2, 1)
    k = 1im * h * gamma
    BLAS.scal!(n, k, _uabs2, 1)
    @devec u[:] = exp(_uabs2) .* u
end    