# split-step Fourier method

function ssfm(u, L, h, t_grid, w_grid,
              alpha, beta, gamma, steep, t_raman,
              fft_plan!, ifft_plan!)

    nz = int(L/h)
    dt = t_grid[end] - t_grid[end - 1]
    _u = similar(u)
    _uabs2 = similar(u)
    _duabs2 = similar(u)
    _du = similar(u)

    d_exp = dispersion_exponent(w_grid, alpha, beta)
    disp_full = exp(h * d_exp)
    disp_half = exp(h/2. * d_exp)
    disp_minus_half = exp(-h/2. * d_exp)

    # scheme 1/2D => [N, D] loop, -1/2D
    D!(u, disp_half, fft_plan!, ifft_plan!)
    @profile @time for i in 1:nz
        N_ssfm_raman!(u, _uabs2, _duabs2, _du, h, dt, gamma, steep, t_raman)
        D!(u, disp_full, fft_plan!, ifft_plan!)
    end
    D!(u, disp_minus_half, fft_plan!, ifft_plan!)
end

function D!(u, disp, fft_plan!, ifft_plan!)
    ifft_plan!(u)
    @devec u[:] = disp .* u
    fft_plan!(u)
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