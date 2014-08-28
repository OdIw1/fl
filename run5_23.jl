using NLSE

n = 2^9
alpha = 0.
beta = (-5.92e-20, 2.98e-34)
gamma = 1.
T0 = 28.0e-15
P0 = 3.01e8
C0 = 0.
wl = 2.64e-6
steep = 0
t_raman = 2.80e-15
#t_raman = 0.
T_window = 50
L = 5.3e-8

function run5_23()
    T = T0 * T_window
    t = t_grid(n, T)
    w = w_grid(n, T)
    u0 = secant_pulse(T0, P0, C0, t)

    fft_tmp = ones(Complex{Float64}, n)
    fft_plan = plan_fft(fft_tmp, (1,), FFTW.MEASURE)
    ifft_plan = plan_ifft(fft_tmp, (1,), FFTW.MEASURE)

    Res = integrate_RK4IP(u0, t, w, fft_plan, ifft_plan, L, 1e-4L, alpha, beta, gamma, steep, t_raman)

    # println("$(u0)")
    return u0
end
