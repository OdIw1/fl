function soliton_vec(N=1)
    n = 2^14
    T0 = 56.7e-14
    C0 = 0.
    theta = pi / 4
    T_window = 40

    alpha = 0.
    betha2 = -1.e-26
    dbetha = 0.2 * 2 * abs(betha2) / T0
    gamma = 1e-2
    steep = 0.

    P0 = N^2 * abs(betha2) / (gamma * T0^2)
    L = 10 * T0^2 / abs(betha2)
    # L = pi / 2 * T0^2 / abs(betha2)
    run_vec(n, T_window, alpha, [betha2], dbetha, gamma, L, T0, P0, C0, theta)
end

function Agr3_6_vec(case=0, theta=0.)
    n = 2^12
    T0 = 85.0e-15
    P0 = 10.e3
    C0 = 0.  
    T_window = 20

    alpha = 0.
    b3 = 8.119e-2
    b2 = case == 0 ? 0. : b3 / (T0 / 1.e-12)
    betha = ps_km2s_m([b2, b3])
    gamma = 0.

    L = 5 * T0^3 / abs(betha[2])
    run_vec(n, T_window, alpha, betha, 0., gamma, L, T0, P0, C0, theta, 1)
end

function Agr3_6_waveplate(case=0, theta=0.::Real)
    n = 2^12
    T0 = 85.e-15
    P0 = 10.e3 
    b3 = 8.119e-2
    b2 = case == 0 ? 0. : b3 / (T0 / 1.e-12)
    betha = ps_km2s_m([b2, b3])
    gamma = 0.
    
    @show L = 5 * T0^3 / abs(betha[2])
    @show Ld, Ln, sol_order = pulse_params(T0, P0, betha, gamma)
    f = Fiber(L, 0., betha, gamma)
    p = Pulse(1, T0, P0, 0. , theta, n, 20T0)

    wp = QuarterWavePlate(pi/4)

    outdir = "/mnt/hgfs/VM_shared/out/"
    fout = FileOutput(outdir)

    scheme = LaserElement[]
    push!(scheme, wp, fout, f, fout)
    run_laser_scheme!(p, scheme)
end

function Agr4_15_waveplate(case=0, wp_angle=pi/4)
    n = 2^12
    T0 = 85.e-15
    P0 = 10.e3 
    b3 = 8.119e-2
    b2 = case == 0 ? 0. : b3 / (T0 / 1.e-12)
    betha = ps_km2s_m([b2, b3])
    gamma = abs(betha[2]) / (P0 * T0^3)
    
    @show L = 5 * T0^3 / abs(betha[2])
    @show Ld, Ln, sol_order = pulse_params(T0, P0, betha, gamma)
    f = Fiber(L, 0., betha, gamma)
    p = Pulse(1, T0, P0, 0., 0., n, 20T0)

    # wp = HalfWavePlate(pi/4)
    wp = QuarterWavePlate(wp_angle)

    outdir = "/mnt/hgfs/VM_shared/out/"
    fout = FileOutput(outdir)

    scheme = LaserElement[]
    push!(scheme, wp, fout, f, fout)
    run_laser_scheme!(p, scheme)
end

function Yarutkina13(n_iter=1, a1=pi/6, a2=0, a3=0)
    # fiber parameters are from 13[Yarutkina, Shtyrina]{Opt.Expr} Numerical Modeling of ...
    wl = 1550e-9
    gain_bw_wl = 50e-9
    gain_bw = bandwidth_wl2fr(wl, gain_bw_wl)
    La = 2.
    Lp = 5.
    t_round = (La + Lp) * 1.47 / 3.e8
    sat_e = 20.e-3 * t_round

    fiber_active = Fiber(La, 0., fs_mm2s_m([76.9, 168.]), 9.32e-3, 10^(5.4/10), gain_bw, sat_e)
    fiber_passive = Fiber(Lp, 10^(0.2/10) * 1.e-3, fs_mm2s_m([4.5, 109]), 2.1e-3)

    PC1 = QuarterWavePlate(a1)

    J2 = HalfWavePlate(a2)
    J3 = QuarterWavePlate(a3)
    J4 = Polarizer()
    PC2 = J4 * J3 * J2

    # seed pulse params
    n = 2^14
    T0 = 1.e-10
    P0 = 1.e-2
    T = 30T0
    p = Pulse(0, T0, P0, 0., 0., n, T)

    F = SpectralFilter(0.8)

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    # outdir = "/mnt/hgfs/VM_shared/out/"
    foutPC1 = FileOutput(outdir, "PC1")
    foutA   = FileOutput(outdir, "A")
    foutP   = FileOutput(outdir, "P")
    foutPC2 = FileOutput(outdir, "PC2")
    # run
    laser = LaserElement[PC1, foutPC1, fiber_active, F, foutA, fiber_passive, foutP, PC2, foutPC2]   
    run_laser_scheme!(p, laser, n_iter)
end


function Yarutkina13scalar(n_iter=1)
    # fiber parameters are from 13[Yarutkina, Shtyrina]{Opt.Expr} Numerical Modeling of ...
    wl = 1550e-9
    gain_bw_wl = 50e-9
    gain_bw = bandwidth_wl2fr(wl, gain_bw_wl)
    La = 2.
    Lp = 100.
    t_round = (La + Lp) * 1.47 / 3.e8
    sat_e = 20.e-3 * t_round

    betha_a = fs_mm2s_m([76.9, 168.])
    gamma_a = 9.32e-3
    betha_p = fs_mm2s_m([4.5, 109])
    gamma_p = 2.1e-3
    fiber_active = Fiber(La, 0., betha_a, gamma_a, 10^(5.4/10), gain_bw, sat_e, 500)
    fiber_passive = Fiber(Lp, 10^(0.2/10) * 1.e-3, betha_p, gamma_p, 1000)

    sa = SaturableAbsorber(0.1, 3.69)
    cp = Coupler(0.9)

    # seed pulse params
    n = 2^13
    T0 = 1.e-10
    P0 = 1.e-2
    T = 5.e-9
    p = Pulse(0, T0, P0, 0., 0., n, T)
    @show pulse_params(T0, P0, betha_a, gamma_a)

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    # outdir = "/mnt/hgfs/VM_shared/out/"
    foutA   = FileOutput(outdir, "A", true)
    foutS   = FileOutput(outdir, "S", true)
    foutP   = FileOutput(outdir, "P", true)

    es = PulseSensor()
    pol = Polarizer()
    # run
    laser = LaserElement[pol, foutP, es, fiber_active, es, foutA, cp, sa, foutS, es, fiber_passive]   
    run_laser_scheme!(p, laser, n_iter)
end