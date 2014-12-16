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

function flat_gain(adaptive_step=true, n_iter=100)
    n = 10
    T = 1
    t = t_grid(n, T)
    w = w_grid(n, T)
    pa = [1. + 0.im for i = 1:n]
    p = Pulse(pa, 0. * pa, t, w)
    alpha = 100.
    gain = 1.e3
    E_sat = 10.e-9
    E = E_sat*(gain/alpha -1)
    print_with_color(:green, "Energy should be $E\n")
    f = Fiber(1., alpha, [0.], 0., gain, 1.e100, E_sat, 1000, adaptive_step)
    E = PulseSensor("F")
    run_laser_scheme!(p, LaserElement[E, f], n_iter)
end


function Yarutkina13_vec(n_iter=1, a1=pi/6, a2=0, a3=0)
    # fiber parameters are from 13[Yarutkina, Shtyrina]{Opt.Expr} Numerical Modeling of ...
    wl = 1550e-9
    gain_bw_wl = 50e-9
    gain_bw = bandwidth_wl2fr(wl, gain_bw_wl)
    La = 2.
    Lp = 30.
    t_round = (La + Lp) * 1.47 / 3.e8
    sat_e = 20.e-3 * t_round

    betha_a = fs_mm2s_m([76.9, 168.])
    gamma_a = 9.32e-3
    betha_p = fs_mm2s_m([4.5, 109])
    gamma_p = 2.1e-3

    fiber_active = Fiber(La, 0., betha_a, gamma_a, 10^(5.4/10), gain_bw, sat_e, 2000, true)
    fiber_passive = FiberPassive(Lp, 10^(0.2/10) * 1.e-3, betha_p, gamma_p, 5000, true)

    PC1 = QuarterWavePlate(a1)

    J2 = HalfWavePlate(a2)
    J3 = QuarterWavePlate(a3)
    J4 = Polarizer()
    PC2 = J4 * J3 * J2

    coupler = Coupler(0.9)

    # seed pulse params
    n = 2^15
    T0 = 1.e-9
    P0 = 1.e-10
    T = 5.e-9
    #p = Pulse(1, T0, P0, 0., 0., n, T)
    p = NoisePulse(1.e-10, 3.e-11, n, T)
    @show bandwidth_wl(n, T, wl)

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    # outdir = "/mnt/hgfs/VM_shared/out/"
    foutPC1 = FileOutput(outdir, "PC1")
    foutA   = FileOutput(outdir, "A")
    foutP   = FileOutput(outdir, "P")
    foutPC2 = FileOutput(outdir, "PC2")

    S = PulseSensor()
    
    laser = LaserElement[PC1, foutPC1, S, fiber_active, foutA, S, fiber_passive, foutP, S,
                         PC2, foutPC2, S, coupler]   
    run_laser_scheme!(p, laser, n_iter)
end

function Chong08(adaptive_step=false, n_iter=9999)
    # 08[Chong, Renninger, Wise]{J.Opt.Soc.Am.B,25,2} Properties of normal-dispersion
    #   femtosecond lasers
    wl0 = 1030.e-9 # Ytterbium laser
    bw_wl = 50.e-9 # half of FWHM in article
    bw_fr = bandwidth_wl2fr_derivative(wl0, bw_wl)

    Lp1 = 3.
    La = 0.6
    Lp2 = 1.
    steps_per_meter = 500

    t_round = (Lp1 + La + Lp2) * 1.47 / 3.e8
    betha = fs_cm2s_m([230.])
    gamma = 4.7e-3
    
    Fp1 = FiberPassive(Lp1, 0., betha, gamma, int(steps_per_meter * Lp1), adaptive_step)
    Fp2 = FiberPassive(Lp2, 0., betha, gamma, int(steps_per_meter * Lp2), adaptive_step)
    
    E_sat = 0.25e-9 # varied from 0.25 nJ to 6 nJ
    gain = 10^(30./10) * La
    Fa = Fiber(La, 0., betha, gamma, gain, bw_fr, E_sat,
               int(steps_per_meter * La), adaptive_step)

    SA = SaturableAbsorber(0.7, 1.0e3) # varied from 0.1 to 2.4 kW
    coupler = Coupler((1-0.7)*(1-0.1)) # p.142, 70% out and 10% loss after that
    SF = GaussianSpectralFilter(wl0, 10.e-9) # 8 to 25 nm
    polarizer = Polarizer()

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    o1   = FileOutput(outdir, "1", ONLY_X)
    o2   = FileOutput(outdir, "2", ONLY_X)
    o3   = FileOutput(outdir, "3", ONLY_X)
    o4   = FileOutput(outdir, "4", ONLY_X)
    o5   = FileOutput(outdir, "5", ONLY_X)
    o6   = FileOutput(outdir, "6", ONLY_X)

    E1 = PulseSensor("SMF1")
    E2 = PulseSensor("gain")
    E3 = PulseSensor("SMF2")
    E4 = PulseSensor("SA")
    E5 = PulseSensor("SF")
    E6 = PulseSensor("coupler")

    laser = LaserElement[polarizer, Fp1, E1, o1, Fa, E2, o2, Fp2, E3, o3,
                         SA, E4, o4, SF, E5, o5, coupler, E6, o6]

    n = 2^14
    T = 5.e-11
    # p = Pulse(1, 1.e-12, 1.e-10, 0., 0., n, T)
    p = WhiteNoisePulse(1.e-10, n, T)
    @show bandwidth_wl(n, T, wl0)

    run_laser_scheme!(p, laser, n_iter)
end

function Abdelalim08(gain_dB=65,adaptive_step=false, n_iter=9999)
    # 08[Abdelalim et al.]{Opt.Expr,17,4} Properties and stability limits of
    # an optimized mode-locked Yb-doped femtosecond fiber laser
    wl0 = 1030.e-9 # Ytterbium laser
    bw_wl = 40.e-9 
    bw_fr = bandwidth_wl2fr_derivative(wl0, bw_wl)

    La = 0.6
    steps_per_meter = 5000

    t_round = La * 1.47 / 3.e8
    alpha = 0.04
    betha = fs_cm2s_m([240.])
    gamma = 5.0e-3 
   
    E_sat = 3.0e-9
    gain = gain_dB !=0. ? 10^(gain_dB/10) : 0
    Fa = Fiber(La, alpha, betha, gamma, gain, bw_fr, E_sat,
               int(steps_per_meter * La), adaptive_step)

    SA = SaturableAbsorber(0.5, 200.0e3)
    coupler = Coupler(0.1) # 10 dB loss
    SF = GaussianSpectralFilter(wl0, 10.e-9)
    polarizer = Polarizer()

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    o1   = FileOutput(outdir, "1", ONLY_X)
    o2   = FileOutput(outdir, "2", ONLY_X)
    o3   = FileOutput(outdir, "3", ONLY_X)

    E1 = PulseSensor("gain")
    E2 = PulseSensor("SF")
    E3 = PulseSensor("coupler")

    laser = LaserElement[polarizer, Fa, E1, o1, SA, E2, o2, coupler, E3, o3]

    n = 2^15
    T = 10.e-11
    # p = Pulse(1, 1.e-12, 1.e-10, 0., 0., n, T)
    p = WhiteNoisePulse(1.e-20, n, T)
    @show bandwidth_wl(n, T, wl0)

    run_laser_scheme!(p, laser, n_iter)
end

function Smirnov12(adaptive_step=false, n_iter=9999)
    # 12 [Smirnov, Kobtsev et al.]{Opt.Expr}
    # Three key regimes of single pulse generation per
    # round trip of all-normal-dispersion fiber lasers
    # mode-locked with nonlinear polarization
    # rotation
    wl0 = 1030.e-9 # Ytterbium laser
    bw_wl = 50.e-9 # half of FWHM in article
    bw_fr = bandwidth_wl2fr_derivative(wl0, bw_wl)

    Lp = 33.
    La = 7.

    steps_per_meter = 500

    t_round = (Lp + La) * 1.47 / 3.e8
    betha = ps_km2s_m([23])
    gamma = 5.0e-3
    
    Fp = FiberPassive(Lp, 0., betha, gamma, int(steps_per_meter * Lp), adaptive_step)
    
    P_sat = 52e-3
    E_sat = P_sat * t_round
    gain = 10^(5.4/10) * La
    Fa = Fiber(La, 0., betha, gamma, gain, bw_fr, E_sat,
               int(steps_per_meter * La), adaptive_step)

    SA = SaturableAbsorber(0.5, 50.0e3) # 
    coupler = Coupler(0.8) #
    SF = GaussianSpectralFilter(wl0, 100.e-9) # bylo 50
    polarizer = Polarizer()

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    o1   = FileOutput(outdir, "1", ONLY_X)
    o2   = FileOutput(outdir, "2", ONLY_X)
    o3   = FileOutput(outdir, "3", ONLY_X)
    o4   = FileOutput(outdir, "4", ONLY_X)
    o5   = FileOutput(outdir, "5", ONLY_X)

    E1 = PulseSensor("SMF1")
    E2 = PulseSensor("gain")
    E3 = PulseSensor("SA")
    E4 = PulseSensor("SF")
    E5 = PulseSensor("coupler")

    laser = LaserElement[polarizer, Fp, E1, o1, Fa, E2, o2, 
                         SA, E3, o3, SF, E4, o4, coupler, E5, o5]

    n = 2^14
    T = 1000.e-12
    p = Pulse(1, 1.e-10, 1.e-10, 0., 0., n, T)
    # p = WhiteNoisePulse(1.e-10, n, T)
    @show bandwidth_wl(n, T, wl0)

    run_laser_scheme!(p, laser, n_iter)
end

function Yarutkina13_scalar(adaptive_step=false, n_iter=9999)
    # fiber parameters are from 13[Yarutkina, Shtyrina]{Opt.Expr} Numerical Modeling of ...
    wl0 = 1550.e-9
    gain_bw_wl = 50.e-9
    gain_bw = bandwidth_wl2fr_derivative(wl0, gain_bw_wl)

    La = 2.
    Lp = 30.
    t_round = (La + Lp) * 1.47 / 3.e8
    sat_e = 20.e-3 * t_round

    steps_per_meter = 500

    betha_a = fs_mm2s_m([76.9, 168.])
    betha_p = fs_mm2s_m([4.5, 109])
    # seems like gain and absorption is provided in field, not intensity related units, 
    # so gain and alpha may 2 factor
    Fa = Fiber(La, 0., betha_a, 9.32e-3, 2*10^(5.4/10), gain_bw, sat_e, 
               int(steps_per_meter * Lp), adaptive_step)
    Fp = FiberPassive(Lp, 2*10^(0.2/10) * 1.e-3, betha_p, 2.1e-3,
                      int(steps_per_meter * Lp), adaptive_step)

    SA = SaturableAbsorber(0.3, 3.69e3)
    coupler = Coupler(0.1)
    SF = GaussianSpectralFilter(wl0, 100.e-9) # bylo 50
    polarizer = Polarizer()

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    o1   = FileOutput(outdir, "1", ONLY_X)
    o2   = FileOutput(outdir, "2", ONLY_X)
    o3   = FileOutput(outdir, "3", ONLY_X)
    o4   = FileOutput(outdir, "4", ONLY_X)
    o5   = FileOutput(outdir, "5", ONLY_X)

    E1 = PulseSensor("SMF1")
    E2 = PulseSensor("gain")
    E3 = PulseSensor("SA")
    E4 = PulseSensor("SF")
    E5 = PulseSensor("coupler")

    laser = LaserElement[polarizer, Fa, E1, o1, Fp, E2, o2, 
                         SA, E3, o3, SF, E4, o4, coupler, E5, o5]  

    n = 2^14
    T0 = 1.e-10
    P0 = 1.e-10
    T = 5.e-9
    p = Pulse(1, 1.e-9, 1.e-10, 0., 0., n, T)
    # p = WhiteNoisePulse(1.e-10, n, T)
    
    @show bandwidth_wl(n, T, wl0)
    run_laser_scheme!(p, laser, n_iter)
end

function Felleher14(adaptive_step=false, n_iter=9999)
    # fiber parameters are from 14[Felleher et al.]{Opt.Lett} Chirp pulse
    # formation dynamics in ultra-long mode-locked fiber lasers
    wl0 = 1550.e-9
    gain_bw_wl = 40.e-9
    gain_bw = bandwidth_wl2fr_derivative(wl0, gain_bw_wl)

    La = 2.
    Lp = 750.
    t_round = (La + Lp) * 1.47 / 3.e8
    SCALE = (La + Lp) / 1202.
    sat_e = 200.e-12 * SCALE

    steps_per_meter = 500

    betha = ps_m2s_m([0.018])
    gamma = 3.e-3
    Fa = Fiber(La, 0., betha, gamma, 10^(30./10), gain_bw, sat_e, 
               1000, adaptive_step)
    Fp = FiberPassive(Lp, 0., betha, gamma,
                      5000, adaptive_step)

    SA = SaturableAbsorber(0.4, 6.e1 / SCALE, 0.45)
    coupler = Coupler(1 - 1./ 10^(3./10))
    SF = GaussianSpectralFilter(wl0, 40.e-9) # THIS IS IMPORTANT
    polarizer = Polarizer()

    outdir = mkpath_today("/home/s_koval/vm_shared")
    o1   = FileOutput(outdir, "1", ONLY_X)
    o2   = FileOutput(outdir, "2", ONLY_X)
    o3   = FileOutput(outdir, "3", ONLY_X)
    o4   = FileOutput(outdir, "4", ONLY_X)
    o5   = FileOutput(outdir, "5", ONLY_X)

    E1 = PulseSensor("SMF1")
    E2 = PulseSensor("gain")
    E3 = PulseSensor("SA")
    E4 = PulseSensor("SF")
    E5 = PulseSensor("coupler")

    laser = LaserElement[polarizer, Fp, E1, o1, Fa, E2, o2, 
                         SA, E3, o3, SF, E4, o4, coupler, E5, o5]  

    n = 2^15
    T = 20.e-10
    # p = Pulse(1, 1.e-9, 1.e-10, 0., 0., n, T)
    p = WhiteNoisePulse(1.e-10, n, T)
    
    @show bandwidth_wl(n, T, wl0)
    run_laser_scheme!(p, laser, n_iter)
end