alpha=-0.0;
beta2=-5.92e-20;%1.1e-27;
%beta2=0;
beta3=2.98e-34;
gamma=1;%.1;
T0=28.0e-15;
P0=3.01e8;
shape=0;
chirp0=0;
wavelength=2.64e-6;
%steep=(2*pi*3e8/wavelength)^-1;% T0=steep/0.01
steep=0;
t_raman=2.80e-15;
%t_raman=0;
chirp=0;
time_window=50; %22;
L=5.3e-8;

tic;
[t w u0 U0 phi0 uL UL phiL tPlot zPlot uPlot UPlot]=rk4ip_mex(alpha, beta2, beta3,...
    gamma,steep,t_raman,T0,P0,shape,chirp0,time_window, L,12,1.e-6,1.e-6);
toc
