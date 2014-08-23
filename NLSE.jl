module NLSE
export create_t_grid, create_w_grid

function t_grid(n, T)
	# In MATLAB version 'T' was a scaling factor: time window was (-T0 T, T0 T),
	# where T0 is pulse FWHM;  in Julia time window is (-T, T)
	# also 'n' is no longer binary log of the number of grid points, it is the number itself
	dt = 2 * T / n
	(-n/2 : n/2 - 1) * dt
end

function w_grid(n, T)
	dw = pi / T
	dw * [(0 : n/2 -1), (-n/2 : -1)]
end

function secant_pulse(t, T0)


