include("rotationFunctions.jl")
include("LTV_observers.jl")
include("Lyapunov_observers.jl")
include("initialize.jl")

function sys(x, tau, R, K)
	w = K.J \ (R' * x[5:7])
	q = x[1:4]
	dq = 0.5 * [0.0 -w'; w -S(w)] * q + K.a * (1.0 / (q' * q) - 1.0) * q
	dL = R * tau

	return [dq; dL]
end

function f(x, p, t)
	R = q2R(x[p[1].list[1:4]])
	tau = [sin(t); cos(pi*t); 0]
	z = p[1].J \ (R' * x[p[1].list[5:7]]) + p[1].b
	
	return [sys(x[p[1].list], tau, R, p[1]);
			Lyapunov_minimal(x[p[2].list], tau, R, p[2]);
			Lyapunov_biased(x[p[3].list], tau, R, z, p[3]);
			Lyapunov_kinematic(x[p[4].list], R, z, p[4]);
			LTV_minimal(x[p[5].list], R, tau, p[5]);
			LTV_biased(x[p[6].list], R, z, tau, p[6]);
			LTV_kinematic(x[p[7].list], R, z, p[7])]
end

function solve_and_plot(tspan, p, w0)
	figs = [plot(); plot()]
	for k in (0, 1)
		x0 = initial_conditions(w0, wh = k * 0.95 * w0, bh = -50 * k * p[1].b)
		prob = ODEProblem(f, x0, tspan, p)
		sol = solve(prob, Vern9(), reltol=1e-15, abstol=1e-15, dense=false)

		font = Plots.font("Computer Modern", 24)
		time, E_L, E_w, E_R = proces_output(sol, p)
		labels = ["Lyapunov minimal" "Lyapunov biased" "Lyapunov kinematic" "LTV minimal" "LTV biased" "LTV kinematic"]
		fig_w = plot(time, E_w, yaxis=:log, lab="", xlabel="Time", ylabel=L"$\|\tilde{\omega}\|$", xtickfont=font,ytickfont=font,xguidefont=font,yguidefont=font)
		fig_R = plot(time, E_R, yaxis=:log, lab="", ylabel=L"$\|u^{eb}\|$", xtickfont=font,ytickfont=font,xguidefont=font,yguidefont=font)
		fig_labels = plot([0,0], zeros(2,6), lab=labels, axis=false, grid=false, legendfont=font)
		fn = plot(axis=false, grid=false)
		figs[k+1] = plot(fn,fig_R,fig_labels,fn,fig_w,fn,fn,fn,fn, layout=grid(3, 3, widths=(0.0,0.6,0.4), heights=(0.5,0.5,0.0)), size=(1440,1060))
	end
	return figs
end
