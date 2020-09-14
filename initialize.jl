function xM2M(d)
	Symmetric(reshape(1:d^2,d,d) - repeat((collect(0:d-1)' .* (collect(0:d-1)' .+ 1)) .>> 1, d),:L)[:] .+ d;
end

function M2xM(d)
	UnitLowerTriangular(trues(d,d));
end

function param(J)
	I3 = Matrix(1.0I, 3, 3);
	I6 = Matrix(1.0I, 6, 6);
	I9 = Matrix(1.0I, 9, 9);
	I12 = Matrix(1.0I, 12, 12);
	Jh = J;
	p1 = (J=J, a=1, b=[1.0;2.0;3.0], list=1:7, R=Matrix(1.0I,3,2), Ntau=0e-6, Nz=0e-8, NR=1e-4, dt=1e-4, seed=0);
	p2 = (J=Jh, a=1, k1=1, k2=2, k3=0, v1=I3[:,1], v2=I3[:,2], v3=I3[:,3], D=2*I3, G=20*I3, list=8:14);
	p3 = (J=Jh, a=1, k1=1, k2=2, k3=0, v1=I3[:,1], v2=I3[:,2], v3=I3[:,3], D=I6, X=I3, G=10*I6, list=15:24);
	p4 = (J=Jh, a=1, k1=1, k2=2, k3=0, v1=I3[:,1], v2=I3[:,2], v3=I3[:,3], D=10*I3, G=40*I3, list=25:31);
	Gj = [S([1 0 0]) S([0 1 0]) S([0 0 1])]' / Jh;
	G = [S([1 0 0]) S([0 1 0]) S([0 0 1])]';
	W = Matrix(I,6,9);
	p5 = (J=Jh, G=Gj, W=W, V=I9, Winv=I6, d=1, xM2M=xM2M(9), M2xM=M2xM(9), list=32:85);
	p6 = (J=Jh, G=Gj, W=W, V=I12, Winv=I9, d=1, xM2M=xM2M(12), M2xM=M2xM(12), list=86:175);
	p7 = (J=Jh, G=G, W=W, V=I9, Winv=I6, d=1, xM2M=xM2M(9), M2xM=M2xM(9), list=176:229);
	
	return (p1, p2, p3, p4, p5, p6, p7);
end

function initial_conditions(w0; q0 = [0.5;0.5;0.5;0.5], wh = zeros(3), qh = [1.;0;0;0], bh = zeros(3))
	p = param(Matrix(I,3,3));
	Rh = p[5].W * q2R(qh)[:];
	M9 = Matrix(1.0I,9,9);
	M12 = Matrix(1.0I,12,12);

	return [q0; w0;
			qh; wh;
			qh; wh; bh;
			qh; bh;
			Rh; wh; M9[p[5].M2xM];
			Rh; wh; bh; M12[p[6].M2xM];
			Rh; bh; M9[p[7].M2xM]];
end

function proces_output(sol, p, step=1, tstart=0)
	len = length(sol.t)
	list = sum(sol.t.<=tstart) : step : len
	time = sol.t[list]
	E_L = zeros(6,length(time))
	E_L[1,:] = sqrt.(sum((sol[p[2].list[5:7],list] - sol[p[1].list[5:7],list]).^2, dims=1))
	E_L[2,:] = sqrt.(sum((sol[p[3].list[5:7],list] - sol[p[1].list[5:7],list]).^2, dims=1))
	E_L[4,:] = sqrt.(sum((sol[p[5].list[7:9],list] - sol[p[1].list[5:7],list]).^2, dims=1))
	E_L[5,:] = sqrt.(sum((sol[p[6].list[7:9],list] - sol[p[1].list[5:7],list]).^2, dims=1))
	E_w = zeros(6,length(time))
	if isdefined(sol, :W)
		E_R = zeros(7,length(time))
		Random.seed!(p[1].seed)
	else
		E_R = zeros(6,length(time))
	end
	for n = 1:length(time)
		m = list[n]
		q = sol[p[1].list[1:4],m]
		R = q2R(q)
		w = p[1].J \ (R' * sol[p[1].list[5:7],m])
		if isdefined(sol, :W)
			noise = randn(12)
			d = p[1].Nz * noise[4:6]
			E_R[7,n] = R_angle(q, wahba(R * p[1].R + p[1].NR * [noise[7:9] noise[10:12]], p[1].R))
		else
			d = zeros(3)
		end
		E_L[3,n] = norm(R * p[1].J * (p[1].b - sol[p[4].list[5:7],m] + d))
		E_L[6,n] = norm(R * p[1].J * (p[1].b - sol[p[7].list[7:9],m] + d))
		for k = 1:3
			E_R[k,n] = quat_angle(q, sol[p[k+1].list[1:4],m])
			E_R[k+3,n] = R_angle(q, wahba(reshape(sol[p[k+4].list[1:6],m],3,2),p[1].R))
		end
		E_w[1,n] = norm(p[2].J \ (R' * sol[p[2].list[5:7],list[n]]) - w)
		E_w[2,n] = norm(p[3].J \ (R' * sol[p[3].list[5:7],list[n]]) - w)
		E_w[3,n] = norm(p[1].b - sol[p[4].list[5:7],m] + d)
		E_w[4,n] = norm(p[5].J \ (R' * sol[p[5].list[7:9],list[n]]) - w)
		E_w[5,n] = norm(p[6].J \ (R' * sol[p[6].list[7:9],list[n]]) - w)
		E_w[6,n] = norm(p[1].b - sol[p[7].list[7:9],m] + d)
	end
	E_L[E_L .< eps()] .= eps()
	E_w[E_w .< eps()] .= eps()
	E_R[E_R .< eps()] .= eps()
	
	return time, E_L', E_w', E_R'
end

function proces_output_Erjen(sol, step=1, tstart=0)
	len = length(sol.t)
	list = sum(sol.t.<=tstart) : step : len
	time = sol.t[list]
	E_L = zeros(6,length(time))
	E_L[1:3,:] = sol[p[2].list[5:7],list] - sol[p[1].list[5:7],list]
	E_L[4:6,:] = sol[19:21,list]
	E_w = zeros(6,length(time))
	E_q = zeros(8,length(time))
	for n = 1:length(time)
		m = list[n]
		q = sol[p[1].list[1:4],m]
		qh = sol[p[2].list[1:4],m]
		R = q2R(q)
		w = p[1].J \ (R' * sol[p[1].list[5:7],m])
		E_w[1:3,n] = p[2].J \ (R' * sol[p[2].list[5:7],m]) - w
		E_w[4:6,n] = p[2].J \ (R' * sol[19:21,m])
		E_q[1:4,n] = [q[1] * qh[1] + q[2:4]' * qh[2:4]; q[1] * qh[2:4] - qh[1] * q[2:4] + cross(q[2:4], qh[2:4])] / norm(q) / norm(qh)
		E_q[5:8,n] = sol[15:18,m] / norm(sol[15:18,m])
	end
	
	return time, E_L', E_w', E_q'
end
