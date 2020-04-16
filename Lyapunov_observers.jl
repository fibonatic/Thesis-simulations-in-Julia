function Lyapunov_minimal(x, tau, R, K)
	q = x[1:4]
	Rh = q2R(q)
	C = K.k1 * cross(R' * K.v1, Rh' * K.v1) + 
		K.k2 * cross(R' * K.v2, Rh' * K.v2) + 
		K.k3 * cross(R' * K.v3, Rh' * K.v3)
	w = K.J \ (R' * x[5:7]) + K.D * C
	
	dqh = 0.5 * [0.0 -w'; w -S(w)] * q + K.a * (1.0 / (q' * q) - 1.0) * q
	dLh = R * tau + K.G * R * (K.J \ C)
	
	return [dqh; dLh]
end

function Lyapunov_biased(x, tau, R, z, K)
	q = x[1:4]
	Rh = q2R(q)
	C = K.k1 * cross(R' * K.v1, Rh' * K.v1) + 
		K.k2 * cross(R' * K.v2, Rh' * K.v2) + 
		K.k3 * cross(R' * K.v3, Rh' * K.v3)
	z_til = K.J \ (R' * x[5:7]) + x[8:10] - z
	e = [C; z_til]
	eye = Matrix(I,3,3)
	del = -K.G * [R / K.J zeros(3,3); eye eye] * [K.X * C + K.D[4:6,4:6] * z_til; C]
	w = z - x[8:10] + K.D[1:3,1:3] * C - (2 * K.D[1:3,4:6] + K.X') * z_til
	
	dqh = 0.5 * [0.0 -w'; w -S(w)] * q + K.a * (1.0 / (q' * q) - 1.0) * q
	dLh = R * tau + del[1:3]
	dbh = del[4:6]
	
	return [dqh; dLh; dbh]
end

function Lyapunov_kinematic(x, R, z, K)
	q = x[1:4]
	Rh = q2R(q)
	C = K.k1 * cross(R' * K.v1, Rh' * K.v1) + 
		K.k2 * cross(R' * K.v2, Rh' * K.v2) + 
		K.k3 * cross(R' * K.v3, Rh' * K.v3)
	w = z - x[5:7] + K.D * C
	
	dqh = 0.5 * [0.0 -w'; w -S(w)] * q + K.a * (1.0 / (q' * q) - 1.0) * q
	dbh = -K.G * C
	
	return [dqh; dbh]
end
