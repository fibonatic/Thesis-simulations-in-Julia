function LTV_minimal(x, R, tau, K)
	M = reshape(x[K.xM2M],9,9)
	F = K.W * kron(Matrix(I,3,3), R) * K.G * R'
	cor = M * [K.Winv; zeros(3,6)] * (K.W * R[:] - x[1:6])
	A = zeros(9,9); A[1:6,7:9] = F
	cwc = zeros(9,9); cwc[1:6,1:6] = K.Winv

	drh = F * x[7:9] + cor[1:6]
	dLh = R * tau + cor[7:9]
	dM = A * M + M * A' - M * cwc * M + K.V + K.d * M

	return [drh; dLh; dM[K.M2xM]]
end

function LTV_biased(x, R, z, tau, K)
	M = reshape(x[K.xM2M],12,12)
	I3 = Matrix(I,3,3)
	C = [Matrix(I,6,12); zeros(3,6) K.J \ R' I3]
	F = K.W * kron(I3, R) * K.G * R'
	cor = M * C' * K.Winv * ([K.W * R[:]; z] - C * x[1:12])
	A = zeros(12,12); A[1:6,7:9] = F

	drh = F * x[7:9] + cor[1:6]
	dLh = R * tau + cor[7:9]
	dbh = cor[10:12]
	dM = A * M + M * A' - M * C' * K.Winv * C * M + K.V + K.d * M

	return [drh; dLh; dbh; dM[K.M2xM]]
end

function LTV_kinematic(x, R, z, K)
	M = reshape(x[K.xM2M],9,9)
	F = K.W * kron(Matrix(I,3,3), R) * K.G
	cor = M * [K.Winv * (K.W * R[:] - x[1:6]); zeros(3)]
	A = zeros(9,9); A[1:6,7:9] = -F
	cwc = zeros(9,9); cwc[1:6,1:6] = K.Winv

	drh = F * (z - x[7:9]) + cor[1:6]
	dbh = cor[7:9]
	dM = A * M + M * A' - M * cwc * M + K.V + K.d * M

	return [drh; dbh; dM[K.M2xM]];
end
