using LinearAlgebra

function S(v)
	[   0  -v[3]  v[2];
      v[3]    0  -v[1];
     -v[2]  v[1]    0]
end

function wahba(U, V)
	Y = svd(V * U')
	return Y.U * Diagonal([1.0, 1.0, det(Y.U) * det(Y.V)]) * Y.Vt
end

function q2R(q)
	r = 2 * q / (q' * q)
	return [1.0-q[3]*r[3]-q[4]*r[4] q[2]*r[3]-q[1]*r[4] q[2]*r[4]+q[1]*r[3];
			q[2]*r[3]+q[1]*r[4] 1.0-q[2]*r[2]-q[4]*r[4] q[3]*r[4]-q[1]*r[2];
			q[2]*r[4]-q[1]*r[3] q[3]*r[4]+q[1]*r[2] 1.0-q[2]*r[2]-q[3]*r[3]]
end

function quat_angle(q1, q2)
	dq = q2[1] * q1 + [0 q2[2:4]'; -q2[2:4] S(q2[2:4])] * q1
	return abs(2 * atan(norm(dq[2:4]) / dq[1]))
end

function R_angle(q, R)
	dR = R * q2R(q)
	return abs(atan(norm(dR - dR'), sqrt(2) * (tr(dR) - 1)))
end

function q_mul(q1, q2)
	return [q1[1] * q2[1] - q1[2:4]' * q2[2:4]; q1[1] * q2[2:4] + q2[1] * q1[2:4] + 
			cross(q1[2:4], q2[2:4])] / norm(q1) / norm(q2)
end

function q_div(q1, q2)
	return q_mul(q1, [q2[1]; -q2[2:4]])
end
