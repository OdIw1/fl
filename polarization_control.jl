function R(a)
	JonesMatrix([cos(a)  -sin(a); sin(a) -cos(a)])
end

function W2(a)
    half_wave = [-1.im 0; 0 1.im]
	JonesMatrix(R(a) * half_wave * R(-a))
end

function W4(a)
    quarter_wave = 1. / sqrt(2) * [(1. - 1.im)  0; 0 (1. + 1.im)]
	JonesMatrix(R(a) * quarter_wave * R(-a))
end

function apply_Jones_matrix!(M::JonesMatrix, uX, uY)
    a = M.m
    for i = 1:length(uX)
        uX_ = uX[i]                 
        uY_ = uY[i]
        uX[i] = a[1,1] * uX_ + a[1,2] * uY_
        uY[i] = a[2,1] * uX_ + a[2,2] * uY_
    end
end

apply_Jones_matrix!(M::JonesMatrix, p::Pulse) = apply_Jones_matrix!(M, p.uX, p.uY)