function d(a::Vector, b::Vector, L)

    c = fill(0.0, length(a))

    for i in 1:length(a)

        dx = a[i] - b[i]

        if dx > L/2
            dx = dx - L
        elseif dx < -L/2
            dx = dx + L
        end

        c[i] = dx

    end
    c
end

function d(a::NTuple{2, Float64}, b::NTuple{2, Float64}, system_size::NTuple{2, Float64})
    dx = abs(a[1]-b[1])
    dy = abs(a[2]-b[2])

    if dx > system_size[1]/2
        dx = system_size[1] - dx
    end
    if dy > system_size[2]/2
        dy = system_size[2] - dy
    end

    return sqrt(dx^2 + dy^2)
end

d(a::NTuple{2, Float64}, b::NTuple{2, Float64}) = sqrt((a[1]-b[1])^2+(a[2]-b[2])^2)
e_(a::NTuple{2, Float64}, b::NTuple{2, Float64}) = (a.-b)./d(a, b)
e_v(v_a::NTuple{2, Float64}, v_b::NTuple{2, Float64}) = normalize(v_a .- v_b)
normalize(a::NTuple{2, Float64}) = a./abs(a)
l(l1::Float64, l2::Float64) = l1/2 + l2/2
⋅(u::NTuple{2, Float64}, v::NTuple{2, Float64}) = u[1]*v[1]+u[2]*v[2]

d(x1, y1, x2, y2) = sqrt((x1-x2)^2+(y1-y2)^2)
function e_(x1, y1, x2, y2)

    (x1-x2)/d(x1, y1, x2, y2), (y1-y2)/d(x1, y1, x2, y2)

end


⋅(u, v) = u[1]*v[1]+u[2]*v[2]


Base.abs(a, b) = sqrt(a^2+b^2)
Base.abs(a::NTuple{2, Float64}) = sqrt(a[1]^2+a[2]^2)
