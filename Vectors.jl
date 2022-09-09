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

Base.abs(a, b) = sqrt(a^2+b^2)
