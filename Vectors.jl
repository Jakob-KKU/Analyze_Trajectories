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
d(a::Float64, b::Float64) = sqrt((a-b)^2)

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

function Rate_Of_Approach(x_a::NTuple{2, Float64}, x_b::NTuple{2, Float64}, v_a::NTuple{2, Float64}, v_b::NTuple{2, Float64})

    -1 .*(v_a .- v_b)⋅(x_a.-x_b)/d(x_a, x_b)

end

function Rate_Of_Approach(x_a::Float64, x_b::Float64, v_a::Float64, v_b::Float64)

    -1*(v_a - v_b)*(x_a-x_b)/d(x_a, x_b)

end

Intersection(fr1, fr2) = fr1[1] > fr2[end] || fr2[1] > fr1[end] ? false : true

ϕ(a::NTuple{2, Float64}, b::NTuple{2, Float64}) = acos(clamp((a⋅b)/(abs(a)*abs(b)),-1.0, 1.0))

function Min_R(df_i, df_f)

    x_i = (df_i.x[1], df_i.y[1])
    r_min = 99999.9

    for row in eachrow(df_f)
        x_j = (row.x, row.y)
        r_min = min(r_min, d(x_i, x_j))
    end

    r_min

end

function Min_R_1D(df_i, df_f)

    x_i = df_i.x[1]
    r_min = 99999.9

    for row in eachrow(df_f)
        x_j = row.x
        r_min = min(r_min, d(x_i, x_j))
    end

    r_min

end

IN(r::Vector, r_soc, dia) = [IN(r_i, r_soc, dia) for r_i in r]

IN(r::Float64, r_soc, dia) = (r_soc - dia)/(r-dia)

AV(ttc::Vector, T) = T./ttc

AV(ttc::Float64, T) = T/ttc

function ρ_Global(df)

    A = (maximum(df.x)-minimum(df.x))*(maximum(df.y)-minimum(df.y))

    rho = 0.0

    gdf = groupby(df, :Frame)

    for df_ in gdf

    rho += length(df_.ID)

    end

    rho = rho/(A*length(gdf))

end

function Min_R_ϕ(df_i, df_f, ϕ_)

    x_i = (df_i.x[1], df_i.y[1])
    v_i = (df_i.v_x[1], df_i.v_y[1])
    r_min = 99999.9

    for row in eachrow(df_f)

        x_j = (row.x, row.y)

        if ϕ(v_i, x_j .- x_i) > ϕ_
            r_min = min(r_min, d(x_i, x_j))
        end
    end

    r_min

end
