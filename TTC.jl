function TTC(x_a::NTuple{2, Float64}, x_b::NTuple{2, Float64}, v_a::NTuple{2, Float64}, v_b::NTuple{2, Float64},
        l_a = 0.3, l_b = 0.3)

    cos_α = e_(x_a, x_b)⋅e_v(v_a, v_b)
    A = (cos_α^2-1)*d(x_a,x_b)^2+l(l_a, l_b)^2

    if A < 0.0 || -cos_α*d(x_a,x_b)-sqrt(A) < 0.0 || abs(v_a .- v_b) == 0.0
        99999.9
    else
        (-cos_α*d(x_a,x_b)-sqrt(A))/abs(v_a .- v_b)
    end
end

function TTC(x1, y1, x2, y2, vx1, vy1, vx2, vy2, l_a = 0.3, l_b = 0.3)

    cos_α = e_(x1, y1, x2, y2)⋅e_(vx1, vy1, vx2, vy2)
    A = (cos_α^2-1)*d(x1, y1, x2, y2)^2+l(l_a, l_b)^2

    if A < 0.0 || -cos_α*d(x1, y1, x2, y2)-sqrt(A) < 0.0 || d(vx1, vy1, vx2, vy2) < 0.0
        99999.9
    else
        (-cos_α*d(x1, y1, x2, y2)-sqrt(A))/d(vx1, vy1, vx2, vy2)
    end
end

function TTC(x1, x2, v1, v2, l_a = 0.3, l_b = 0.3)

    ttc_ = ((x2-x1)-l(l_a, l_b))/(v1-v2)

    if ttc_ > 0.0 && v2-v1 != 0.0

        ttc_

    else

        99999.9

    end

end

function Min_TTC_1D(df_i, index, df_f)

    x_i = df_i.x[index]
    v_i = df_i.v_x[index]

    ttc_min = 99999.9

    for row in eachrow(df_f)
        ttc_min = min(ttc_min, TTC(x_i, row.x, v_i, row.v_x))
    end

    ttc_min

end

function Min_TTC(df_i, index, df_f)

    x_i = (df_i.x[index], df_i.y[index])
    v_i = (df_i.v_x[index], df_i.v_y[index])

    ttc_min = 99999.9

    for row in eachrow(df_f)
        ttc_min = min(ttc_min, TTC(x_i, (row.x, row.y), v_i, (row.v_x, row.v_y)))
    end

    ttc_min

end


function Min_TTC_1D(df_i, df_f, l)

    x_i = df_i.x[1]
    v_i = df_i.v_x[1]

    ttc_min = 99999.9

    for row in eachrow(df_f)
        ttc_min = min(ttc_min, TTC(x_i, row.x, v_i, row.v_x, l, l))
    end

    ttc_min

end

function Min_TTC(df_i, df_f, l)

    x_i = (df_i.x[1], df_i.y[1])
    v_i = (df_i.v_x[1], df_i.v_y[1])

    ttc_min = 99999.9

    for row in eachrow(df_f)
        ttc_min = min(ttc_min, TTC(x_i, (row.x, row.y), v_i, (row.v_x, row.v_y), l, l))
    end

    ttc_min

end
