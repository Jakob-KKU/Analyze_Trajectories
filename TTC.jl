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

        x_j = row.x
        v_j = row.v_x

        ttc_min = min(ttc_min, TTC(x_i, x_j, v_i, v_j, l, l))
    end

    ttc_min

end

function Min_TTC(df_i, df_f, l)

    x_i = (df_i.x[1], df_i.y[1])
    v_i = (df_i.v_x[1], df_i.v_y[1])

    ttc_min = 99999.9

    for row in eachrow(df_f)

        x_j = (row.x, row.y)
        v_j = (row.v_x, row.v_y)

        ttc_min = min(ttc_min, TTC(x_i, x_j, v_i, v_j, l, l))

    end

    ttc_min

end

function Min_TTC_1D_VariableRadius(df_i, df_f, R_soc)

    x_i = df_i.x[1]
    v_i = df_i.v_x[1]

    ttc_min = 99999.9

    for row in eachrow(df_f)

        x_j = row.x
        v_j = row.v_x

        diameter = 2*min(R_soc, d(x_i, x_j)/4)

        ttc_min = min(ttc_min, TTC(x_i, x_j, v_i, v_j, diameter, diameter))
    end

    ttc_min

end

function Min_TTC_VariableRadius(df_i, df_f, R_soc)

    x_i = (df_i.x[1], df_i.y[1])
    v_i = (df_i.v_x[1], df_i.v_y[1])

    ttc_min = 99999.9

    for row in eachrow(df_f)

        x_j = (row.x, row.y)
        v_j = (row.v_x, row.v_y)

        diameter = 2*min(R_soc, d(x_i, x_j)/4)

        ttc_min = min(ttc_min, TTC(x_i, x_j, v_i, v_j, diameter, diameter))

    end

    ttc_min

end


function Min_TTC_1D_Perception(df_i, df_f, l)

    x_i = df_i.x[1]
    v_i = df_i.v_x[1]

    ttc_min = 99999.9

    for row in eachrow(df_f)

        x_j = row.x
        v_j = row.v_x

        if Perceivable_1D(x_i, x_j, df_f, l) == true
            ttc_min = min(ttc_min, TTC(x_i, x_j, v_i, v_j, l, l))
        end
    end

    ttc_min

end

function Min_TTC_Perception(df_i, df_f, l)

    x_i = (df_i.x[1], df_i.y[1])
    v_i = (df_i.v_x[1], df_i.v_y[1])

    ttc_min = 99999.9

    for row in eachrow(df_f)

        x_j = (row.x, row.y)
        v_j = (row.v_x, row.v_y)

        #can we actually perceive agent j?
        if Perceivable(x_i, x_j, df_f, l) == true

            ttc_min = min(ttc_min, TTC(x_i, x_j, v_i, v_j, l, l))

        end
    end

    ttc_min

end

function Perceivable(x_i, x_j, df_f, l)

    perceivable_ = true

    for row in eachrow(df_f)

        x_k = (row.x, row.y)

        if x_k != x_j && Intersect_Line_Circle(x_i, x_j, x_k, l) == true

            perceivable_ = false
            break

        end
    end

    perceivable_

end

function Perceivable_1D(x_i::Float64, x_j::Float64, df_f, l)

    perceivable_ = true

    for row in eachrow(df_f)

        x_k = row.x

        if x_k != x_j && Intersect_Line_Circle((x_i, 0.0), (x_j, 0.0), (x_k, 0.0), l) == true

            perceivable_ = false
            break

        end
    end

    perceivable_

end

function Intersect_Line_Circle(x_i, x_j, x_k, l)

    r_i = x_i[1]+im*x_i[2]
    r_j = x_j[1]+im*x_j[2]
    r_k = x_k[1]+im*x_k[2]

    z = (r_k - r_i )/(r_j - r_i)

    distance = 0.0 < real(z) < 1.0 ? abs(imag(z)*(r_j - r_i)) : min(abs(r_j - r_k), abs(r_i - r_k))

    if distance <= l/2
        return true
    else
        return false
    end
end
