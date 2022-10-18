function TimeGap(x_a::NTuple{2, Float64}, x_b::NTuple{2, Float64}, v_a::NTuple{2, Float64}, l_a = 0.3, l_b = 0.3)

    cos_α = e_(x_a, x_b)⋅normalize(v_a)
    A = (cos_α^2-1)*d(x_a,x_b)^2+l(l_a, l_b)^2

    if A < 0.0 || -cos_α*d(x_a,x_b)-sqrt(A) < 0.0
        999.9
    else
        (-cos_α*d(x_a,x_b)-sqrt(A))/abs(v_a)
    end
end

function TimeGap(x1, x2, v1, l_a = 0.3, l_b = 0.3)

    tg_ = ((x2-x1)-l(l_a, l_b))/v1

    if tg_ < 0
        999.9
    else
        tg_
    end

end

function Min_TG_1D(df_i, index, df_f)

    x_i = df_i.x[index]
    v_i =  df_i.v_x[index]

    TG_min = 999.9

    for row in eachrow(df_f)
        TG_min = min(TG_min, TimeGap(x_i, row.x,  v_i))
    end

    TG_min

end

function Min_TG(df_i, index, df_f)

    x_i = (df_i.x[index], df_i.y[index])
    v_i = (df_i.v_x[index], df_i.v_y[index])

    TG_min = 999.9

    for row in eachrow(df_f)
        TG_min = min(TG_min, TimeGap(x_i, (row.x, row.y),  v_i))
    end

    TG_min

end

function Min_TG_1D(df_i, df_f)

    x_i = df_i.x[1]
    v_i =  df_i.v_x[1]

    TG_min = 999.9

    for row in eachrow(df_f)
        TG_min = min(TG_min, TimeGap(x_i, row.x,  v_i))
    end

    TG_min

end

function Min_TG_1D_vC(df_i, df_f)

    x_i = df_i.x[1]
    v_i =  1.0

    TG_min = 999.9

    for row in eachrow(df_f)
        TG_min = min(TG_min, TimeGap(x_i, row.x,  v_i))
    end

    TG_min

end

function Min_TG(df_i, df_f)

    x_i = (df_i.x[1], df_i.y[1])
    v_i = (df_i.v_x[1], df_i.v_y[1])

    TG_min = 999.9

    for row in eachrow(df_f)
        TG_min = min(TG_min, TimeGap(x_i, (row.x, row.y),  v_i))
    end

    TG_min

end
