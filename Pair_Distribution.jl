function Calc_Dist(x, dx=0.1)

    x_min = minimum(x)
    x_ = collect(x_min:dx:maximum(x)+dx)
    p_x = fill(0.0, length(x_))

    for x_i in x

        p_x[Int(floor((x_i-x_min)/dx))+1] += 1

    end

    p_x ./length(x), x_

end

function Calc_Dist(x, x_bins::Vector, dx)

    x_min = minimum(x_bins)
    x_max = maximum(x_bins)

    ct = 0


    p_x = fill(0.0, length(x_bins))

    for x_i in x

        if x_min < x_i < x_max

            p_x[Int(floor((x_i-x_min)/dx))+1] += 1
            ct += 1

        end

    end

    p_x ./ct

end

function Num_of_Mutual_Dist(gdf, Δf)

    len = 0

    for df_f in gdf[1:Δf:end]

        ids = length(unique(df_f.ID))

        len+= Int(ids*(ids-1)/2)

    end

    len

end

function Add_Obs!(obs, df_i, df_j, fi, fj, line)

    obs[line, 1:4] .= df_i.x[fi], df_i.y[fi], df_i.v_x[fi], df_i.v_y[fi]
    obs[line, 5:8] .= df_j.x[fj], df_j.y[fj], df_j.v_x[fj], df_j.v_y[fj]
    obs[line, 9] = d(df_i.x[fi], df_i.y[fi] , df_j.x[fj], df_j.y[fj])
    obs[line, 10] = TTC((df_i.x[fi], df_i.y[fi]), (df_j.x[fj], df_j.y[fj]),
        (df_i.v_x[fi], df_i.v_y[fi]), (df_j.v_x[fj], df_j.v_y[fj]))
    obs[line, 11] = Rate_Of_Approach((df_i.x[fi], df_i.y[fi]), (df_j.x[fj], df_j.y[fj]),
        (df_i.v_x[fi], df_i.v_y[fi]), (df_j.v_x[fj], df_j.v_y[fj]))

end

function Calc_DF_Interaction(df, Δf)

    gdf = groupby(df, :Frame)
    N_val = Num_of_Mutual_Dist(gdf, Δf)
    obs = Matrix{Float64}(undef, N_val, 11)
    line = 1

    for df_f in gdf[1:Δf:end]

        gdf_2 = groupby(df_f, :ID)

        for (i, df_i) in enumerate(gdf_2)

            for df_j in gdf_2[i+1:end]

                Add_Obs!(obs, df_i, df_j, 1, 1, line)
                line += 1

            end
        end

    end

    DataFrame(x1 = obs[:, 1], y1 = obs[:, 2], v_x1 = obs[:, 3], v_y1 = obs[:, 4]
                 ,x2 = obs[:, 5], y2 = obs[:, 6], v_x2 = obs[:, 7], v_y2 = obs[:, 8]
                 , r = obs[:, 9], ttc = obs[:, 10], roa = obs[:, 11])

end

function Calc_DF_Independent(df, Δf, N_val)

    gdf_id = groupby(df, :ID)
    obs_ind = Matrix{Float64}(undef, N_val, 11)
    line = 1

    while line < N_val

        for df_i in gdf_id

            df_j = gdf_id[rand(1:length(gdf_id))]

            if Intersection(df_i.Frame, df_j.Frame) == false

                fr_i, fr_j = rand(1:length(df_i.x)), rand(1:length(df_j.x))

                Add_Obs!(obs_ind, df_i, df_j, fr_i, fr_j, line)

                line = min(line + 1, N_val)

            end

        end
    end


    DataFrame(x1 = obs_ind[:, 1], y1 = obs_ind[:, 2], v_x1 = obs_ind[:, 3], v_y1 = obs_ind[:, 4]
                 ,x2 = obs_ind[:, 5], y2 = obs_ind[:, 6], v_x2 = obs_ind[:, 7], v_y2 = obs_ind[:, 8]
                 , r = obs_ind[:, 9], ttc = obs_ind[:, 10], roa = obs_ind[:, 11])
end

function Add_Obs_1d!(obs, df_i, df_j, fi, fj, line)

    obs[line, 1:4] .= df_i.x[fi], df_i.y[fi], df_i.v_x[fi], df_i.v_y[fi]
    obs[line, 5:8] .= df_j.x[fj], df_j.y[fj], df_j.v_x[fj], df_j.v_y[fj]
    obs[line, 9] = d(df_i.x[fi], df_j.x[fj])
    obs[line, 10] = TTC(df_i.x[fi], df_j.x[fj], df_i.v_x[fi], df_j.v_x[fj])
    obs[line, 11] = Rate_Of_Approach(df_i.x[fi], df_j.x[fj], df_i.v_x[fi], df_j.v_x[fj])
    obs[line, 12] = min(TimeGap(df_i.x[fi], df_j.x[fj], df_i.v_x[fi]), TimeGap(df_j.x[fj], df_i.x[fi], df_j.v_x[fj]))


end

function Calc_DF_Interaction_1d(df, Δf)

    gdf = groupby(df, :Frame)
    N_val = Num_of_Mutual_Dist(gdf, Δf)
    obs = Matrix{Float64}(undef, N_val, 12)
    line = 1

    for df_f in gdf[1:Δf:end]

        gdf_2 = groupby(df_f, :ID)

        for (i, df_i) in enumerate(gdf_2)

            for df_j in gdf_2[i+1:end]

                Add_Obs_1d!(obs, df_i, df_j, 1, 1, line)
                line += 1

            end
        end

    end

    DataFrame(x1 = obs[:, 1], y1 = obs[:, 2], v_x1 = obs[:, 3], v_y1 = obs[:, 4]
                 ,x2 = obs[:, 5], y2 = obs[:, 6], v_x2 = obs[:, 7], v_y2 = obs[:, 8]
                 , r = obs[:, 9], ttc = obs[:, 10], roa = obs[:, 11], tg = obs[:, 12])

end

function Calc_DF_Independent_1d(df, Δf, N_val)

    gdf_id = groupby(df, :ID)
    obs_ind = Matrix{Float64}(undef, N_val, 12)
    line = 1

    while line < N_val

        for df_i in gdf_id

            df_j = gdf_id[rand(1:length(gdf_id))]

            if Intersection(df_i.Frame, df_j.Frame) == false

                fr_i, fr_j = rand(1:length(df_i.x)), rand(1:length(df_j.x))

                Add_Obs_1d!(obs_ind, df_i, df_j, fr_i, fr_j, line)

                line = min(line + 1, N_val)

            end

        end
    end


    DataFrame(x1 = obs_ind[:, 1], y1 = obs_ind[:, 2], v_x1 = obs_ind[:, 3], v_y1 = obs_ind[:, 4]
                 ,x2 = obs_ind[:, 5], y2 = obs_ind[:, 6], v_x2 = obs_ind[:, 7], v_y2 = obs_ind[:, 8]
                 , r = obs_ind[:, 9], ttc = obs_ind[:, 10], roa = obs_ind[:, 11],  tg = obs_ind[:, 12])
end
