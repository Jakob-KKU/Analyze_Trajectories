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

function Calc_PD(x, x_scrambled ,x_bins::Vector, dx)

    x_min = minimum(x_bins)
    x_max = maximum(x_bins)


    pd = fill(0.0, length(x_bins))
    p_x = fill(0.0, length(x_bins))
    p_xs = fill(0.0, length(x_bins))

    ct = 0
    ct_s = 0

    for x_i in x

        if x_min < x_i < x_max+dx

            p_x[Int(floor((x_i-x_min)/dx))+1] += 1
            ct += 1

        end

    end

    for x_i in x_scrambled

        if x_min < x_i < x_max+dx

            p_xs[Int(floor((x_i-x_min)/dx))+1] += 1
            ct_s += 1

        end

    end

    #(p_x) ./ (p_xs)
    #(p_x./ct) ./ (p_xs./ct_s)
    (p_x./length(x)) ./ (p_xs./length(x_scrambled))


end

function Num_of_Mutual_Dist(gdf, Δf)

    len = 0

    for df_f in gdf[1:Δf:end]

        ids = length(unique(df_f.ID))

        len+= Int(ids*(ids-1)/2)

    end

    len

end

function Calc_DF_PairDist(df, Δf)

    gdf = groupby(df, :Frame)
    N_val = Num_of_Mutual_Dist(gdf, Δf)
    obs = Matrix{Float64}(undef, N_val, 13)
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
                 , x2 = obs[:, 5], y2 = obs[:, 6], v_x2 = obs[:, 7], v_y2 = obs[:, 8]
                 , r = obs[:, 9], ttc = obs[:, 10], roa = obs[:, 11], phi_1 = obs[:, 12], phi_2 = obs[:, 13])

end

function Calc_DF_PairDist_PairsOnly(df, Δf, f_min, d_mean, d_max)


    pair_ids = Calc_Pair_IDs(df, f_min, d_mean, d_max)

    gdf = groupby(df, :Frame)
    N_val = Num_of_Mutual_Dist(gdf, Δf)
    obs = Matrix{Float64}(undef, N_val, 13)
    line = 1

    for df_f in gdf[1:Δf:end]

        gdf_2 = groupby(df_f, :ID)

        for (i, df_i) in enumerate(gdf_2)

            if df_i.ID[1] ∈ pair_ids

                for df_j in gdf_2[i+1:end]

                    Add_Obs!(obs, df_i, df_j, 1, 1, line)
                    line += 1

                end
            end
        end

    end

    obs = obs[1:line, :]

    DataFrame(x1 = obs[:, 1], y1 = obs[:, 2], v_x1 = obs[:, 3], v_y1 = obs[:, 4]
                 , x2 = obs[:, 5], y2 = obs[:, 6], v_x2 = obs[:, 7], v_y2 = obs[:, 8]
                 , r = obs[:, 9], ttc = obs[:, 10], roa = obs[:, 11], phi_1 = obs[:, 12], phi_2 = obs[:, 13])

end

function Calc_DF_PairDist_NoPairs(df, Δf, f_min, d_mean, d_max)


    pair_ids = Calc_Pair_IDs(df, f_min, d_mean, d_max)

    gdf = groupby(df, :Frame)
    N_val = Num_of_Mutual_Dist(gdf, Δf)
    obs = Matrix{Float64}(undef, N_val, 13)
    line = 1

    for df_f in gdf[1:Δf:end]

        gdf_2 = groupby(df_f, :ID)

        for (i, df_i) in enumerate(gdf_2)

            if df_i.ID[1] ∉ pair_ids

                for df_j in gdf_2[i+1:end]

                    Add_Obs!(obs, df_i, df_j, 1, 1, line)
                    line += 1

                end
            end
        end

    end

    obs = obs[1:line, :]

    DataFrame(x1 = obs[:, 1], y1 = obs[:, 2], v_x1 = obs[:, 3], v_y1 = obs[:, 4]
                 , x2 = obs[:, 5], y2 = obs[:, 6], v_x2 = obs[:, 7], v_y2 = obs[:, 8]
                 , r = obs[:, 9], ttc = obs[:, 10], roa = obs[:, 11], phi_1 = obs[:, 12], phi_2 = obs[:, 13])

end

function Calc_DF_PairDist_1d(df, Δf)

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

function Add_Obs!(obs, df_i, df_j, fi, fj, line)

    x_i = (df_i.x[fi], df_i.y[fi])
    x_j = (df_j.x[fj], df_j.y[fj])
    v_i = (df_i.v_x[fi], df_i.v_y[fi])
    v_j = (df_j.v_x[fj], df_j.v_y[fj])

    obs[line, 1:4] .= df_i.x[fi], df_i.y[fi], df_i.v_x[fi], df_i.v_y[fi]
    obs[line, 5:8] .= df_j.x[fj], df_j.y[fj], df_j.v_x[fj], df_j.v_y[fj]
    obs[line, 9] = d(df_i.x[fi], df_i.y[fi] , df_j.x[fj], df_j.y[fj])

    obs[line, 10] = TTC(x_i, x_j, v_i, v_j, 0.3, 0.3)
    obs[line, 11] = Rate_Of_Approach(x_i, x_j, v_i, v_j)
    obs[line, 12] = ϕ(v_i, x_j .- x_i)
    obs[line, 13] = ϕ(v_j, x_i .- x_j)

end

function Add_Obs_1d!(obs, df_i, df_j, fi, fj, line)

    obs[line, 1:4] .= df_i.x[fi], df_i.y[fi], df_i.v_x[fi], df_i.v_y[fi]
    obs[line, 5:8] .= df_j.x[fj], df_j.y[fj], df_j.v_x[fj], df_j.v_y[fj]
    obs[line, 9] = d(df_i.y[fi], df_j.y[fj])
    obs[line, 10] = TTC(df_i.x[fi], df_j.x[fj], df_i.v_x[fi], df_j.v_x[fj], 0.2, 0.2)
    obs[line, 11] = Rate_Of_Approach(df_i.x[fi], df_j.x[fj], df_i.v_x[fi], df_j.v_x[fj])
    obs[line, 12] = min(TimeGap(df_i.x[fi], df_j.x[fj], df_i.v_x[fi]), TimeGap(df_j.x[fj], df_i.x[fi], df_j.v_x[fj]))


end

function Scramble_Time!(df)

    gdf = groupby(df, :ID)

    for i in 1:length(gdf)

        #gdf[i].Frame = [gdf[i].Frame[x] for x in rand(1:length(gdf[i].Frame), length(gdf[i].Frame))]

        shuffle!(gdf[i].Frame)

    end

end
