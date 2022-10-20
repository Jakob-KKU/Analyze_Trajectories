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

function Calc_All_Distances_SameFrames(df, Δfr = 1)

    gdf = groupby(df, :Frame)

    d_ = Vector{Float64}()

    for df_f in gdf[1:Δfr:end]

        append!(d_, Calc_All_Distances_singleFrame(df_f))

    end

    d_

end

function Calc_All_Distances_singleFrame(df)

    gdf = groupby(df, :ID)
    ids = length(unique(df.ID))

    d_ = fill(0.0, Int(ids*(ids-1)/2))
    ind = 1


    for (i, df_i) in enumerate(gdf)

        for df_j in gdf[i+1:end]

            d_[ind] = d(df_i.x[1], df_i.y[1] , df_j.x[1], df_j.y[1])
            ind += 1

        end
    end

    d_
end

function Calc_Distances_NotSameFrames(df, N=10000)

    d_ = fill(0.0, N)

    gdf = groupby(df, :ID)

    ids = 1:length(unique(df.ID))

    for i in 1:N

        id1, id2 = rand(ids), rand(ids)

        while length(intersect(gdf[id1].Frame, gdf[id2].Frame)) > 0
           id2 = rand(ids)
        end

        f1, f2 = rand(gdf[id1].Frame), rand(gdf[id2].Frame)

        x1, y1 = gdf[id1][gdf[id1].Frame .== f1, :].x[1], gdf[id1][gdf[id1].Frame .== f1, :].y[1]
        x2, y2 = gdf[id2][gdf[id2].Frame .== f2, :].x[1], gdf[id2][gdf[id2].Frame .== f2, :].y[1]

        d_[i] = d(x1, y1, x2, y2)

    end

    d_

end
