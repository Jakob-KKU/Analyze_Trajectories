function Is_Pair(df1, df2, f_min, d_mean, d_max)

    if Intersection_Length(df1.Frame, df2.Frame) > f_min

        f1, f2 = Intersection_Interval(df1.Frame, df2.Frame)

        df_i = filter(row -> f1 < row.Frame < f2, df1)
        df_j = filter(row -> f1 < row.Frame < f2, df2)

        #if length(df_i.x) != length(df_j.x)
            #println("Back and Forth movement out of Measurement area")
            #false

        #if length(df_i.x) == 0 || length(df_j.x) == 0
        #    println(f1 ," and ", f2)
        #end

        #else

            d_ij = [d((df_i.x[i], df_i.y[i]), (df_j.x[i], df_j.y[i])) for i in 1:length(df_i.x)]

            if maximum(d_ij) < d_max && mean(d_ij) < d_mean
                true
            else
                false
            end

        #end

    else

        false

    end

end

function Calc_Pair_IDs(df, f_min, d_mean, d_max)

    gdf = groupby(df, :ID)
    pair_ids = Vector{Int}()


    for i in 1:length(gdf)

        for j in i+1:length(gdf)

            #println(i, " and ", j)

            if Is_Pair(gdf[i], gdf[j], f_min, d_mean, d_max) == true
                push!(pair_ids, gdf[j].ID[1])
                push!(pair_ids, gdf[i].ID[1])

            end

        end
    end

    pair_ids
end

function Calc_Neglect_Pair_IDs(df, f_min, d_mean, d_max)

    gdf = groupby(df, :ID)
    neglect_ids = Vector{Int}()


    for i in 1:length(gdf)

        for j in i+1:length(gdf)

            if Is_Pair(gdf[i], gdf[j], f_min, d_mean, d_max) == true
                push!(neglect_ids, j)
            end

        end
    end

    neglect_ids
end



########################
# does not work right
function Init_Pairs!(df, f_min, d_mean, d_max)

    gdf = groupby(df, :ID)
    df[!, :Pairs] =  fill(Vector{Int}(), nrow(df))

    for i in 1:length(gdf)

        for j in i+1:length(gdf)

            if Is_Pair(gdf[i], gdf[j], f_min, d_mean, d_max) == true

                #println((i, j))

                #Init_Pair!(gdf, i, j)
                append!(gdf[i].Pairs[1], j)
                append!(gdf[j].Pairs[1], i)

                #gdf[i].Pairs[:] .= j
                #gdf[j].Pairs[:] .= i
                #Init_Pair!(gdf, j, i)

            end

        end

    end

end

function Init_Pair!(gdf, i, j)

    if length(gdf[i].Pairs) > 1
        append!(gdf[i].Pairs, j)
    else
        gdf[i].Pairs = [j]
    end

    if length(gdf[j].Pairs) > 1
        append!(gdf[j].Pairs, i)
    else
        gdf[j].Pairs = [i]
    end

end
