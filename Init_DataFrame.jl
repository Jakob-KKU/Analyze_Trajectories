function Init_Velocities!(df::DataFrame, k::Int, Δt::Float64, L)

    n = nrow(df)
    v_x = fill(0.0, n)
    v_y = fill(0.0, n)

    ids = unique(df[!, :ID])
    index = 0


    for i in ids

        df_ = df[(df.ID .== i), :]
        n_ = nrow(df_)

        v_x[index+k+1:index+n_-k] = d(df_[!, :x][2*k+1:end], df_[!, :x][1:end-2*k], L[1])./(2*k*Δt)
        v_y[index+k+1:index+n_-k] = d(df_[!, :y][2*k+1:end], df_[!, :y][1:end-2*k], L[2])./(2*k*Δt)

        index += n_

    end

    df[!, :v_x] = v_x
    df[!, :v_y] = v_y;

end

function Init_Global_Density!(df::DataFrame, L)

    n = nrow(df)
    index = 0


    frames = unique(df[!, :Frame])
    ρ_global = fill(0.0, n)


    for fr in frames

        df_ = df[(df.Frame .== fr), :]
        n_ = nrow(df_)

        ρ_global[index+1:index+n_] .= n_/(L[1]*L[2])
        index += n_

    end

    df[!, :ρ_global] = ρ_global

end

function Init_Vornoi_Area!(df, L)

    gdf = groupby(df, :Frame)
    rect = Rectangle(Point2(0.0, 0.0), Point2(L))


    for df_ in gdf

        points = [Point2(df_.x[i], df_.y[i]) for i in 1:nrow(df_)];
        tess = voronoicells(points, rect);
        df_[!, :A] = voronoiarea(tess)

    end

    gdf
end


############## probably incorrect !!!! ############
function Init_Accelerations!(df, k, Δt)

    n = nrow(df)
    a_x = fill(0.0, n)
    a_y = fill(0.0, n)

    ids = unique(df[!, :ID])
    index = 0

    for i in ids

        df_ = df[(df.ID .== i), :]
        n_ = nrow(df_)

        a_x[index+2*k+1:index+n_-2*k] = (df_[!, :v_x][2*k+2:end-k] .- df_[!, :v_x][k+1:end-2*k-1])./(2*k*Δt)
        a_y[index+2*k+1:index+n_-2*k] = (df_[!, :v_y][2*k+2:end-k] .- df_[!, :v_y][k+1:end-2*k-1])./(2*k*Δt)

        index += n_

    end

    df[!, :a_x] = a_x
    df[!, :a_y] = a_y;
end
