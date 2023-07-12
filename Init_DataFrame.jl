function Init_Velocities_FW!(df::DataFrame, Δt::Float64)

    df[!, :v_x] = fill(0.0, nrow(df))
    df[!, :v_y] = fill(0.0, nrow(df))

    gdf = groupby(df, :ID)


    for df_i in gdf

        if length(df_i.x) > 2

        df_i.v_x[1:end-1] = (df_i.x[2:end].-df_i.x[1:end-1])./((df_i.Frame[2:end].-df_i.Frame[1:end-1])*Δt)
        df_i.v_y[1:end-1] = (df_i.y[2:end].-df_i.y[1:end-1])./((df_i.Frame[2:end].-df_i.Frame[1:end-1])*Δt)

        df_i.v_x[end] = df_i.v_x[end-1]
        df_i.v_y[end] = df_i.v_y[end-1]

        end

    end


end

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

function Init_Velocities!(df::DataFrame, k::Int, Δt::Float64)

    n = nrow(df)
    df[!, :v_x] = fill(0.0, n)
    df[!, :v_y] = fill(0.0, n)

    gdf = groupby(df, :ID)

    for df_i in gdf

        df_i.v_x[k+1:end-k] = (df_i.x[2*k+1:end].-df_i.x[1:end-2*k])./(2*k*Δt)
        df_i.v_y[k+1:end-k] = (df_i.y[2*k+1:end].-df_i.y[1:end-2*k])./(2*k*Δt)

        if k == 1

            df_i.v_x[1] = df_i.v_x[2]
            df_i.v_x[end] = df_i.v_x[end-1]

            df_i.v_y[1] = df_i.v_y[2]
            df_i.v_y[end] = df_i.v_y[end-1]

        end

    end

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

function Init_TTC!(df::DataFrame,)
    #Init columns
    df[!, :TTC] = fill(0.0, nrow(df))

    #group by frames
    gdf = groupby(df, :Frame)

    #Init TTC and TG
    for df_f in gdf

        TTCs = fill(0.0, nrow(df_f))
        j = 1

        for id in unique(df_f.ID)

            df_f_i = df_f[df_f.ID .== id, :]
            df_f_ = df_f[df_f.ID .!= id, :]

            TTCs[j] = Min_TTC(df_f_i, df_f_)

            j += 1

        end

        df_f.TTC = TTCs

    end

end

function Init_Min_R!(df::DataFrame)
    #Init columns
    df[!, :r] = fill(0.0, nrow(df))

    #group by frames
    gdf = groupby(df, :Frame)

    #Init TTC and TG
    for df_f in gdf

        rs = fill(0.0, nrow(df_f))
        j = 1

        for id in unique(df_f.ID)

            df_f_i = df_f[df_f.ID .== id, :]
            df_f_ = df_f[df_f.ID .!= id, :]

            rs[j] = Min_R(df_f_i, df_f_)
            j += 1

        end

        df_f.r = rs

    end

end

function Init_TG_TTC!(df::DataFrame, l)
    #Init columns
    df[!, :TTC] = fill(0.0, nrow(df))
    df[!, :TG] = fill(0.0, nrow(df))

    #group by frames
    gdf = groupby(df, :Frame)

    #Init TTC and TG
    for df_f in gdf

        TGs, TTCs = fill(0.0, nrow(df_f)), fill(0.0, nrow(df_f))
        j = 1

        for id in unique(df_f.ID)

            df_f_i = df_f[df_f.ID .== id, :]
            df_f_ = df_f[df_f.ID .!= id, :]

            TGs[j] = Min_TG_vC(df_f_i, df_f_, l)
            TTCs[j] = Min_TTC(df_f_i, df_f_, l)
            #R_soc = 0.2
            #TTCs[j] = Min_TTC_VariableRadius(df_f_i, df_f_, R_soc)
            j += 1


        end

        df_f.TG = TGs
        df_f.TTC = TTCs

    end

end

function Init_Min_R_1D!(df::DataFrame)
    #Init columns
    df[!, :r] = fill(0.0, nrow(df))

    #group by frames
    gdf = groupby(df, :Frame)

    #Init TTC and TG
    for df_f in gdf

        rs = fill(0.0, nrow(df_f))
        j = 1

        for id in unique(df_f.ID)

            df_f_i = df_f[df_f.ID .== id, :]
            df_f_ = df_f[df_f.ID .!= id, :]

            rs[j] = Min_R_1D(df_f_i, df_f_)
            j += 1

        end

        df_f.r = rs

    end

end

function Init_TG_TTC_1D!(df::DataFrame, l)
    #Init columns
    df[!, :TTC] = fill(0.0, nrow(df))
    df[!, :TG] = fill(0.0, nrow(df))

    #group by frames
    gdf = groupby(df, :Frame)

    #Init TTC and TG
    for df_f in gdf

        TGs, TTCs = fill(0.0, nrow(df_f)), fill(0.0, nrow(df_f))
        j = 1

        for id in unique(df_f.ID)

            df_f_i = df_f[df_f.ID .== id, :]
            df_f_ = df_f[df_f.ID .!= id, :]

            TGs[j] = Min_TG_1D_vC(df_f_i, df_f_, l)
            TTCs[j] = Min_TTC_1D(df_f_i, df_f_, l)

            #R_soc = 0.25
            #TTCs[j] = Min_TTC_1D_VariableRadius(df_f_i, df_f_, R_soc)
            j += 1

        end

        df_f.TG = TGs
        df_f.TTC = TTCs

    end

end

function Reset_Undef(x::Vector)

    for (i, x_) in enumerate(x)

        if x_ == 999.9

            x[i] = 9999999999999999999.9

        end

    end

    x

end

function Init_Min_R_F!(df::DataFrame, ϕ_)
    #Init columns
    df[!, :r_F] = fill(0.0, nrow(df))

    #group by frames
    gdf = groupby(df, :Frame)

    #Init TTC and TG
    for df_f in gdf

        rs = fill(0.0, nrow(df_f))
        j = 1

        for id in unique(df_f.ID)

            df_f_i = df_f[df_f.ID .== id, :]
            df_f_ = df_f[df_f.ID .!= id, :]

            rs[j] = Min_R_ϕ(df_f_i, df_f_, ϕ_)
            j += 1

        end

        df_f.r_F = rs

    end

end
