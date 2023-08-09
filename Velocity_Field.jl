#Initialize velocity field of trajectories df in the grid df_G with Gaussian Kernels of width ξ
function Init_Velocity_Field!(df_Grid::DataFrame, df::DataFrame, ξ)

    #Init new rows in Grid
    df_Grid[!, :v_x] = fill(0.0, nrow(df_Grid))
    df_Grid[!, :v_y] = fill(0.0, nrow(df_Grid))

    #Filter Trajectories to those within the Grid
    x_min, x_max, y_min, y_max = minimum(df_Grid.x), maximum(df_Grid.x), minimum(df_Grid.y), maximum(df_Grid.y)
    filter!(row -> (x_min-2*ξ .≤ row.x .≤ x_max+2*ξ) .& (y_min-2*ξ .≤ row.y .≤ y_max+2*ξ), df)

    #Iterate over all Frames in df
    gdf = groupby(df, :Frame)

    for df_frame in gdf

        Init_Velocity_Field_SingleFrame!(df_Grid, df_frame, ξ)

    end

    #Normalize v_x and v_y
    df_Grid.v_x = df_Grid.v_x ./ length(gdf)
    df_Grid.v_y = df_Grid.v_y ./ length(gdf)

end

#Initialize velocity field of single Frame
function Init_Velocity_Field_SingleFrame!(df_Grid::DataFrame, df::SubDataFrame, ξ)

    for row in eachrow(df_Grid)

        v_Grid_x, v_Grid_y = Calc_Velocity_Field(df, row.x, row.y, ξ)

        row.v_x += v_Grid_x
        row.v_y += v_Grid_y

    end

end

#Calculate velocity field of trajectories df in the grid df_G with Gaussian Kernels of width ξ
function Calc_Velocity_Field(df_Grid::DataFrame, df::DataFrame, ξ)

    v_x, v_y = fill(0.0, nrow(df_Grid)), fill(0.0, nrow(df_Grid))

    #Filter Trajectories to those within the Grid
    x_min, x_max, y_min, y_max = minimum(df_Grid.x), maximum(df_Grid.x), minimum(df_Grid.y), maximum(df_Grid.y)
    filter!(row -> (x_min-2*ξ .≤ row.x .≤ x_max+2*ξ) .& (y_min-2*ξ .≤ row.y .≤ y_max+2*ξ), df)

    #Iterate over all Frames in df
    gdf = groupby(df, :Frame)

    for df_frame in gdf

        v_x_f, v_y_f = Calc_Velocity_Field_SingleFrame(df_Grid, df_frame, ξ)

        v_x = v_x .+ v_x_f
        v_y = v_y .+ v_y_f

    end

    #Normalize and return v_x and v_y
    v_x ./ length(gdf), v_y ./ length(gdf)

end

#Initialize velocity field of single Frame
function Calc_Velocity_Field_SingleFrame(df_Grid::DataFrame, df::SubDataFrame, ξ)

    v_x, v_y = fill(0.0, nrow(df_Grid)), fill(0.0, nrow(df_Grid))


    for (i, row) in enumerate(eachrow(df_Grid))

        v_Grid_x, v_Grid_y = Calc_Velocity_Field(df, row.x, row.y, ξ)

        v_x[i] += v_Grid_x
        v_y[i] += v_Grid_y

    end

    v_x, v_y

end


function Calc_Velocity_Field(df::SubDataFrame, x, y, ξ)

    v_x, v_y, norm = 0.0, 0.0, 0.0

    for row in eachrow(df)

        wheight = Gaussian_Wheight(d((row.x, row.y), (x, y)), ξ)

        v_x += wheight*row.v_x
        v_y += wheight*row.v_y
        norm += wheight

    end

    if norm == 0.0
        0.0, 0.0
    else
        v_x/norm, v_y/norm
    end

end

Gaussian_Wheight(dist::Float64, ξ::Float64) = dist ≤ 2*ξ ? exp(-dist^2/(2*ξ^2)) : 0.0

function Round_Velocity(v_x, v_y, ϵ = 0.01)

    if abs((v_x, v_y)) ≤ ϵ

        0.0, 0.0

    else

        v_x, v_y

    end

end

function Round_Velocities(v_xs, v_ys, ϵ = 0.01)

    v_xs_round, v_ys_round = fill(0.0, length(v_xs)), fill(0.0, length(v_xs))

    for i ∈ 1:length(v_xs)

        v_xs_round[i], v_ys_round[i] = Round_Velocity(v_xs[i], v_ys[i], ϵ)

    end

    v_xs_round, v_ys_round

end
