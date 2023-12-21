function Calc_ρ_Field(df_Grid, df, ξ)

    ρ_s = fill(0.0, nrow(df_Grid))

    #Filter Trajectories to those within the Grid
    x_min, x_max, y_min, y_max = minimum(df_Grid.x), maximum(df_Grid.x), minimum(df_Grid.y), maximum(df_Grid.y)
    filter!(row -> (x_min-2*ξ .≤ row.x .≤ x_max+2*ξ) .& (y_min-2*ξ .≤ row.y .≤ y_max+2*ξ), df)

    gdf = groupby(df, :Frame)

    for df_i ∈ gdf

        ρ_s = ρ_s .+ Calc_ρ_Field_SingleFrame(df_Grid, df_i, ξ)

    end

    ρ_s./length(gdf)

end

function Calc_ρ_Field_SingleFrame(df_G, df, ξ)

    ρ_s = fill(0.0, nrow(df_G))

    for (i, row) ∈ enumerate(eachrow(df_G))


        ρ_s[i] = ρ_s[i] + Calc_ρ_Field(df, row.x[1], row.y[1], ξ)

    end

    ρ_s

end

function Calc_ρ_Field(df, x, y, ξ)

   ρ, norm = 0.0, 0.0

    for row in eachrow(df)

        wheight = Gaussian_Wheight(d((row.x, row.y), (x, y)), ξ)

        ρ += wheight*row.ρ
        norm += wheight

    end

    if norm == 0.0
        0.0
    else
        ρ/norm
    end

end

function Init_Vornoi_Density!(df)

    n = nrow(df)
    df[!, :ρ] = fill(0.0, n)

    ε = 0.025
    gdf = groupby(df, :Frame)
    rect = Rectangle(Point2(minimum(df.x)-ε, minimum(df.y)-ε), Point2(maximum(df.x)+ε, maximum(df.y)+ε))

    for df_ in gdf

        points = [Point2(df_.x[i], df_.y[i]) for i in 1:nrow(df_)];
        tess = voronoicells(points, rect);
        df_.ρ = 1 ./voronoiarea(tess)

    end
end

TriangleHelp(p1, p2, p3) = (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])

function Point_In_Triangle(point, tri_1, tri_2, tri_3)

    d1 = TriangleHelp(point, tri_1, tri_2);
    d2 = TriangleHelp(point, tri_2, tri_3);
    d3 = TriangleHelp(point, tri_3, tri_1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos)

end

function Point_In_VoronoiCell(Cell::Vector{Point2{Float64}}, position_cell ,point)

    a = false
    i = 1

    while i < length(Cell) && a == false

        a = Point_In_Triangle(point, position_cell, Cell[i], Cell[i+1])
        i=i+1
    end

    if a == false

        a =  Point_In_Triangle(point, position_cell, Cell[length(Cell)], Cell[1])

    end


    return a

end

function VoronoiDensity(point, df::SubDataFrame, tess)

    ρ_p = 0.0
    i = 1
    n = nrow(df)

    while i <= n && ρ_p == 0.0

        if Point_In_VoronoiCell(tess.Cells[i], (df.x[i], df.y[i]), point) == true
            ρ_p = df.ρ[i]
        end

        i+=1

    end

    ρ_p

end

function Calc_DensityGrid(df_Grid::DataFrame, df::SubDataFrame)

    n = nrow(df_Grid)
    ρ_s = fill(0.0, n)
    ct = fill(0, n)

    ε = 0.025

    rect = Rectangle(Point2(minimum(df.x)-ε, minimum(df.y)-ε), Point2(maximum(df.x)+ε, maximum(df.y)+ε))
    points = [Point2(df.x[i], df.y[i]) for i in 1:nrow(df)];
    tess = voronoicells(points, rect);

    for i in 1:n
        ρ_s[i] = VoronoiDensity((df_Grid.x[i], df_Grid.y[i]), df, tess)

        if ρ_s[i] != 0.0
            ct[i] += 1
        end

    end

    ρ_s, ct

end

function Init_DensityGrid!(df_Grid::DataFrame, df::DataFrame)

    gdf = groupby(df, :Frame)
    cts = fill(0, nrow(df_Grid))


    for df_i in gdf
        ρ_s, cts_ = Calc_DensityGrid(df_Grid, df_i)
        df_Grid.ρ = df_Grid.ρ .+ ρ_s
        cts = cts .+ cts
    end

    for i in 1:nrow(df_Grid)
        if cts[i] != 0
            df_Grid.ρ[i] = df_Grid.ρ[i]/cts[i]
        end
    end

end

function Calc_DensityGrid(df_Grid::DataFrame, df::DataFrame)

    gdf = groupby(df, :Frame)
    ρ_s, cts = fill(0.0, nrow(df_Grid)), fill(0, nrow(df_Grid))

    for df_i in gdf
        ρ_s_, cts_ = Calc_DensityGrid(df_Grid, df_i)
        ρ_s = ρ_s .+ ρ_s_
        cts = cts .+ cts_
    end

    for i in 1:nrow(df_Grid)
        if cts[i] != 0
            ρ_s[i] = ρ_s[i]/cts[i]
        end
    end

    ρ_s

end

function Find_Cell(point, df, tess)

    for (i, cell) in enumerate(tess.Cells)

        if Point_In_VoronoiCell(cell, (df.x[i], df.y[i]), point) == true
            return i
        end

    end
    return false
end

function Find_All_Overlapping_Cells(df_Grid_INT, df, tess)

    index_ = fill(0, nrow(df_Grid_INT))

    for i in 1:nrow(df_Grid_INT)

        point = [df_Grid_INT.x[i], df_Grid_INT.y[i]]
        index_[i] = Find_Cell(point, df, tess)

    end

    unique_index = unique(index_)
    count_index = [length(filter(y-> y == i, index_)) for i in unique_index]

    unique_index, count_index

end

Recalc_Density(ρ, N, Δx, Δy) = 1/(1/ρ - N*(Δx*Δy))

function Recalc_Density!(df, N, Δx, Δy)
    df.ρ = 1 ./(1 ./df.ρ .- N*(Δx*Δy))
end


function Init_Voronoi_Density_Intruder!(df::SubDataFrame, rect, df_Grid_INT, Δx, Δy)

    ϵ = 0.025
    points = [Point2(df.x[j], df.y[j]) for j in 1:nrow(df)]
    tess = voronoicells(points, rect)
    df.ρ = 1 ./voronoiarea(tess)

    unique_index, count_index = Find_All_Overlapping_Cells(df_Grid_INT, df, tess)

    for (i, index_) in enumerate(unique_index)

        df.ρ[index_] = Recalc_Density(df.ρ[index_], count_index[i], Δx, Δy)

    end

end

function Init_Voronoi_Density_Intruder!(df::DataFrame, df_Grid_INT, Δx, Δy)

    df[!, :ρ] = fill(0.0, nrow(df))

    ε = 0.025
    rect = Rectangle(Point2(minimum(df.x)-ε, minimum(df.y)-ε), Point2(maximum(df.x)+ε, maximum(df.y)+ε))

    gdf = groupby(df, :Frame)

    for df_ in gdf

        Init_Voronoi_Density_Intruder!(df_, rect, df_Grid_INT, Δx, Δy)

    end

end;
