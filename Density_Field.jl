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

end;
