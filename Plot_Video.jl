function Plot_Boundaries!(x, y, size, color)

    GR.setmarkertype(GR.MARKERTYPE_ASTERISK)
    GR.setmarkersize(size)
    GR.setmarkercolorind(color)
    GR.polymarker(x, y)

end

function Boundaries_Rectangle_Positions(system_size)

    x = vcat(zeros(length(LinRange(0:0.1:system_size[2]))), LinRange(0:0.1:system_size[1]), fill(system_size[1], length(LinRange(0:0.1:system_size[2]))), LinRange(0:0.1:system_size[1]))
    y = vcat(LinRange(0:0.1:system_size[2]) , fill(system_size[2], length(LinRange(0:0.1:system_size[1]))), LinRange(0:0.1:system_size[2]), zeros(length(LinRange(0:0.1:system_size[1]))))

    return x, y
end

function Plot_Agents!(df::DataFrame, size, color)

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size)
    GR.setmarkercolorind(color)
    GR.polymarker(df.x, df.y)

end

function Plot_Agents!(df::DataFrame, IDs, size, color)

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size)
    GR.setmarkercolorind(color)

    for id in IDs
       GR.polymarker(df[df.ID .== id, :].x, df[df.ID .== id, :].y)
    end

end

function Plot_Agents!(df, size, color, size_2, color_2, dist)

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size)
    GR.setmarkercolorind(color)
    GR.polymarker(df.x, df.y)

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size_2)
    GR.setmarkercolorind(color_2)

    e_x = df.v_x ./ sqrt.(df.v_x .^2 .+ df.v_y .^2)
    e_y = df.v_y ./ sqrt.(df.v_x .^2 .+ df.v_y .^2)

    GR.polymarker(df.x .+ dist .* e_x, df.y .+ dist .* e_y)

end

function Plot_Agents!(df, IDs, size, color, size_2, color_2, dist)

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size)
    GR.setmarkercolorind(color)


    for id in IDs
       GR.polymarker(df[df.ID .== id, :].x, df[df.ID .== id, :].y)
    end


    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size_2)
    GR.setmarkercolorind(color_2)

    e_x = df.v_x ./ sqrt.(df.v_x .^2 .+ df.v_y .^2)
    e_y = df.v_y ./ sqrt.(df.v_x .^2 .+ df.v_y .^2)

    for id in IDs

        e_x = df[df.ID .== id, :].v_x ./ sqrt.(df[df.ID .== id, :].v_x .^2 .+ df[df.ID .== id, :].v_y .^2)
        e_y = df[df.ID .== id, :].v_y ./ sqrt.(df[df.ID .== id, :].v_x .^2 .+ df[df.ID .== id, :].v_y .^2)


        GR.polymarker(df[df.ID .== id, :].x .+ dist .* e_x, df[df.ID .== id, :].y .+ dist .* e_y)
    end

end

function Plot_Agents_TTC!(df::SubDataFrame, size)



    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size)

    gdf = groupby(df, :ID)


    for df_i in gdf

        red = max(min(1.0, 1.0 - df_i.TTC[1]/10), 0.0)
        GR.setcolorrep(1, red, 0.1, 0.1)
        GR.setmarkercolorind(1)
        GR.polymarker(df_i.x, df_i.y)

    end

end

function Plot_Agents_TTC!(df, size, size_2, color_2, dist)

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size)

    gdf = groupby(df, :ID)

    for df_i in gdf

        red = max(min(1.0, 1.0 - df_i.TTC[1]/10), 0.0)
        GR.setcolorrep(1, red, 0.1, 0.1)
        GR.setmarkercolorind(1)
        GR.polymarker(df_i.x, df_i.y)

    end

    GR.setmarkertype(GR.MARKERTYPE_SOLID_CIRCLE)
    GR.setmarkersize(size_2)
    GR.setmarkercolorind(color_2)

    e_x = df.v_x ./ sqrt.(df.v_x .^2 .+ df.v_y .^2)
    e_y = df.v_y ./ sqrt.(df.v_x .^2 .+ df.v_y .^2)

    GR.polymarker(df.x .+ dist .* e_x, df.y .+ dist .* e_y)

end
