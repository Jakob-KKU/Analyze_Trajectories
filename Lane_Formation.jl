function Travel_Velocities(Files, path, x_min, x_max, Δt, N, k, L)

    #Calculate MeanSpeed BaseLine
    v = fill(0.0, length(Files), N)

    for (j, f) in enumerate(Files)

        data = CSV.File(string(path, f); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data);
        rename!(df, 1 => "ID", 2 => "Frame", 3 => "x", 4 => "y");

        #Initialize the velocities
        Init_Velocities!(df, k, Δt);

        #Filter to measurement Area
        filter!(row -> x_min < row.x < x_max, df)

        gdf = groupby(df, :ID)
        i=1

        for df_i in gdf

            v[j, i] = (x_max-x_min)/(nrow(df_i)*Δt)
            i+=1

        end

    end

    v
end


function Init_Lane_Group!(df::DataFrame)

    df[!, :LF_ID] .= 0
    gdf_ID = groupby(df, :ID)

    for df_i in gdf_ID

        if mean(df_i.v_x) > 0.0
            df_i.LF_ID .= 1
        else
            df_i.LF_ID .= -1
        end

    end
end

## ORDER PARAMETER LANEFORMATION

# Iterate over agents
function ϕ_LF(df::DataFrame, l, Δt)

    frames = sort!(unique(df.Frame))
    op_LF = fill(0.0, length(frames))


    for (i, f) in enumerate(frames)
        op_LF[i] = ϕ_LF_t(df[(df.Frame .== f), :], l)
    end

    op_LF
end

function ϕ_LF_t(df::DataFrame, r)

    op_LF = 0.0

    for row in eachrow(df)

        L_i = floor(row.y/r)

        N_same = nrow(df[(df.ID .!= row.ID) .& (floor.(df.y/r) .== L_i) .& (row.LF_ID .* df.LF_ID .== 1), :])
        N_diff = nrow(df[(df.ID .!= row.ID) .& (floor.(df.y/r) .== L_i) .& (row.LF_ID .* df.LF_ID .== -1), :])

       op_LF+=ϕ_LF(N_same, N_diff)

    end

    op_LF/N
end

ϕ_LF(N_same, N_diff) = (N_same + N_diff) == 0 ? 0.0 : (N_same - N_diff)^2/(N_same + N_diff)^2

#Iterate over Rows
function ϕ_LF(df::DataFrame, l, L)

    frames = sort!(unique(df.Frame))
    op_LF = fill(0.0, length(frames))


    for (i, f) in enumerate(frames)
        op_LF[i] = ϕ_LF_t(df[(df.Frame .== f), :], l, L)
    end

    op_LF
end

function ϕ_LF_t(df::DataFrame, l, L)

    op_LF = 0.0

    for y in l+0.1:l:L-0.1

        N_L = nrow(df[(y-l .≤ df.y .≤ y) .& (df.LF_ID .== 1), :])
        N_R = nrow(df[(y-l .≤ df.y .≤ y) .& (df.LF_ID .== -1), :])

        op_LF+=ϕ_LF(N_L, N_R)

    end

    op_LF/length(l+0.1:l:L-0.1)
end


#Iterate over Rows & Overlap
function ϕ_LF(df::DataFrame, l, l_a, Δt, L)

    frames = sort!(unique(df.Frame))
    op_LF = fill(0.0, length(frames))


    for (i, f) in enumerate(frames)
        op_LF[i] = ϕ_LF_t(df[(df.Frame .== f), :], l, l_a, L)
    end

    op_LF
end


function ϕ_LF_t(df::DataFrame, l, l_a, L)

    op_LF = 0.0

    for y in l/2+0.1:l:L-0.1

        N_L = nrow(df[(InLane.(df.y, y, l_a, l) .== true) .& (df.LF_ID .== 1), :])
        N_R = nrow(df[(InLane.(df.y, y, l_a, l) .== true) .& (df.LF_ID .== -1), :])

        op_LF+=ϕ_LF(N_L, N_R)

    end

    op_LF/length(l/2+0.1:l:L-0.1)
end



function InLane(y, y_l, l, l_a)

    if abs(y + l_a/2 - y_l) ≤ l/2 || abs(y - l_a/2 - y_l) ≤ l/2
        true
    else
        false
    end
end
