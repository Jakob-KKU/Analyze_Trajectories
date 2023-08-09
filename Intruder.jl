ϕ_perp(x, y) = acos(abs(y)/sqrt(x^2+y^2))

function Find_Intruder(df)

    ind = zeros(0)

    for i in unique(df.ID)
        if Is_Intruder(df[df.ID .== i, :]) == true
            push!(ind, i)
        end
    end

    ind

end

function Is_Intruder(df::DataFrame)

    if maximum(df.x) > 8.0 && minimum(df.x) < 2.0
        true
    elseif maximum(df.y) > 8.0 && minimum(df.y) < 2.0
        true
    else
        false
    end

end

ϕ(r::Float64, ξ::Float64) = r ≤ 2*ξ ? exp(-r^2/(2*ξ^2)) : 0.0

function ρ_Field(df, r, L, ξ)

   A, norm = 0.0, 0.0

    for row in eachrow(df)

        if row.ID == 160

            ϕ_ = ϕ(d((row.x, row.y), r, L), 0.17)


        else

            ϕ_ = ϕ(d((row.x, row.y), r, L), ξ)

        end

        A += ϕ_*row.A
        norm += ϕ_

    end

    if norm == 0.0
        0.0
    else
        norm/A
    end

end


function Calc_ρ_Field(df_G, df, L, ξ, dx, dy)

    ID_Intruder = maximum(df.ID)
    frames = minimum(df.Frame):10:maximum(df.Frame)

    #Split data Frame from Intruder and Agents
    df_Obs = df[(df.ID .== ID_Intruder), :]
    #df = df[(df.ID .< ID_Intruder), :]

    A = fill(0.0, nrow(df_G));

    for fr in frames

        P_Obs = df_Obs[(df_Obs.Frame .== fr), :]

        #Calculate relative Positions in current Frame and Filter to relevant positions
        df_fr = df[
            (df.Frame .== fr) .&
            (-dx-2*ξ .≤ (df.x .- P_Obs.x) .≤ dx+2*ξ) .&
            (-dy-2*ξ .≤ (df.y .- P_Obs.y) .≤ dy+2*ξ), :]

        df_fr[!, :x] = df_fr[!, :x] .- P_Obs.x
        df_fr[!, :y] = df_fr[!, :y] .- P_Obs.y;

        for (i, row) in enumerate(eachrow(df_G))

           A[i] = A[i] + ρ_Field(df_fr, (row.x[1], row.y[1]), L, ξ)

        end

    end

    A./length(frames)

end
;
