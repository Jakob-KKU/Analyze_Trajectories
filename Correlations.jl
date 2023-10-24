function IN_AV_Crosscorr_allTimes(df::DataFrame, Δt, τ_max)

    τ_s = collect(0:Δt:τ_max)
    corr_ = fill(0.0, length(τ_s))

    for (i, τ) in enumerate(τ_s)

        corr_[i] = IN_AV_Crosscorr(df, τ, Δt)

    end

    τ_s, corr_

end

function IN_AV_Crosscorr(df::DataFrame, τ, Δt)

    corr_ = 0.0

    gdf = groupby(df, :ID)

    for df_i in gdf

        corr_ += IN_AV_Crosscorr(df_i, τ, Δt)

    end

    if length(gdf) == 0
        0.0
    else
        corr_/(length(gdf))
    end

end


function IN_AV_Crosscorr(df::SubDataFrame, τ, Δt)

    corr_ = 0.0

    Δf = Int(round(τ/Δt))

    AV_mean, IN_mean = mean(df.AV), mean(df.IN)

    for i in 1:nrow(df)-Δf

        corr_ += (df.AV[i]-AV_mean)*(df.IN[i+Δf]-IN_mean)

    end

    if std(df.AV) == 0 || std(df.IN) == 0 || (nrow(df)-Δf) <= 0
        0.0
    else
        corr_ / ((nrow(df)-Δf)*(std(df.AV)*std(df.IN)))
    end

end

function AutoCorrelation_allTimes(df::DataFrame, x::Symbol, Δt, τ_max)

    τ_s = collect(0:Δt:τ_max)
    corr_ = fill(0.0, length(τ_s))

    for (i, τ) in enumerate(τ_s)

        corr_[i] = Autocorrelation(df, x, τ, Δt)

    end

    τ_s, corr_

end


function Autocorrelation(df::DataFrame, x::Symbol, τ, Δt)

    corr_ = 0.0

    gdf = groupby(df, :ID)

    for df_i in gdf

        corr_ += Autocorrelation(df_i, x, τ, Δt)

    end

    if length(gdf) == 0
        0.0
    else
        corr_/(length(gdf))
    end

end


function Autocorrelation(df::SubDataFrame, x::Symbol, τ, Δt)

    corr_ = 0.0

    Δf = Int(round(τ/Δt))

    x_s = df[!,x]

    x_mean = mean(x_s)

    for i in 1:nrow(df)-Δf

        corr_ += (x_s[i]-x_mean)*(x_s[i+Δf]-x_mean)

    end

    if std(x_s) == 0 || (nrow(df)-Δf) <= 0
        0.0
    else
        corr_ / ((nrow(df)-Δf)*std(x_s)^2)
    end

end
