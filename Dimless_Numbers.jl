function Calc_IN_AV(path, Files, r_min, t_min, t_max, r_soc, dia, T)

    IN_ = fill(0.0, length(Files))
    σ_IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    σ_AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)
        #INs = IN_Add(df, r_soc, dia)
        INs = IN_Max(df, r_soc, dia)

        #INs = IN(df[df.r .> r_min, :].r, r_soc, dia)
        #INs = 1 ./df[t_max .> df.TG .> t_min, :].TG
        #r_Front = df[t_max .> df.TG .> t_min, :].TG .+ dia
        #INs = IN(r_Front, r_soc, dia)
        #AVs = AV_Add(df, T, dia)

        AVs = AV(df[t_max .> df.TTC .> t_min, :].TTC, T)

        IN_[i], σ_IN_[i] = mean(INs), std(INs)
        AV_[i], σ_AV_[i] = mean(AVs), std(AVs)

        rho_[i] = ρ_Global(df)

    end

    IN_ ,σ_IN_, AV_, σ_AV_, rho_

end


function Calc_IN_AV_NEGLECT_PAIRS(path, Files, r_min, t_min, t_max, r_soc, dia, T, f_min, d_mean, d_max)

    IN_ = fill(0.0, length(Files))
    σ_IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    σ_AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)

        neglect_IDs = Calc_Neglect_Pair_IDs(df, f_min, d_mean, d_max)
        df = filter(row -> row.ID ∉ neglect_IDs, df);

        #INs = IN_Add(df, r_soc, dia)
        INs = IN_Max(df, r_soc, dia)

        #INs = IN(df[df.r .> r_min, :].r, r_soc, dia)
        #INs = 1 ./df[t_max .> df.TG .> t_min, :].TG
        #r_Front = df[t_max .> df.TG .> t_min, :].TG .+ dia
        #INs = IN(r_Front, r_soc, dia)

        AVs = AV(df[t_max .> df.TTC .> t_min, :].TTC, T)
        #AVs = AV_Add(df, T, dia)


        IN_[i], σ_IN_[i] = mean(INs), std(INs)
        AV_[i], σ_AV_[i] = mean(AVs), std(AVs)

        rho_[i] = ρ_Global(df)

    end

    IN_ ,σ_IN_, AV_, σ_AV_, rho_

end

function Calc_IN_AV_1D(path, Files, r_min, t_min, t_max, r_soc, dia, T)

    IN_ = fill(0.0, length(Files))
    σ_IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    σ_AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)

        #rename!(df,4 => "x", 3 => "y")

        #INs = IN(df[df.r .> r_min, :].r, r_soc, dia)
        #INs = IN_Add_1d(df, r_soc, dia)
        INs = IN_Max_1d(df, r_soc, dia)
        #INs = 1 ./df[t_max .> df.TG .> t_min, :].TG
        #r_Front = df[t_max .> df.TG .> t_min, :].TG .+ dia
        #INs = IN(r_Front, r_soc, dia)

        #AVs = AV_Add_1d(df, T, dia)

        AVs = AV(df[t_max .> df.TTC .> t_min, :].TTC, T)
        #AVs = AV(df[t_max .> df.TTC .> t_min, :].TTC, T)

        IN_[i], σ_IN_[i] = mean(INs), std(INs)
        AV_[i], σ_AV_[i] = mean(AVs), std(AVs)

        rho_[i] = ρ_Global(df)

    end

    IN_ ,σ_IN_, AV_, σ_AV_, rho_

end

IN(r::Vector, r_soc, dia) = [IN(r_i, r_soc, dia) for r_i in r]

IN(r::Float64, r_soc, dia) = ((r_soc - dia)/(r-dia))^2

AV(ttc::Vector, T) = (T./ttc)#.^2# * exp.(-1 .* ttc./T)

AV(ttc::Float64, T) = (T/ttc)#^2# * exp(-1 * ttc/T)

function AV_Add(x, v, df, T, dia)

    av_ = 0.0

    for i in 1:length(df.x)

        ttc_ = TTC(x, (df.x[i], df.y[i]), v, (df.v_x[i], df.v_y[i]), dia, dia)

        if ttc_ > 0.01 && ttc_ < 999

            av_ += AV(ttc_, T)

        end

    end

    av_

end

function AV_Add(df, T, dia)

    avs = fill(0.0, nrow(df))
    index = 1

    gdf = groupby(df, :Frame)

    for df_ in gdf

        gdf_ = groupby(df_, :ID)

        for df_2 in gdf_

            x = (df_2.x[1], df_2.y[1])
            v = (df.v_x[1], df.v_y[1])

            avs[index] = AV_Add(x, v, df_, T, dia)
            index += 1

        end

    end

    avs

end

function AV_Add_1d(x, v, df, T, dia)

    av_ = 0.0

    for i in 1:length(df.x)

        ttc_ = TTC(x, df.x[i], v, df.v_x[i], dia, dia)

        if ttc_ > 0.01 && ttc_ < 999

            av_ += AV(ttc_, T)

        end

    end

    av_

end


function AV_Add_1d(df, T, dia)
    avs = fill(0.0, nrow(df))
    index = 1

    gdf = groupby(df, :Frame)

    for df_ in gdf

        gdf_ = groupby(df_, :ID)

        for df_2 in gdf_

            x = df_2.x[1]
            v = df_2.v_x[1]

            avs[index] = AV_Add_1d(x, v, df_, T, dia)
            index += 1


        end

    end

    avs

end

function IN_Add(x, df, r_soc, dia)

    in_ = 0.0

    for i in 1:length(df.x)

        d_ = d(x, (df.x[i], df.y[i]))

        if d_ > dia + 0.01

            in_ += IN(d_, r_soc, dia)

        end

    end

    in_

end


function IN_Add(df, r_soc, dia)
    ins = fill(0.0, nrow(df))
    index = 1

    gdf = groupby(df, :Frame)

    for df_ in gdf

        gdf_ = groupby(df_, :ID)

        for df_2 in gdf_

            x = (df_2.x[1], df_2.y[1])

            ins[index] = IN_Add(x, df_, r_soc, dia)
            index += 1


        end

    end

    ins

end

function IN_Max(x, df, r_soc, dia)

    in_ = 0.0

    for i in 1:length(df.x)

        d_ = d(x, (df.x[i], df.y[i]))

        if d_ > dia + 0.01

            in_ = max(in_, IN(d_, r_soc, dia))

        end

    end

    in_

end

function IN_Max(df, r_soc, dia)
    ins = fill(0.0, nrow(df))
    index = 1

    gdf = groupby(df, :Frame)

    for df_ in gdf

        gdf_ = groupby(df_, :ID)

        for df_2 in gdf_

            x = (df_2.x[1], df_2.y[1])

            ins[index] = IN_Max(x, df_, r_soc, dia)
            index += 1


        end

    end

    ins

end

function IN_Add_1d(x, df, r_soc, dia)

    in_ = 0.0

    for i in 1:length(df.x)

        d_ = d(x, df.x[i])

        if d_ > dia + 0.01

            in_ += IN(d_, r_soc, dia)

        end

    end

    in_

end


function IN_Add_1d(df, r_soc, dia)
    ins = fill(0.0, nrow(df))
    index = 1

    gdf = groupby(df, :Frame)

    for df_ in gdf

        gdf_ = groupby(df_, :ID)

        for df_2 in gdf_

            x = df_2.x[1]

            ins[index] = IN_Add_1d(x, df_, r_soc, dia)
            index += 1


        end

    end

    ins

end

function IN_Max_1d(x, df, r_soc, dia)

    in_ = 0.0

    for i in 1:length(df.x)

        d_ = d(x, df.x[i])

        if d_ > dia + 0.01

            in_ = max(IN(d_, r_soc, dia), in_)

        end

    end

    in_

end


function IN_Max_1d(df, r_soc, dia)
    ins = fill(0.0, nrow(df))
    index = 1

    gdf = groupby(df, :Frame)

    for df_ in gdf

        gdf_ = groupby(df_, :ID)

        for df_2 in gdf_

            x = df_2.x[1]

            ins[index] = IN_Max_1d(x, df_, r_soc, dia)
            index += 1


        end

    end

    ins

end

function Mean_AV_DataSet(df, T, t_min, t_max)

    AV_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        av_help = AV(df_f, T, t_min, t_max)

        if isnan(av_help) == false

            AV_mean += av_help
            counter += 1

        end

    end

    AV_mean/counter
end

AV(df::SubDataFrame, T, t_min, t_max) = mean(T./filter(row -> t_max > row.TTC > t_min, df).TTC)


function Mean_IN_DataSet(df, r_min, r_soc, l_min)

    IN_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        in_help = IN(df_f, r_min, r_soc, l_min)

        if isnan(in_help) == false

            IN_mean += in_help
            counter += 1

        end

    end

    IN_mean/counter
end

IN(df::SubDataFrame, r_min, r_soc, l_min) = mean(((r_soc - l_min)./(filter(row -> row.r > r_min, df).r .- l_min)).^2)

function Calc_IN_AV(path, Files, r_min, t_min, t_max, r_soc, l_min, T)

    IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)

        IN_[i] = Mean_IN_DataSet(df, r_min, r_soc, l_min)
        AV_[i] = Mean_AV_DataSet(df, T, t_min, t_max)
        rho_[i] = ρ_Global(df)

    end

    IN_ , AV_, rho_

end

function Calc_IN_AV_NEGLECT_PAIRS(path, Files, r_min, t_min, t_max, r_soc, l_min, T, f_min, d_mean, d_max)

    IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)

        #exclude all pairs
        neglect_IDs = Calc_Neglect_Pair_IDs(df, f_min, d_mean, d_max)
        df = filter(row -> row.ID ∉ neglect_IDs, df);

        #calculate the average values averaged over the frames
        IN_[i] = Mean_IN_DataSet(df, r_min, r_soc, l_min)
        AV_[i] = Mean_AV_DataSet(df, T, t_min, t_max)
        rho_[i] = ρ_Global(df)

    end

    IN_, AV_, rho_

end

function Calc_IN_AV_1D(path, Files, r_min, t_min, t_max, r_soc, l_min, T)

    IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)

        #calculate the average values averaged over the frames
        IN_[i] = Mean_IN_DataSet(df, r_min, r_soc, l_min)
        AV_[i] = Mean_AV_DataSet(df, T, t_min, t_max)
        rho_[i] = ρ_Global(df)

    end

    IN_ , AV_, rho_

end
