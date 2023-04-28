IN(r::Vector, r_soc, dia) = [IN(r_i, r_soc, dia) for r_i in r]

IN(r::Float64, r_soc, dia) = ((r_soc - dia)/(r-dia))^2

IN_MIN(df::SubDataFrame, r_min, r_soc, l_min) = mean(((r_soc - l_min)./(filter(row -> row.r > r_min, df).r .- l_min)).^2)

AV(ttc::Vector, T) = (T./ttc)#.^2# * exp.(-1 .* ttc./T)

AV(ttc::Float64, T) = (T/ttc)#^2# * exp(-1 * ttc/T)

#AV_MIN(df::SubDataFrame, T, t_min, t_max) = mean(T./filter(row -> t_max > row.TTC > t_min, df).TTC)
AV_MIN(df::SubDataFrame, T, t_min, t_max) = mean((T./filter(row -> row.TTC > t_min, df).TTC).^2)


function Mean_AV_DataSet_MIN(df, T, t_min, t_max)

    AV_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        av_help = AV_MIN(df_f, T, t_min, t_max)

        if isnan(av_help) == false

            AV_mean += av_help
            counter += 1
            #println(counter)

        end

    end

    if counter == 0
        0.0
    else
        AV_mean/counter
    end
end


function Mean_IN_DataSet_SUM(df, r_min, r_soc, l_min)

    IN_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        in_help = IN_SUM(df_f, r_min, r_soc, l_min)

        if isnan(in_help) == false

            IN_mean += in_help
            counter += 1

        end

    end

    IN_mean/counter
end

function IN_SUM(df::SubDataFrame, r_min, r_soc, l_min)

    in_ = 0.0
    gdf = groupby(df, :ID)

    for df_i in gdf

        x = (df_i.x, df_i.y)
        in_+=IN_SUM(x, df::SubDataFrame, r_min, r_soc, l_min)

    end

    in_/length(gdf)

end


function IN_SUM(x, df::SubDataFrame, r_min, r_soc, l_min)

    in_ = 0.0

    dx = sqrt.((df.x .- x[1]).^2+(df.y .-x[2]).^2)

    for dx_i in dx

        if dx_i > r_min && dx_i < 3*r_soc

            in_+= IN(dx_i, r_soc, dia)

        end

    end

    in_

end



function Mean_IN_DataSet_MIN(df, r_min, r_soc, l_min)

    IN_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        in_help = IN_MIN(df_f, r_min, r_soc, l_min)

        if isnan(in_help) == false

            IN_mean += in_help
            counter += 1

        end

    end

    IN_mean/counter
end




function Calc_IN_AV(path, Files, r_min, t_min, t_max, r_soc, l_min, T)

    IN_ = fill(0.0, length(Files))
    AV_ = fill(0.0, length(Files))
    rho_ = fill(0.0, length(Files))

    for (i,file) in enumerate(Files)

        #read trajectory file
        data = CSV.File(string(path, file); comment="#");

        #create Data Frame and rename Columns
        df = DataFrame(data)

        IN_[i] = Mean_IN_DataSet_SUM(df, r_min, r_soc, l_min)
        AV_[i] = Mean_AV_DataSet_MIN(df, T, t_min, t_max)
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
        IN_[i] = Mean_IN_DataSet_SUM(df, r_min, r_soc, l_min)
        AV_[i] = Mean_AV_DataSet_MIN(df, T, t_min, t_max)
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
        IN_[i] = Mean_IN_DataSet_SUM(df, r_min, r_soc, l_min)
        AV_[i] = Mean_AV_DataSet_MIN(df, T, t_min, t_max)
        rho_[i] = ρ_Global(df)

    end

    IN_ , AV_, rho_

end

function Mean_AV_DataSet_SUM(df, T, t_min, t_max)

    AV_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        av_help = AV_SUM(df_f, T, t_min, t_max)

        if isnan(av_help) == false

            AV_mean += av_help
            counter += 1

        end

    end

    if counter == 0
        0.0
    else
        AV_mean/counter
    end
end

#calculate for one frame
function AV_SUM(df::SubDataFrame, T, t_min, t_max)

    av_ = 0.0
    counter = 0
    gdf = groupby(df, :ID)

    for df_i in gdf

        x = (df_i.x[1], df_i.y[1])
        v = (df_i.v_x[1], df_i.v_y[1])

        av_help = AV_SUM(x, v, df, T, t_min, t_max)

        if isnan(av_help) == false

            av_ += av_help
            counter += 1

        end



    end

    av_/counter

end

#calculate for one ID in one Frame
function AV_SUM(x, v, df::SubDataFrame, T, t_min, t_max)

    av_ = 0.0

    for i in 1:length(df.x)

        v_j = (df.v_x[i], df.v_y[i])
        x_j = (df.x[i], df.y[i])

        ttc_ = TTC(x, x_j, v, v_j, 0.2, 0.2)

            if ttc_ > t_min && ttc_ < t_max
                av_+= AV(ttc_, T)
            end

    end

    av_

end

function Mean_AV_DataSet_SUM_1d(df, T, t_min, t_max)

    AV_mean = 0.0
    counter = 0

    gdf = groupby(df, :Frame)

    for df_f in gdf

        av_help = AV_SUM_1d(df_f, T, t_min, t_max)

        if isnan(av_help) == false

            AV_mean += av_help
            counter += 1

        end

    end

    if counter == 0
        0.0
    else
        AV_mean/counter
    end
end

function AV_SUM_1d(df::SubDataFrame, T, t_min, t_max)

    av_ = 0.0
    counter = 0
    gdf = groupby(df, :ID)

    for df_i in gdf

        x = df_i.x[1]
        v = df_i.v_x[1]

        av_help = AV_SUM_1d(x, v, df, T, t_min, t_max)

        if isnan(av_help) == false

            av_ += av_help
            counter += 1

        end



    end

    av_/counter

end


function AV_SUM_1d(x, v, df::SubDataFrame, T, t_min, t_max)

    av_ = 0.0

    for i in 1:length(df.x)

        v_j = df.v_x[i]
        x_j = df.x[i]

        ttc_ = TTC(x, x_j, v, v_j, 0.2, 0.2)

            if ttc_ > t_min && ttc_ < t_max
                av_+= AV(ttc_, T)
            end

    end

    av_

end
