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

        INs = IN(df[df.r_F .> r_min, :].r_F, r_soc, dia)
        #INs = 1 ./df[t_max .> df.TG .> t_min, :].TG
        #r_Front = df[t_max .> df.TG .> t_min, :].TG .+ dia
        #INs = IN(r_Front, r_soc, dia)

        AVs = AV(df[t_max .> df.TTC .> t_min, :].TTC, T)

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

        INs = IN(df[df.r .> r_min, :].r, r_soc, dia)
        #INs = 1 ./df[t_max .> df.TG .> t_min, :].TG
        #r_Front = df[t_max .> df.TG .> t_min, :].TG .+ dia
        #INs = IN(r_Front, r_soc, dia)

        AVs = AV(df[t_max .> df.TTC .> t_min, :].TTC, T)

        IN_[i], σ_IN_[i] = mean(INs), std(INs)
        AV_[i], σ_AV_[i] = mean(AVs), std(AVs)

        rho_[i] = ρ_Global(df)

    end

    IN_ ,σ_IN_, AV_, σ_AV_, rho_

end
