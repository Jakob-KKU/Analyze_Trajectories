using DSP

function Return_BW_Filter(df, freq, fps, order)

    responsetype = Lowpass(freq, fs = fps)
    designmethod = Butterworth(order)

    x = filtfilt(digitalfilter(responsetype, designmethod), df.x)
    y = filtfilt(digitalfilter(responsetype, designmethod), df.y)

    x, y

end

function Apply_BW_Filter!(df, freq, fps, order)

    responsetype = Lowpass(freq, fs = fps)
    designmethod = Butterworth(order)

    gdf = groupby(df, :ID)

    for (i, df_) in enumerate(gdf)

        if length(df_.x) > 12

            try

                df_.x = filtfilt(digitalfilter(responsetype, designmethod), df_.x)
                df_.y = filtfilt(digitalfilter(responsetype, designmethod), df_.y)

            catch err

                if isa(err, BoundsError)

                    df_.y = df_.y
                    df_.x = df_.y
                end

            end

        end



        #    if@warn "Not enought data points to apply filter! Used original data instead."

        #end

        #println(i)

    end

end

function Apply_BW_Filter!(df, norm_freq, order)

    responsetype = Lowpass(norm_freq)
    designmethod = Butterworth(order)

    gdf = groupby(df, :ID)

    for df_ in gdf

        df_.x = filtfilt(digitalfilter(responsetype, designmethod), df_.x)
        df_.y = filtfilt(digitalfilter(responsetype, designmethod), df_.y)

    end

end

α_wheight(t, t_start, t_end, frequency) =  max(exp(- 4.0*frequency*(t - t_start)),exp(- 4.0*frequency*(t_end - t)))


function Apply_BW_Filter_WheightedAverage!(df, frequency, fps, order)

    responsetype = Lowpass(frequency, fs = fps)
    designmethod = Butterworth(order)

    gdf = groupby(df, :ID)

    for (i, df_) in enumerate(gdf)

        if length(df_.x) > 12

            try

                x_filt = filtfilt(digitalfilter(responsetype, designmethod), df_.x)
                y_filt = filtfilt(digitalfilter(responsetype, designmethod), df_.y)

                alphas = α_wheight.(df_.Frame./fps, df_.Frame[1]/fps, df_.Frame[end]/fps, frequency)

                df_.x = (1 .- alphas) .* x_filt .+ alphas .* df_.x
                df_.y = (1 .- alphas) .* y_filt .+ alphas .* df_.y

            catch err

                if isa(err, BoundsError)

                    df_.y = df_.y
                    df_.x = df_.y

                end

            end

        end

    end

end
