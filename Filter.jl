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
