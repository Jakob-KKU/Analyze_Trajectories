using DSP

function Apply_BW_Filter!(df, frq)

    responsetype = Lowpass(frq)
    designmethod = Butterworth(2)

    gdf = groupby(df, :ID)

    for df_ in gdf

        df_.x = filt(digitalfilter(responsetype, designmethod), df_.x)
        df_.y = filt(digitalfilter(responsetype, designmethod), df_.y)

    end

end
