using CSV, DataFrames, Statistics

include("./Init_DataFrame.jl")
include("./Intruder.jl")
include("./Vectors.jl")
include("./TTC.jl")
include("./TimeGap.jl")
include("./Pair_Distribution.jl")
include("./Lane_Formation.jl")
include("./Plot_Video.jl")
include("./Dimless_Numbers.jl")
include("./Pairs.jl")


function Keep_only(Files, suff)

    Files_ = Vector{String}()

    for (i, f) in enumerate(Files)

        if endswith(f, suff) == true
            push!(Files_, f)
        end
    end

    Files_
end
;
