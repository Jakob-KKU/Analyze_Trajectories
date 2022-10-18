using CSV, DataFrames, Statistics

include("./Init_DataFrame.jl")
include("./Intruder.jl")
include("./Vectors.jl")
include("./TTC.jl")
include("./TimeGap.jl")

include("./Lane_Formation.jl")
include("./Plot_Video.jl")


function Keep_only(Files, suff)

    Files_ = Vector{String}()

    for (i, f) in enumerate(Files)

        if endswith(f, ".csv") == true
            push!(Files_, f)
        end
    end

    Files_
end
;
