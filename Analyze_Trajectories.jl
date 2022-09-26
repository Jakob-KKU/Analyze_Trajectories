using CSV, DataFrames, Statistics

include("./Init_DataFrame.jl")
include("./Intruder.jl")
include("./Vectors.jl")
include("./Lane_Formation.jl")
include("./Plot_Video.jl")


function Keep_only!(Files, suff)

    for (i, f) in enumerate(Files)

        if endswith(f, suff) == false
            deleteat!(Files, i)
        end
    end
end
;
