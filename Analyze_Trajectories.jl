using CSV, DataFrames, Statistics, SpecialFunctions

include("./Init_DataFrame.jl")
include("./Density_Field.jl")
include("./Vectors.jl")
include("./TTC.jl")
include("./TimeGap.jl")
include("./Pair_Distribution.jl")
include("./Lane_Formation.jl")
include("./Plot_Video.jl")
include("./Dimless_Numbers.jl")
include("./Pairs.jl")
include("./Filter.jl")
include("./Velocity_Field.jl")


function Keep_only(Files, suff)

    Files_ = Vector{String}()

    for (i, f) in enumerate(Files)

        if endswith(f, suff) == true
            push!(Files_, f)
        end
    end

    Files_
end

function Open_Dir(dir, Path)

    cd(Path)

    if dir in readdir()

        cd(dir)

    else

        println("Created Directory /", dir, ".")
        mkdir(dir)
        cd(dir)

    end
end
;
