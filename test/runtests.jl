using PragmaticAD
using Base.Test

@testset "PragmaticAD" begin

    import PragmaticAD: Tape,afloat_next_index,afloat_count,get_tape
    
    include("tape.jl")
    include("array.jl")

end;
nothing
