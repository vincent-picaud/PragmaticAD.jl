using PragmaticAD
using Base.Test

@testset "PragmaticAD" begin

    import PragmaticAD: Tape,afloat_next_index,afloat_count
    
    include("tape.jl")

end;
nothing
