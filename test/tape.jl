@testset "Tape constructor" begin

    tape=Tape(Float64)

    @test afloat_next_index(tape)==1
    @test afloat_count(tape)==0
end 
