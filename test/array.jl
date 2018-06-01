@testset "Array constructor" begin

    tape=get_tape(Float64)
    asize = afloat_count(tape)
    a=rand(2,3)
    aa=AArray(a)
    asize_next = afloat_count(tape)
    
    @test asize_next == asize + length(a)
    @test aa[2,3]==a[2,3]
end 
