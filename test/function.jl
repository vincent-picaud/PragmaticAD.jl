@testset "Comparison operators" begin

    x=2.0
    ax=AFloat(x+1)

    @test x!=ax
    @test ax==ax
    @test ax<=ax
    @test x<ax
    @test !(ax<x)
    
end 
