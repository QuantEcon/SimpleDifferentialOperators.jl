
@testset "Operators under mixed barrier conditions" begin
    @testset "1-element xs" begin
        # right at the elements at xs
        @test findnearestindex([1], 1) == 1

        # values out of the range
        @test findnearestindex([1], 0.8) == 1
        @test findnearestindex([1], 9) == 1
    end

    @testset "2-element xs" begin
        # right at the elements at xs
        @test findnearestindex([1;3], 1) == 1
        @test findnearestindex([1;3], 3) == 2

        # values in between
        @test findnearestindex([1;3], 1.2) == 1
        @test findnearestindex([1;3], 1.5) == 1
        @test findnearestindex([1;3], 1.8) == 1
        @test findnearestindex([1;3], 2) == 2
        @test findnearestindex([1;3], 2.0) == 2
        @test findnearestindex([1;3], 2.3) == 2
        @test findnearestindex([1;3], 2.8) == 2

        # values out of the range
        @test findnearestindex([1;3], 0.8) == 1
        @test findnearestindex([1;3], 9) == 2
    end

    @testset "3-element xs" begin
        # right at the elements at xs
        @test findnearestindex([1;2;3], 1) == 1
        @test findnearestindex([1;2;3], 2) == 2
        @test findnearestindex([1;2;3], 3) == 3

        # values in between
        @test findnearestindex([1;2;3], 1.2) == 1
        @test findnearestindex([1;2;3], 1.5) == 2
        @test findnearestindex([1;2;3], 1.8) == 2
        @test findnearestindex([1;2;3], 2.3) == 2
        @test findnearestindex([1;2;3], 2.8) == 3
        @test findnearestindex([1;2;3], 3.3) == 3

        # values out of the range
        @test findnearestindex([1;2;3], 0.8) == 1
        @test findnearestindex([1;2;3], 9) == 3
    end
end