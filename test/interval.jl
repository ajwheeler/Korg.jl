@testset "Interval" begin
    function _test_contained_slice(vals::AbstractVector, interval::Korg.Interval)
        idx = Korg.contained_slice(vals, interval)
        first_ind, last_ind = first(idx), last(idx)

        @assert first_ind >= 1 && last_ind <= length(vals)

        result = Korg.contained.(vals, Ref(interval))

        if all(result)
            @test (first_ind == 1) && (last_ind == length(vals))
        elseif any(result)
            @test last_ind >= first_ind
            @test all(.!result[1:first_ind-1])
            @test all(result[first_ind:last_ind])
            @test all(.!result[last_ind+1:length(vals)])
        else
            @test first_ind == last_ind + 1
        end
    end

    # first make sure that the following cases are caught by the constructor:
    @test_throws AssertionError Korg.Interval(5, 5)
    @test_throws AssertionError Korg.Interval(3, 2)
    @test_throws AssertionError Korg.Interval(Inf, Inf)
    @test_throws AssertionError Korg.Interval(-Inf, -Inf)

    # check contained
    sample = Korg.Interval(3, 10)
    @test !Korg.contained(3, sample)
    @test !Korg.contained(10, sample)
    @test Korg.contained(5.0, sample)
    @test Korg.contained(nextfloat(3.0), sample)
    @test Korg.contained(prevfloat(10.0), sample)

    # check contained_slice
    @testset "contained_slice" begin
        # we consider cases where the slice contains just a single element or multiple elements

        # first, try cases where everything is in-bounds
        _test_contained_slice([6.0], sample) # slice of 1 element
        _test_contained_slice([4.0, 5.0, 6.0, 7.0, 8.0, 9.0], sample) # slice of multiple elements

        # next, try cases where some values are out-of bounds
        _test_contained_slice([1.0, 6.0], sample) # slice of 1 element
        _test_contained_slice([1.0, 6.0, 12.0], sample)
        _test_contained_slice([6.0, 12.0], sample)

        _test_contained_slice([1.0, 2.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], sample)
        _test_contained_slice([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.5, 12.0], sample)
        _test_contained_slice([1.0, 2.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.5, 12.0], sample)

        # lastly, consider cases where all values are out of bounds
        _test_contained_slice([1.0], sample)
        _test_contained_slice([100.0], sample)

        _test_contained_slice([1.0, 2.0, 2.5], sample)
        _test_contained_slice([10.5, 12.5, 100.0], sample)
    end
end
