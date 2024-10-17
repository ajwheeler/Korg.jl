@testset "linelist pruning" begin
    linelist = Korg.read_linelist("data/linelists/gfallvac08oct17.stub.dat"; format="kurucz")
    wls = (linelist[1].wl, linelist[end].wl)
    atm = interpolate_marcs(5000, 4.44)
    A_X = format_A_X()

    strong_lines_sorted = Korg.prune_linelist(atm, linelist, A_X, wls; verbose=false)
    strong_lines_unsorted = Korg.prune_linelist(atm, linelist, A_X, wls
                                                ; sort_by_EW=false, verbose=false)
    more_strong_lines = Korg.prune_linelist(atm, linelist, A_X, wls
                                            ; sort_by_EW=false, threshold=0.01, verbose=false)

    @test 0 < length(strong_lines_sorted) < length(more_strong_lines) < length(linelist)

    @test Set(strong_lines_sorted) == Set(strong_lines_unsorted)
    @test issorted(strong_lines_unsorted; by=l -> l.wl)
    @test Set(strong_lines_sorted) ⊆ Set(more_strong_lines)

    merged_lines = Korg.merge_close_lines(strong_lines_sorted)
    @assert issorted(merged_lines; by=first)
end
