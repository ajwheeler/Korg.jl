"""
Functions for fitting to data.

!!! warning

    This submodule is in beta. It's API may change.
"""
module Fit
using Compat: @compat
@compat public fit_spectrum, ews_to_abundances, ews_to_stellar_parameters

using ..Korg, ForwardDiff

include("fit_via_synthesis.jl")
include("fit_via_EWs.jl")

end # module
