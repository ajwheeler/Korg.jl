using DSP: conv
using ForwardDiff

"""
    autodiffable_conv(f, g)

Compute the convolution of two vectors, `f` and `g`. This is a wrapped version of `DSP.conv` that
handles `ForwardDiff.Dual` types via the chain rule.
"""
function autodiffable_conv(f, g)
    conv(f, g)
end

function autodiffable_conv(f::Vector{ForwardDiff.Dual{T, V, P}}, 
                           g::Vector{ForwardDiff.Dual{T, V, P}}) where {T, V, P}
    vf = ForwardDiff.value.(f)
    vg = ForwardDiff.value.(g)
    vconv = conv(vf, vg) # do convolution on the values
                        
    pf = ForwardDiff.partials.(f)
    pg = ForwardDiff.partials.(g)

    # apply the chain rule to get the partials of the convolution wrt the parameters
    # P is the number of partial dims in the dual type (see signature)
    pconv = Matrix{V}(undef, length(vconv), P)
    for i in 1:P 
        pconv[:, i] = conv(vf, [p[i] for p in pg]) + conv([p[i] for p in pf], vg)
    end

    map(vconv, eachrow(pconv)) do v, p
        ForwardDiff.Dual{T}(v, p...)
    end
end

# handle cases where only one array contains duals
autodiffable_conv(f::Vector{ForwardDiff.Dual{T, V, P}}, g) where {T, V, P} = autodiffable_conv(promote(f, g)...)
autodiffable_conv(f, g::Vector{ForwardDiff.Dual{T, V, P}}) where {T, V, P} = autodiffable_conv(promote(f, g)...)