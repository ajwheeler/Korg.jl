"""
    rayleigh(λs, nH_I, nHe_II, nH2)

Absorption coefficient from Rayleigh scattering by neutral H and He.  Formulation taken from Colgan+
2016.  The formulation for H is adapted from Lee 2005, which states that it is applicable redward of 
Lyman alpha. The formulation for He is adapted from Dalgarno 1962 and Dalgarno & Kingston 1960, and 
the formulation for H2 is from Dalgarno and Williams 1963 and is applicable redward of 1300 Å.
"""
function rayleigh(νs::AbstractVector{<:Real}, nH_I, nHe_I, nH2)
    σth = 6.65246e-25 #Thompson scattering cross section [cm^2]

    #(ħω/ 2E_H)^2 in Colgan+ 2016.  The photon energy over 2Ryd
    E_2Ryd_2 = @. (hplanck_eV * νs / (2 * Rydberg_eV ))^2
    E_2Ryd_4 = E_2Ryd_2.^2
    E_2Ryd_6 = E_2Ryd_2.*E_2Ryd_4
    E_2Ryd_8 = E_2Ryd_4.^2

    #Colgan+ 2016 equation 6
    σH_σth = @. 20.24E_2Ryd_4 + 239.2E_2Ryd_6 + 2256E_2Ryd_8
    #Colgan+ 2016 equation 7
    σHe_σth = @. 1.913E_2Ryd_4 + 4.52E_2Ryd_6 + 7.90E_2Ryd_8

    αH_HE = @. (nH_I*σH_σth + nHe_I*σHe_σth) * σth

    #Dalgarno & Williams 1962 equation 3 (which assumes λ in Å)
    invλ2 = (νs ./ (1e8 * c_cgs)) .^2
    invλ4 = invλ2 .^ 2
    invλ6 = invλ2 .* invλ4
    invλ8 = invλ4 .^ 2
    αH2 = @. (8.14e-13*invλ4 + 1.28e-6*invλ6 + 1.61*invλ8) * nH2

    αH_HE .+ αH2
end

"""
    electron_scattering(nₑ, ρ)

Compute the opacity from scattering off of free electrons. This has no wavelength dependence. It
assumes isotropic scattering

# Arguments
- `nₑ::F`: number density of free electrons (in cgs)
- `ρ::F`: the mass density

# Notes
I adopted the formula described in section 5.12 of Kurucz (1970) and the equation in the electron
scattering subsection of Gray (2005); the actual coefficient value comes from the latter. It turns
out that the coefficient in Kurucz (1970) has a typo (it's a factor of 10 too large).
"""
electron_scattering(nₑ::F, ρ::F) where {F<:Real} = 0.6648e-24*nₑ/ρ

