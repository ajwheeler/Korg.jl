"""
    rayleigh(λs, nH_I, nHe_II, nH2)

Absorption coefficient from Rayleigh scattering by neutral H, He, and H2.  Formulations for H and He
are via [Colgan+ 2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...817..116C/abstract).  Formulation
for H2 from
[Dalgarno and Williams 1962](https://ui.adsabs.harvard.edu/abs/1962ApJ...136..690D/abstract).

The Dalgarno and Williams H2 is applicable redward of 1300 Å.  Since Rayleigh scattering breaks down
when the particle size to wavelength ratio gets large, we that all frequencies passed to this
function be equivalent to 1300 Å or greater.

The formulations for H is adapted from
[Lee 2005](https://ui.adsabs.harvard.edu/abs/1962ApJ...136..690D/abstract), which states that it is
applicable redward of Lyman alpha. See Colgan 2016 for details on He.
"""
function rayleigh(νs::AbstractVector{<:Real}, nH_I, nHe_I, nH2)
    @assert c_cgs / maximum(νs) > 1.3e-5
    σth = 6.65246e-25 #Thompson scattering cross section [cm^2]

    #(ħω/ 2E_H)^2 in Colgan+ 2016.  The photon energy over 2Ryd
    E_2Ryd_2 = @. (hplanck_eV * νs / (2 * Rydberg_eV))^2
    E_2Ryd_4 = E_2Ryd_2 .^ 2
    E_2Ryd_6 = E_2Ryd_2 .* E_2Ryd_4
    E_2Ryd_8 = E_2Ryd_4 .^ 2

    #Colgan+ 2016 equation 6
    σH_σth = @. 20.24E_2Ryd_4 + 239.2E_2Ryd_6 + 2256E_2Ryd_8
    #Colgan+ 2016 equation 7
    σHe_σth = @. 1.913E_2Ryd_4 + 4.52E_2Ryd_6 + 7.90E_2Ryd_8

    αH_HE = @. (nH_I * σH_σth + nHe_I * σHe_σth) * σth

    #Dalgarno & Williams 1962 equation 3 (which assumes λ in Å)
    invλ2 = (νs ./ (1e8 * c_cgs)) .^ 2
    invλ4 = invλ2 .^ 2
    invλ6 = invλ2 .* invλ4
    invλ8 = invλ4 .^ 2
    αH2 = @. (8.14e-13 * invλ4 + 1.28e-6 * invλ6 + 1.61 * invλ8) * nH2

    αH_HE .+ αH2
end

"""
    electron_scattering(nₑ)

Compute the linear absorption coefficient, α, from scattering off of free electrons. This has no
wavelength dependence. It assumes isotropic scattering.  (See, e.g. Gray p 160.)

# Arguments

  - `nₑ::F`: number density of free electrons (in cgs)
"""
function electron_scattering(nₑ::F) where {F<:Real}
    8π / 3 * (electron_charge_cgs^2 / (electron_mass_cgs * c_cgs^2))^2 * nₑ
end
