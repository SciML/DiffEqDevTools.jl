using DiffEqDevTools, Test

@test stability_region(constructDormandPrince6(), initial_guess=-3.5) ≈ -3.95413 rtol=1e-3
@test stability_region(constructTsitourasPapakostas6(), initial_guess=-3.5) ≈ -3.95413 rtol=1e-3
@test stability_region(constructRadauIIA5(), initial_guess=12.) ≈ 11.84 rtol=1e-2
