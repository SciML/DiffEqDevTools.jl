module DiffEqDevTools

using DiffEqBase, RecipesBase, OrdinaryDiffEq

import Base: length

const TIMESERIES_ERRORS = Set([:l2,:l∞,:L2,:L∞])
const DENSE_ERRORS = Set([:L2,:L∞])

include("benchmark.jl")
include("convergence.jl")
include("plotrecipes.jl")
include("test_solution.jl")

export ConvergenceSimulation, Shootout, ShootoutSet, TestSolution

#Benchmark Functions
export ode_shootout, ode_shootoutset, ode_workprecision, ode_workprecision_set

export test_convergence, appxtrue!

#Plot Functions
export stability_region

#Tableus
export constructEuler, constructKutta3, constructRK4, constructRK438Rule,
       constructImplicitEuler, constructMidpointRule, constructTrapezoidalRule,
       constructLobattoIIIA4, constructLobattoIIIB2, constructLobattoIIIB4,
       constructLobattoIIIC2, constructLobattoIIIC4, constructLobattoIIICStar2,
       constructLobattoIIICStar4, constructLobattoIIID2, constructLobattoIIID4,
       constructRadauIA3, constructRadauIA5, constructRadauIIA3, constructRadauIIA5,
       constructRalston, constructHeun, constructRKF5, constructBogakiShampine3,
       constructCashKarp, constructRKF8, constructDormandPrince8,
       constructMSRI1,constructFeagin10, constructFeagin12, constructFeagin14,
       constructDormandPrince8_64bit, constructRKF5, constructRungeFirst5,
       constructCassity5, constructLawson5, constructLutherKonen5, constructLutherKonen52,
       constructLutherKonen53, constructPapakostasPapaGeorgiou5, constructPapakostasPapaGeorgiou52,
       constructTsitouras5, constructBogakiShampine5, constructSharpSmart5,
       constructButcher6, constructButcher7, constructDverk, constructClassicVerner6,
       constructClassicVerner7, constructClassicVerner8, constructClassicVerner92,
       constructVernerRobust7, constructEnrightVerner7, constructTanakaYamashitaStable7,
       constructTanakaYamashitaEfficient7, constructSharpSmart7, constructSharpVerner7,
       constructVernerEfficient7, constructCooperVerner8, constructCooperVerner82,
       constructTsitourasPapakostas8, constructdverk78, constructEnrightVerner8,
       constructCurtis8, constructVernerRobust9, constructVernerEfficient9,
       constructSharp9, constructTsitouras9, constructTsitouras92,constructFeagin14Tableau,
       constructFeagin12Tableau, constructOno12, constructCurtis10, constructOno10, constructFeagin10Tableau,
       constructCurtis10, constructBaker10, constructHairer10, constructButcher63,
       constructButcher6, constructButcher62, constructVerner6, constructDormandPrince6,
       constructSharpVerner6, constructVerner9162, constructVerner916, constructVernerRobust6,
       constructVernerEfficient6, constructPapakostas6, constructLawson6,
       constructTsitourasPapakostas6, constructDormandLockyerMcCorriganPrince6,
       constructTanakaKasugaYamashitaYazaki6D, constructTanakaKasugaYamashitaYazaki6C,
       constructTanakaKasugaYamashitaYazaki6B, constructTanakaKasugaYamashitaYazaki6A,
       constructMikkawyEisa, constructChummund6, constructChummund62,
       constructHuta62, constructHuta6, constructRKF4

end # module
