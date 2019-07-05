using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
using DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear,
       prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear

probArr = Vector{ODEProblem}(undef, 2)
bigprobArr = Vector{ODEProblem}(undef, 2)

probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear
bigprobArr[1] = prob_ode_bigfloatlinear
bigprobArr[2] = prob_ode_bigfloat2Dlinear
setprecision(400)
Random.seed!(100)
## Convergence Testing
println("Convergence Test on Linear")
dts = 1 .//2 .^(8:-1:4)
testTol = 0.3
superduperbool = Vector{Bool}(undef, 2)

for i = 1:2 # 1 = num, 2 = ExplicitRK
  global dts
  if i>1
    prob = probArr[2]
    bigprob = bigprobArr[2]
  else
    prob = probArr[1]
    bigprob = bigprobArr[1]
  end


  # Order 2

  tabalg = ExplicitRK(tableau=constructHeun())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-2) < testTol

  tabalg = ExplicitRK(tableau=constructRalston())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-2) < testTol

  tabalg = ExplicitRK(tableau=constructSSPRK22())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-2) < testTol

  # Order 3

  tabalg = ExplicitRK(tableau=constructBogakiShampine3())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-3) < testTol

  tabalg = ExplicitRK(tableau=constructSSPRK33())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-3) < testTol

  tabalg = ExplicitRK(tableau=constructSSPRK43())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-3) < testTol

  # Order 4

  tabalg = ExplicitRK(tableau=constructRKF4())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-4) < testTol

  tabalg = ExplicitRK(tableau=constructSSPRK104())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-4) < testTol

  # Order 5

  dts = 1 .//2 .^(7:-1:4)
  tabalg = ExplicitRK(tableau=constructRKF5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  dts = 1 .//2 .^(7:-1:4)
  tabalg = ExplicitRK(tableau=constructDormandPrince())
  sim = test_convergence(dts,prob,tabalg)
  sim2 = test_convergence(dts,prob,DP5())
  @test (abs(sim.𝒪est[:l∞]-5) < testTol && (maximum(sim[end][end]-sim2[end][end]) < 1e-10))

  tabalg = ExplicitRK(tableau=constructCashKarp())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  dts = 1 .//2 .^(7:-1:4)
  tabalg = ExplicitRK(tableau=constructRungeFirst5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructCassity5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructLawson5())
  sim = test_convergence(dts,prob,tabalg) #10
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructLutherKonen5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructLutherKonen52())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructLutherKonen53())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructPapakostasPapaGeorgiou5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructPapakostasPapaGeorgiou52())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  dts = 1 .//2 .^(6:-1:4)
  tabalg = ExplicitRK(tableau=constructTsitouras5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  dts = 1 .//2 .^(6:-1:4)
  tabalg = ExplicitRK(tableau=constructBogakiShampine5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol

  tabalg = ExplicitRK(tableau=constructSharpSmart5())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5) < testTol


  # Order 6

  dts = 1 .//2 .^(6:-1:4)
  tabalg = ExplicitRK(tableau=constructButcher6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructButcher62())
  sim = test_convergence(dts,prob,tabalg) #20
  @test abs(sim.𝒪est[:l∞]-6) < testTol # Less stringent

  dts = 1 .//2 .^(6:-1:4)
  tabalg = ExplicitRK(tableau=constructButcher63())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructDormandPrince6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol+0.1 # Better on linear

  tabalg = ExplicitRK(tableau=constructSharpVerner6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructVerner916())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructVerner9162())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructVernerRobust6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructVernerEfficient6(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6.6) < testTol

  tabalg = ExplicitRK(tableau=constructPapakostas6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructLawson6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  dts = 1 .//2 .^(3:-1:1)
  tabalg = ExplicitRK(tableau=constructTsitourasPapakostas6())
  sim = test_convergence(dts,prob,tabalg) #30
  @test abs(sim.𝒪est[:l∞]-6.7) < testTol # Better on linear

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructDormandLockyerMcCorriganPrince6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructTanakaKasugaYamashitaYazaki6D())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol


  tabalg = ExplicitRK(tableau=constructTanakaKasugaYamashitaYazaki6C())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructTanakaKasugaYamashitaYazaki6B())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructTanakaKasugaYamashitaYazaki6A())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructMikkawyEisa())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6.53) < testTol # Odd behavior

  tabalg = ExplicitRK(tableau=constructChummund6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructChummund62())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructHuta6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-5.5) < testTol # Low convergence, error noted in Stone notes

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructHuta62())
  sim = test_convergence(dts,prob,tabalg)#40
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructVerner6())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6.7) < testTol # Better on linear

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructDverk())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  tabalg = ExplicitRK(tableau=constructClassicVerner6())
  sim = test_convergence(dts,prob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6) < testTol

  # Order 7

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructButcher7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol

  dts = 1 .//2 .^(5:-1:2)
  tabalg = ExplicitRK(tableau=constructClassicVerner7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol

  tabalg = ExplicitRK(tableau=constructVernerRobust7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol


  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructEnrightVerner7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7.15) < testTol # Better on linear

  tabalg = ExplicitRK(tableau=constructTanakaYamashitaStable7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7.3) < testTol

  tabalg = ExplicitRK(tableau=constructTanakaYamashitaEfficient7(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol

  dts = 1 .//2 .^(8:-1:3)
  tabalg = ExplicitRK(tableau=constructSharpSmart7(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg) #50
  @test abs(sim.𝒪est[:l∞]-7) < testTol

  dts = 1 .//2 .^(3:-1:1)
  tabalg = ExplicitRK(tableau=constructSharpVerner7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-6.5) < testTol # Coefficients aren't accurate enough, drop off error

  tabalg = ExplicitRK(tableau=constructVernerEfficient7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol

  tabalg = ExplicitRK(tableau=constructVerner7())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-7) < testTol

  # Order 8
  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructClassicVerner8())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructVerner8(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructCooperVerner8())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol #Coefficients not accurate enough

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructCooperVerner82())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol #Coefficients not accurate enough

  tabalg = ExplicitRK(tableau=constructTsitourasPapakostas8(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructdverk78(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructEnrightVerner8(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructCurtis8())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  dts = 1 .//2 .^(4:-1:1)
  tabalg = ExplicitRK(tableau=constructRKF8(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8) < testTol

  tabalg = ExplicitRK(tableau=constructDormandPrince8(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8.4) < testTol

  dts = 1 .//2 .^(3:-1:1)
  tabalg = ExplicitRK(tableau=constructDormandPrince8_64bit(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-8.4) < testTol

  # Order 9

  tabalg = ExplicitRK(tableau=constructVernerRobust9(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-9) < testTol

  tabalg = ExplicitRK(tableau=constructVernerEfficient9(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-9) < testTol

  dts = 1 .//2 .^(3:-1:1)
  tabalg = ExplicitRK(tableau=constructSharp9(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-9) < testTol #Only works to Float64 precision

  dts = 1 .//2 .^(2:-1:1)
  tabalg = ExplicitRK(tableau=constructTsitouras9(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10.5) < testTol #Only works to Float64

  dts = 1 .//2 .^(3:-1:1)
  tabalg = ExplicitRK(tableau=constructTsitouras92(BigFloat))
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-9) < testTol  #Only works to Float64

  ## Order 10

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructCurtis10())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10) < testTol

  tabalg = ExplicitRK(tableau=constructOno10())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10) < testTol

  dts = 1 .//2 .^(5:-1:1)
  tabalg = ExplicitRK(tableau=constructFeagin10Tableau())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10) < testTol

  tabalg = ExplicitRK(tableau=constructCurtis10())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10) < testTol

  dts = 1 .//2 .^(6:-1:1)
  tabalg = ExplicitRK(tableau=constructBaker10())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10.8) < testTol


  tabalg = ExplicitRK(tableau=constructHairer10())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-10.8) < testTol

  ## Order 12

  dts = 1 .//2 .^(6:-1:1)
  tabalg = ExplicitRK(tableau=constructFeagin12Tableau())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-12.6) < testTol

  tabalg = ExplicitRK(tableau=constructOno12())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-11.6) < testTol

  ## Order 14

  dts = 1 .//2 .^(6:-1:1)
  tabalg = ExplicitRK(tableau=constructFeagin14Tableau())
  sim = test_convergence(dts,bigprob,tabalg)
  @test abs(sim.𝒪est[:l∞]-15.5) < testTol
end
