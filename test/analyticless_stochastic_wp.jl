using StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary, Base.Test

prob = prob_sde_additivesystem
prob = SDEProblem(prob.f,prob.g,prob.u0,(0.0,1.0))

reltols = 1.0./10.0.^(1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [Dict(:alg=>SRIW1())
          Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) + 1))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) + 1))
          Dict(:alg=>SRIW1(),:dts=>1.0./5.0.^((1:length(reltols)) + 1),:adaptive=>false)
          Dict(:alg=>SRA1(),:dts=>1.0./5.0.^((1:length(reltols)) + 1),:adaptive=>false)
          Dict(:alg=>SRA1())
          ]
names = ["SRIW1","EM","RKMil","SRIW1 Fixed","SRA1 Fixed","SRA1"]
test_dt = 0.1
wp = WorkPrecisionSet(prob,abstols,reltols,setups,test_dt;
                                     numruns=5,names=names,error_estimate=:l2)

se = get_sample_errors(prob,numruns=1000)
se = get_sample_errors(prob,numruns=[5;10;25;50])

prob2 = SDEProblem((t,u,du)->prob.f(t,u,du),prob.g,prob.u0,(0.0,1.0))
test_dt = 2/10^5
appxsol_setup = Dict(:alg=>SRIW1(),:abstol=>1e-7,:reltol=>1e-7)
wp = WorkPrecisionSet(prob2,abstols,reltols,setups,test_dt;
                                     appxsol_setup = appxsol_setup,
                                     numruns=5,names=names,error_estimate=:weak_final)

se2 = get_sample_errors(prob2,test_dt,appxsol_setup = appxsol_setup,
                       numruns=5)
se2 = get_sample_errors(prob2,test_dt,appxsol_setup = appxsol_setup,
                       numruns=[5;10;25;50])

@test all(se-se2 .< 1e-1)
