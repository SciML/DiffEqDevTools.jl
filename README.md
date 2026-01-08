# DiffEqDevTools.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/DiffEqDevTools.jl/workflows/CI/badge.svg)](https://github.com/SciML/DiffEqDevTools.jl/actions?query=workflow%3ACI)
[![Coverage Status](https://coveralls.io/repos/github/SciML/DiffEqDevTools.jl/badge.svg)](https://coveralls.io/github/SciML/DiffEqDevTools.jl)
[![codecov.io](http://codecov.io/github/SciML/DiffEqDevTools.jl/coverage.svg?branch=master)](http://codecov.io/github/SciML/DiffEqDevTools.jl?branch=master)

DiffEqDevTools.jl is a component package in the DifferentialEquations ecosystem. It provides tools for:

- **Convergence testing**: Verify that numerical methods achieve their theoretical convergence orders
- **Benchmarking**: Compare performance and accuracy of different algorithms
- **Error approximation**: Calculate errors when analytical solutions are unavailable

This package is primarily used for testing and development of differential equation solvers, but can also be useful for end users who want to rigorously validate solver performance for their problems.

## Installation

```julia
using Pkg
Pkg.add("DiffEqDevTools")
```

Or for the development version:

```julia
Pkg.add(url="https://github.com/SciML/DiffEqDevTools.jl")
```

## Features

### Convergence Testing

Test that a numerical method achieves its expected order of convergence:

```julia
using OrdinaryDiffEq, DiffEqDevTools

# Define a problem with known analytical solution
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)

# Test convergence across multiple timesteps
dts = 1 ./ 2 .^ (6:10)
sim = test_convergence(dts, prob, Tsit5())

# Check the estimated convergence order
println("Estimated order: ", sim.𝒪est[:final])  # Should be close to 5
```

For problems without analytical solutions, use `analyticless_test_convergence`:

```julia
# Use a high-accuracy reference solution
dts = 1 ./ 2 .^ (6:9)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)
sim = analyticless_test_convergence(dts, prob, Tsit5(), test_setup)
```

### Work-Precision Diagrams

Compare the efficiency (work vs. accuracy) of different algorithms:

```julia
using Plots

abstols = 1 ./ 10 .^ (3:10)
reltols = 1 ./ 10 .^ (3:10)

setups = [
    Dict(:alg => Tsit5()),
    Dict(:alg => Vern7()),
    Dict(:alg => Rodas5())
]

wp = WorkPrecisionSet(prob, abstols, reltols, setups)
plot(wp)  # Generates work-precision diagram
```

### Algorithm Shootouts

Quickly compare algorithms at a fixed timestep:

```julia
setups = [
    Dict(:alg => Euler()),
    Dict(:alg => Midpoint()),
    Dict(:alg => RK4())
]

shoot = Shootout(prob, setups, dt = 1 / 2^4)
println("Winner: ", shoot.winner)
```

### TestSolution for Error Calculation

Create reference solutions for computing errors:

```julia
# Generate high-accuracy reference
ref_sol = solve(prob, Vern9(), abstol = 1e-14, reltol = 1e-14)
test_sol = TestSolution(ref_sol)

# Use it to compute errors for other solutions
sol = solve(prob, Tsit5(), adaptive = false, dt = 0.1)
err_sol = appxtrue(sol, test_sol)
println("Final error: ", err_sol.errors[:final])
```

## Stochastic Differential Equations

DiffEqDevTools also supports convergence testing for SDEs and other stochastic problems:

```julia
using StochasticDiffEq

f(du, u, p, t) = (du .= 1.01u)
g(du, u, p, t) = (du .= 0.87u)
prob = SDEProblem(f, g, [1.0], (0.0, 1.0))

dts = (1 / 2) .^ (7:-1:4)
test_dt = 1 / 2^8
sim = analyticless_test_convergence(dts, prob, SRIW1(), test_dt, trajectories = 100)
```

## Exported Functions

**Main Functions:**
- `test_convergence`: Test convergence with analytical solution
- `analyticless_test_convergence`: Test convergence without analytical solution
- `WorkPrecision`, `WorkPrecisionSet`: Create work-precision diagrams
- `Shootout`, `ShootoutSet`: Algorithm comparison at fixed timesteps
- `TestSolution`: Wrapper for reference solutions
- `appxtrue`, `appxtrue!`: Compute errors against reference solutions

**Tableau Functions:**
- `deduce_Butcher_tableau`: Extract Butcher tableau from an algorithm
- `stability_region`: Plot stability region of a method
- `check_tableau`: Verify order conditions of a Runge-Kutta tableau
- Constructors for various RK tableaus (e.g., `constructEuler`, `constructRK4`, etc.)

## Documentation

For more detailed documentation, please see the [DifferentialEquations.jl documentation](https://docs.sciml.ai/DiffEqDocs/stable/). For questions or issues, please visit the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby) or open an issue on [GitHub](https://github.com/SciML/DiffEqDevTools.jl/issues).

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests. When contributing, please ensure that:
- Code follows the existing style
- Tests are added for new functionality
- Documentation is updated as needed

## Citation

If you use DiffEqDevTools.jl in your research, please cite the DifferentialEquations.jl project. See [CITATION.bib](CITATION.bib) for the BibTeX entry.
