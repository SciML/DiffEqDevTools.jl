# DiffEqDevTools.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/DiffEqDevTools.jl/workflows/CI/badge.svg)](https://github.com/SciML/DiffEqDevTools.jl/actions?query=workflow%3ACI)
[![Coverage Status](https://coveralls.io/repos/github/SciML/DiffEqDevTools.jl/badge.svg)](https://coveralls.io/github/SciML/DiffEqDevTools.jl)
[![codecov.io](http://codecov.io/github/SciML/DiffEqDevTools.jl/coverage.svg?branch=master)](http://codecov.io/github/SciML/DiffEqDevTools.jl?branch=master)

DiffEqDevTools.jl is a component package in the DifferentialEquations ecosystem. It holds the
convergence testing, benchmarking, and error approximation functionality which is
used by the other component packages in order to test for correctness. Users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
