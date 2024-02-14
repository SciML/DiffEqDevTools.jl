"""
    deduce_Butcher_tableau(erk, T=Float64)

Deduce and return the Butcher coefficients `A, b, c` by solving some
specific ordinary differential equations using the explicit Runge-Kutta
method `erk`. The type `T` will be used for computations and is the
`eltype` of `A`, `b`, and `c`.
"""
function deduce_Butcher_tableau(erk, T = Float64)
    # get the number of stages s
    step_counter_s = 0
    function f_for_s(du, u, p, t)
        step_counter_s += 1
        du .= zero(eltype(du))
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = [zero(T)]
    ode = ODEProblem(f_for_s, u0, tspan)
    step_counter_s = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)
    s = step_counter_s

    # get the nodes c[i]
    step_counter_c = 0
    c = zeros(T, s)

    function f_for_c(du, u, p, t)
        step_counter_c += 1
        if step_counter_c <= s
            c[step_counter_c] = t
        end
        du .= zero(eltype(du))
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = [zero(T)]
    ode = ODEProblem(f_for_c, u0, tspan)
    step_counter_c = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)

    # get the weights b[i]
    step_counter_b = 0
    b = zeros(T, s)
    function f_for_b(du, u, p, t)
        step_counter_b += 1
        for idx in 1:length(du)
            du[idx] = step_counter_b == idx
        end
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = zeros(T, s)
    ode = ODEProblem(f_for_b, u0, tspan)
    step_counter_b = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)
    for i in 1:s
        b[i] = integrator.u[i]
    end

    # get the coefficients A[i,j]
    step_counter_A = 0
    A = zeros(T, s, s)
    function f_for_A(du, u, p, t)
        step_counter_A += 1
        for idx in 1:length(du)
            A[step_counter_A, idx] = u[idx]
            du[idx] = step_counter_A == idx
        end
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = zeros(T, s)
    ode = ODEProblem(f_for_A, u0, tspan)
    step_counter_A = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)

    A, b, c
end

#=

Classic Verner 8(9) excluded since it's just hard to write from the paper.
If anyone wants to add them I'd accept them.

=#

"""
Gauss-Legendre Order 2.
"""
function constructGL2(T::Type = Float64)
    c = [1 // 2]
    A = [1 // 2]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 2))
end

"""
Gauss-Legendre Order 4.
"""
function constructGL4(T::Type = Float64)
    c = [(3 - sqrt(3)) / 6; (3 + sqrt(3)) / 6]
    A = [1/4 (3 - 2 * sqrt(3))/12;
         (3 + 2 * sqrt(3))/12 1/4]
    α = [1 / 2; 1 / 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 4))
end

"""
Gauss-Legendre Order 6.
"""
function constructGL6(T::Type = Float64)
    c = [(5 - sqrt(15)) / 10; 1 / 2; (5 + sqrt(15)) / 10]
    A = [5/36 (10 - 3 * sqrt(15))/45 (25 - 6 * sqrt(15))/180;
         (10 + 3 * sqrt(15))/72 2/9 (10 - 3 * sqrt(15))/72;
         (25 + 6 * sqrt(15))/180 (10 + 3 * sqrt(15))/45 5/36]
    α = [5 / 18; 4 / 9; 5 / 18]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 6))
end

"""
Heun's Order 2 method.
"""
function constructHeun(T::Type = Float64)
    A = [0 0
         1 0]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1,
        stability_size = 2.0))
end
"""
Ralston's Order 2 method.
"""
function constructRalston(T::Type = Float64)
    A = [0 0
         2//3 0]
    c = [0; 2 // 3]
    α = [1 // 4; 3 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 2, stability_size = 2.0))
end

"""
Euler's method.
"""
function constructEuler(T::Type = Float64)
    A = Matrix{T}(undef, 1, 1)
    A[1] = 0
    c = [0]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 1, stability_size = 2.0))
end

"""
Kutta's Order 3 method.
"""
function constructKutta3(T::Type = Float64)
    A = [0 0 0
         1//2 0 0
         -1 2 0]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 3, stability_size = 2.5127453266183286))
end

"""
Classic RK4 method.
"""
function constructRK4(T::Type = Float64)
    A = [0 0 0 0
         1//2 0 0 0
         0 1//2 0 0
         0 0 1 0]
    c = [0; 1 // 2; 1 // 2; 1]
    α = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 4, stability_size = 2.785293563405282))
end

"""
Classic RK4 3/8's rule method.
"""
function constructRK438Rule(T::Type = Float64)
    A = [0 0 0 0
         1//3 0 0 0
         -1//3 1 0 0
         1 -1 1 0]
    c = [0; 1 // 3; 2 // 3; 1]
    α = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 4, stability_size = 2.785293563405282))
end

"""
Ralston's Order 4 method with minimum truncation error.
"""
function constructRalston4(T::Type = Float64)
    sqrt5 = sqrt(convert(T, 5))
    a21 = 4 // 10
    a31 = (-2889 + 1428 * sqrt5) / 1024
    a32 = (3785 - 1620 * sqrt5) / 1024
    a41 = (-3365 + 2094 * sqrt5) / 6040
    a42 = (-975 - 3046 * sqrt5) / 2552
    a43 = (467040 + 203968 * sqrt5) / 240845
    A = [0 0 0 0
         a21 0 0 0
         a31 a32 0 0
         a41 a42 a43 0]
    b2 = 4 // 10
    b3 = (14 - 3 * sqrt5) / 16
    c = [0, b2, b3, 1]
    b1 = (263 + 24 * sqrt5) / 1812
    b2 = (125 - 1000 * sqrt5) / 3828
    b3 = 1024 * (3346 + 1623 * sqrt5) / 5924787
    b4 = (30 - 4 * sqrt5) / 123
    α = [b1, b2, b3, b4]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return ExplicitRKTableau(A, c, α, 4, stability_size = 2.7852935634052822)
end

"""
Explicit SSP method of order 2 using 2 stages.
"""
function constructSSPRK22(T::Type = Float64)
    A = [0 0
         1 0]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 2, stability_size = -2.0))
end

"""
Explicit SSP method of order 3 using 3 stages.
"""
function constructSSPRK33(T::Type = Float64)
    A = [0 0 0
         1 0 0
         1//4 1//4 0]
    c = [0; 1; 1 // 2]
    α = [1 // 6; 1 // 6; 2 // 3]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 3, stability_size = 2.5127453266183286))
end

"""
Explicit SSP method of order 3 using 4 stages.
"""
function constructSSPRK43(T::Type = Float64)
    A = [0 0 0 0
         1//2 0 0 0
         1//2 1//2 0 0
         1//6 1//6 1//6 0]
    c = [0; 1 // 2; 1; 1 // 2]
    α = [1 // 6; 1 // 6; 1 // 6; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 3, stability_size = 5.149486147774043))
end

"""
Explicit SSP method of order 4 using 10 stages.
"""
function constructSSPRK104(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0 0 0
         1//6 0 0 0 0 0 0 0 0 0
         1//6 1//6 0 0 0 0 0 0 0 0
         1//6 1//6 1//6 0 0 0 0 0 0 0
         1//6 1//6 1//6 1//6 0 0 0 0 0 0
         1//15 1//15 1//15 1//15 1//15 0 0 0 0 0
         1//15 1//15 1//15 1//15 1//15 1//6 0 0 0 0
         1//15 1//15 1//15 1//15 1//15 1//6 1//6 0 0 0
         1//15 1//15 1//15 1//15 1//15 1//6 1//6 1//6 0 0
         1//15 1//15 1//15 1//15 1//15 1//6 1//6 1//6 1//6 0]
    c = [0; 1 // 6; 1 // 3; 1 // 2; 2 // 3; 1 // 3; 1 // 2; 2 // 3; 5 // 6; 1]
    α = [1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10; 1 // 10;
         1 // 10]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 4, stability_size = 13.917047464637367))
end

"""
Implicit Euler Method
"""
function constructImplicitEuler(T::Type = Float64)
    A = Matrix{T}(undef, 1, 1)
    A[1] = 1
    c = [1]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 1))
end

"""
Order 2 Midpoint Method
"""
function constructMidpointRule(T::Type = Float64)
    A = Matrix{T}(undef, 1, 1)
    A[1] = 1 // 2
    c = [1 // 2]
    α = [1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 2))
end

"""
Order 2 Trapezoidal Rule (LobattoIIIA2)
"""
function constructTrapezoidalRule(T::Type = Float64)
    A = [0 0
         1//2 1//2]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1))
end

"""
LobattoIIIA Order 4 method
"""
function constructLobattoIIIA4(T::Type = Float64)
    A = [0 0 0
         5//24 1//3 -1//24
         1//6 2//3 1//6]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIIB Order 2 method
"""
function constructLobattoIIIB2(T::Type = Float64)
    A = [1//2 0
         1//2 0]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1))
end

"""
LobattoIIIB Order 4 method

"""
function constructLobattoIIIB4(T::Type = Float64)
    A = [1//6 -1//6 0
         1//6 1//3 0
         1//6 5//6 0]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIIC Order 2 method

"""
function constructLobattoIIIC2(T::Type = Float64)
    A = [1//2 -1//2
         1//2 1//2]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    αEEst = [1; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 2, αEEst = αEEst, adaptiveorder = 1))
end

"""
LobattoIIIC Order 4 method

"""
function constructLobattoIIIC4(T::Type = Float64)
    A = [1//6 -1//3 1//6
         1//6 5//12 -1//12
         1//6 2//3 1//6]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIIC* Order 2 method

"""
function constructLobattoIIICStar2(T::Type = Float64)
    A = [0 0
         1 0]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 2, stability_size = -2.0))
end

"""
LobattoIIIC* Order 4 method

"""
function constructLobattoIIICStar4(T::Type = Float64)
    A = [0 0 0
         1//4 1//4 0
         0 1 0]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
LobattoIIID Order 2 method

"""
function constructLobattoIIID2(T::Type = Float64)
    A = [1//2 1//2
         -1//2 1//2]
    c = [0; 1]
    α = [1 // 2; 1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 2))
end

"""
LobattoIIID Order 4 method

"""
function constructLobattoIIID4(T::Type = Float64)
    A = [1//6 0 -1//6
         1//12 5//12 0
         1//2 1//3 1//6]
    c = [0; 1 // 2; 1]
    α = [1 // 6; 2 // 3; 1 // 6]
    αEEst = [-1 // 2; 2; -1 // 2]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ImplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 2))
end

"""
RadauIA Order 3 method

"""
function constructRadauIA3(T::Type = Float64)
    A = [1//4 -1//4
         1//4 5//12]
    c = [0; 2 // 3]
    α = [1 // 4; 3 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 3))
end

"""
RadauIA Order 5 method

"""
function constructRadauIA5(T::Type = Float64)
    A = [1//9 (-1 - sqrt(6))/18 (-1 + sqrt(6))/18
         1//9 11 // 45+7 * sqrt(6) / 360 11 // 45-43 * sqrt(6) / 360
         1//9 11 // 45+43 * sqrt(6) / 360 11 // 45-7 * sqrt(6) / 360]
    c = [0; 3 // 5 - sqrt(6) / 10; 3 // 5 + sqrt(6) / 10]
    α = [1 // 9; 4 // 9 + sqrt(6) / 36; 4 // 9 - sqrt(6) / 36]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 5))
end

"""
RadauIIA Order 3 method

"""
function constructRadauIIA3(T::Type = Float64)
    A = [5//12 -1//12
         3//4 1//4]
    c = [1 // 3; 1]
    α = [3 // 4; 1 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 3))
end

"""
RadauIIA Order 5 method

"""
function constructRadauIIA5(T::Type = Float64)
    sq6 = sqrt(6)
    A = [11 // 45-7sq6 / 360 37 // 225-169sq6 / 1800 -2 // 225+sq6 / 75
         37 // 225+169sq6 / 1800 11 // 45+7sq6 / 360 -2 // 225-sq6 / 75
         4 // 9-sq6 / 36 4 // 9+sq6 / 36 1//9]
    c = [2 // 5 - sq6 / 10; 2 / 5 + sq6 / 10; 1]
    α = [4 // 9 - sq6 / 36; 4 // 9 + sq6 / 36; 1 // 9]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ImplicitRKTableau(A, c, α, 5))
end

"""
Runge-Kutta-Fehlberg Order 4/5 method.
"""
function constructRKF5(T::Type = Float64)
    A = [0 0 0 0 0 0
         1//4 0 0 0 0 0
         3//32 9//32 0 0 0 0
         1932//2197 -7200//2197 7296//2197 0 0 0
         439//216 -8 3680//513 -845//4104 0 0
         -8//27 2 -3544//2565 1859//4104 -11//40 0]
    c = [0; 1 // 4; 3 // 8; 12 // 13; 1; 1 // 2]
    α = [16 // 135; 0; 6656 // 12825; 28561 // 56430; -9 // 50; 2 // 55]
    αEEst = [25 // 216; 0; 1408 // 2565; 2197 // 4104; -1 // 5; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
        stability_size = 3.6777066213218945))
end

"""

Runge's First Order 5 method

"""
function constructRungeFirst5(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)
    c[2] = 1 // 5
    c[3] = 2 // 5
    c[4] = 1
    c[5] = 3 // 5
    c[6] = 4 // 5
    A[2, 1] = 1 // 5
    A[3, 1] = 0
    A[3, 2] = 2 // 5
    A[4, 1] = 9 // 4
    A[4, 2] = -5
    A[4, 3] = 15 // 4
    A[5, 1] = -63 // 100
    A[5, 2] = 9 // 5
    A[5, 3] = -13 // 20
    A[5, 4] = 2 // 25
    A[6, 1] = -6 // 25
    A[6, 2] = 4 // 5
    A[6, 3] = 2 // 15
    A[6, 4] = 8 // 75
    A[6, 5] = 0
    α[1] = 17 // 144
    α[2] = 0
    α[3] = 25 // 36
    α[4] = 1 // 72
    α[5] = -25 // 72
    α[6] = 25 // 48
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 3.2170478666401054))
end

"""

Cassity's Order 5 method
"""
function constructCassity5(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)
    c[2] = 1 // 7
    c[3] = 5 // 14
    c[4] = 9 // 14
    c[5] = 6 // 7
    c[6] = 1
    A[2, 1] = 1 // 7
    A[3, 1] = -367 // 4088
    A[3, 2] = 261 // 584
    A[4, 1] = 41991 // 2044
    A[4, 2] = -2493 // 73
    A[4, 3] = 57 // 4
    A[5, 1] = -108413 // 196224
    A[5, 2] = 58865 // 65408
    A[5, 3] = 5 // 16
    A[5, 4] = 265 // 1344
    A[6, 1] = -204419 // 58984
    A[6, 2] = 143829 // 58984
    A[6, 3] = 171 // 202
    A[6, 4] = 2205 // 404
    A[6, 5] = -432 // 101
    α[1] = 1 // 9
    α[2] = 7 // 2700
    α[3] = 413 // 810
    α[4] = 7 // 450
    α[5] = 28 // 75
    α[6] = -101 // 8100
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 2.1686021275264866))
end

"""

Lawson's 5th order scheme

An Order Five Runge Kutta Process with Extended Region of Stability, J. Douglas Lawson,
 Siam Journal on Numerical Analysis, Vol. 3, No. 4, (Dec., 1966) pages 593-597

"""
function constructLawson5(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)
    c[2] = 1 // 12
    c[3] = 1 // 4
    c[4] = 1 // 2
    c[5] = 3 // 4
    c[6] = 1
    A[2, 1] = 1 // 12
    A[3, 1] = -1 // 8
    A[3, 2] = 3 // 8
    A[4, 1] = 3 // 5
    A[4, 2] = -9 // 10
    A[4, 3] = 4 // 5
    A[5, 1] = 39 // 80
    A[5, 2] = -9 // 20
    A[5, 3] = 3 // 20
    A[5, 4] = 9 // 16
    A[6, 1] = -59 // 35
    A[6, 2] = 66 // 35
    A[6, 3] = 48 // 35
    A[6, 4] = -12 // 7
    A[6, 5] = 8 // 7
    α[1] = 7 // 90
    α[2] = 0
    α[3] = 16 // 45
    α[4] = 2 // 15
    α[5] = 16 // 45
    α[6] = 7 // 90
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 3.7343596072347225))
end

"""

Luther and Konen's First Order 5
Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,
 Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

"""
function constructLutherKonen5(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)
    c[2] = 1 / 2
    c[3] = 1 / 2 - 1 / 10 * 5^(1 / 2)
    c[4] = 1 / 2
    c[5] = 1 / 2 + 1 / 10 * 5^(1 / 2)
    c[6] = 1
    A[2, 1] = 1 / 2
    A[3, 1] = 1 / 5
    A[3, 2] = 3 / 10 - 1 / 10 * 5^(1 / 2)
    A[4, 1] = 1 / 4
    A[4, 2] = 1 / 4
    A[4, 3] = 0
    A[5, 1] = 1 / 20 - 1 / 20 * 5^(1 / 2)
    A[5, 2] = -1 / 5
    A[5, 3] = 1 / 4 + 3 / 20 * 5^(1 / 2)
    A[5, 4] = 2 / 5
    A[6, 1] = 1 / 4 * 5^(1 / 2) - 1 / 4
    A[6, 2] = 1 / 2 * 5^(1 / 2) - 1 / 2
    A[6, 3] = 5 / 4 - 1 / 4 * 5^(1 / 2)
    A[6, 4] = -2
    A[6, 5] = 5 / 2 - 1 / 2 * 5^(1 / 2)
    α[1] = 1 / 12
    α[2] = 0
    α[3] = 5 / 12
    α[4] = 0
    α[5] = 5 / 12
    α[6] = 1 / 12
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 3.217047866640106))
end

"""

Luther and Konen's Second Order 5
Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,
 Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

"""
function constructLutherKonen52(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)
    c[2] = 2 / 5
    c[3] = 1 / 2
    c[4] = 1
    c[5] = 1 / 2 - 1 / 10 * 15^(1 / 2)
    c[6] = 1 / 2 + 1 / 10 * 15^(1 / 2)
    A[2, 1] = 2 / 5
    A[3, 1] = 3 / 16
    A[3, 2] = 5 / 16
    A[4, 1] = 1 / 4
    A[4, 2] = -5 / 4
    A[4, 3] = 2
    A[5, 1] = 3 / 20 - 1 / 100 * 15^(1 / 2)
    A[5, 2] = -1 / 4
    A[5, 3] = 3 / 5 - 2 / 25 * 15^(1 / 2)
    A[5, 4] = -1 / 100 * 15^(1 / 2)
    A[6, 1] = -3 / 20 - 1 / 20 * 15^(1 / 2)
    A[6, 2] = -1 / 4
    A[6, 3] = 3 / 5
    A[6, 4] = 3 / 10 - 1 / 20 * 15^(1 / 2)
    A[6, 5] = 1 / 5 * 15^(1 / 2)
    α[1] = 0
    α[2] = 0
    α[3] = 4 / 9
    α[4] = 0
    α[5] = 5 / 18
    α[6] = 5 / 18
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 2.6515956444339794))
end

"""

Luther and Konen's Third Order 5
Some Fifth-Order Classical Runge Kutta Formulas, H.A.Luther and H.P.Konen,
 Siam Review, Vol. 3, No. 7, (Oct., 1965) pages 551-558.

"""
function constructLutherKonen53(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)
    c[2] = 3 / 25
    c[3] = 5 / 18
    c[4] = 45 / 89
    c[5] = 3 / 5 - 1 / 10 * 6^(1 / 2)
    c[6] = 3 / 5 + 1 / 10 * 6^(1 / 2)
    A[2, 1] = 3 / 25
    A[3, 1] = -85 / 1944
    A[3, 2] = 625 / 1944
    A[4, 1] = 610425 / 1409938
    A[4, 2] = -961875 / 1409938
    A[4, 3] = 532170 / 704969
    A[5, 1] = 7411 / 37500 - 673 / 18750 * 6^(1 / 2)
    A[5, 2] = 0
    A[5, 3] = 27621 / 228125 * 6^(1 / 2) - 6561 / 228125
    A[5, 4] = 1180229 / 2737500 - 126736 / 684375 * 6^(1 / 2)
    A[6, 1] = -5351 / 62500 - 7087 / 281250 * 6^(1 / 2)
    A[6, 2] = 0
    A[6, 3] = 2736423 / 1140625 + 786753 / 1140625 * 6^(1 / 2)
    A[6, 4] = 73736589 / 86687500 + 101816534 / 195046875 * 6^(1 / 2)
    A[6, 5] = -30448 / 11875 - 12903 / 11875 * 6^(1 / 2)
    α[1] = 1 / 9
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 4 / 9 + 1 / 36 * 6^(1 / 2)
    α[6] = 4 / 9 - 1 / 36 * 6^(1 / 2)
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 3.679935798458283))
end

"""

 S.N. Papakostas and G. PapaGeorgiou higher error more stable

 A Family of Fifth-order Runge-Kutta Pairs, by S.N. Papakostas and G. PapaGeorgiou,
 Mathematics of Computation,Volume 65, Number 215, July 1996, Pages 1165-1181.

"""
function constructPapakostasPapaGeorgiou5(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)
    c[2] = 64 // 315
    c[3] = 115 // 381
    c[4] = 762 // 935
    c[5] = 25 // 28
    c[6] = 1
    c[7] = 1
    A[2, 1] = 64 // 315
    A[3, 1] = 480815 // 6193536
    A[3, 2] = 462875 // 2064512
    A[4, 1] = 344904825219069 // 345923838700000
    A[4, 2] = -2360077908867 // 601606676000
    A[4, 3] = 40439332863108 // 10810119959375
    A[5, 1] = 12078745127989699 // 5009699157786624
    A[5, 2] = -791781731775 // 81669668864
    A[5, 3] = 39297175833216951 // 4694461413969824
    A[5, 4] = -10508413393960625 // 54097233987826176
    A[6, 1] = 2251421737440863 // 828701767536000
    A[6, 2] = -39895842357 // 3782730880
    A[6, 3] = 34564628685305112534 // 3916944972468643375
    A[6, 4] = 12051135495733565 // 36492943960723917
    A[6, 5] = -26808346215168 // 82592376030125
    A[7, 1] = 2405713 // 26289000
    A[7, 2] = 0
    A[7, 3] = 63896466577779 // 141024193000600
    A[7, 4] = 454128848141375 // 589615117674696
    A[7, 5] = -1359311744 // 2892576375
    A[7, 6] = 256979 // 1656648
    α[1] = 2405713 // 26289000
    α[2] = 0
    α[3] = 63896466577779 // 141024193000600
    α[4] = 454128848141375 // 589615117674696
    α[5] = -1359311744 // 2892576375
    α[6] = 256979 // 1656648
    αEEst[1] = 1818563883019 // 20194131951000
    αEEst[2] = 0
    αEEst[3] = 5513862498202899713 // 12036555896794210600
    αEEst[4] = 324806515311046773125 // 452918159177876804664
    αEEst[5] = -126112324722496 // 317422653663375
    αEEst[6] = 137695258717 // 1272569071032
    αEEst[7] = 1 // 42
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4, fsal = true,
        stability_size = 3.306567892634947))
end

"""

 S.N. Papakostas and G. PapaGeorgiou less stable lower error
 Strictly better than DP5

 A Family of Fifth-order Runge-Kutta Pairs, by S.N. Papakostas and G. PapaGeorgiou,
 Mathematics of Computation,Volume 65, Number 215, July 1996, Pages 1165-1181.

"""
function constructPapakostasPapaGeorgiou52(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)

    c[2] = 35 // 159
    c[3] = 42 // 131
    c[4] = 131 // 143
    c[5] = 21 // 22
    c[6] = 1
    c[7] = 1
    A[2, 1] = 35 // 159
    A[3, 1] = 7476 // 85805
    A[3, 2] = 20034 // 85805
    A[4, 1] = 2438549411 // 1983961980
    A[4, 2] = -3707256508 // 716430715
    A[4, 3] = 25077455105 // 5158301148
    A[5, 1] = 105337889067 // 64388030080
    A[5, 2] = -1698584121 // 245755840
    A[5, 3] = 6869523776931 // 1096562558080
    A[5, 4] = -927215289 // 26981535520
    A[6, 1] = 67512025387 // 32454592380
    A[6, 2] = -20051384 // 2293935
    A[6, 3] = 10587214001321 // 1373901639516
    A[6, 4] = 731293420 // 8209319229
    A[6, 5] = -144610048 // 1077690663
    A[7, 1] = 669707 // 6932520
    A[7, 2] = 0
    A[7, 3] = 2215522905683 // 4570867891800
    A[7, 4] = 349043981 // 116904400
    A[7, 5] = -2234144 // 575505
    A[7, 6] = 9363 // 7120
    α[1] = 669707 // 6932520
    α[2] = 0
    α[3] = 2215522905683 // 4570867891800
    α[4] = 349043981 // 116904400
    α[5] = -2234144 // 575505
    α[6] = 9363 // 7120
    αEEst[1] = 2243660497 // 23535905400
    αEEst[2] = 0
    αEEst[3] = 7589131232781673 // 15518096492661000
    αEEst[4] = 1104461697911 // 396890438000
    αEEst[5] = -6925033984 // 1953839475
    αEEst[6] = 3529851 // 3021550
    αEEst[7] = 1 // 112

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4, fsal = true,
        stability_size = 3.4272251630453394))
end

"""
Runge–Kutta pairs of orders 5(4) using the minimal set of simplifying assumptions,
 by Ch. Tsitouras, TEI of Chalkis, Dept. of Applied Sciences, GR34400, Psahna, Greece.
"""
function constructTsitouras5(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)
    c[2] = convert(T, 161 // 1000)
    c[3] = convert(T, 327 // 1000)
    c[4] = convert(T, 9 // 10)
    c[5] = convert(T,
        big".9800255409045096857298102862870245954942137979563024768854764293221195950761080302604")
    c[6] = convert(T, 1)
    c[7] = convert(T, 1)
    A[2, 1] = convert(T, 161 // 1000)
    A[3, 1] = convert(T,
        big"-.8480655492356988544426874250230774675121177393430391537369234245294192976164141156943e-2")
    A[3, 2] = convert(T,
        big".3354806554923569885444268742502307746751211773934303915373692342452941929761641411569")
    A[4, 1] = convert(T,
        big"2.897153057105493432130432594192938764924887287701866490314866693455023795137503079289")
    A[4, 2] = convert(T,
        big"-6.359448489975074843148159912383825625952700647415626703305928850207288721235210244366")
    A[4, 3] = convert(T,
        big"4.362295432869581411017727318190886861027813359713760212991062156752264926097707165077")
    A[5, 1] = convert(T,
        big"5.325864828439256604428877920840511317836476253097040101202360397727981648835607691791")
    A[5, 2] = convert(T,
        big"-11.74888356406282787774717033978577296188744178259862899288666928009020615663593781589")
    A[5, 3] = convert(T,
        big"7.495539342889836208304604784564358155658679161518186721010132816213648793440552049753")
    A[5, 4] = convert(T,
        big"-.9249506636175524925650207933207191611349983406029535244034750452930469056411389539635e-1")
    A[6, 1] = convert(T,
        big"5.861455442946420028659251486982647890394337666164814434818157239052507339770711679748")
    A[6, 2] = convert(T,
        big"-12.92096931784710929170611868178335939541780751955743459166312250439928519268343184452")
    A[6, 3] = convert(T,
        big"8.159367898576158643180400794539253485181918321135053305748355423955009222648673734986")
    A[6, 4] = convert(T,
        big"-.7158497328140099722453054252582973869127213147363544882721139659546372402303777878835e-1")
    A[6, 5] = convert(T,
        big"-.2826905039406838290900305721271224146717633626879770007617876201276764571291579142206e-1")
    A[7, 1] = convert(T,
        big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1")
    A[7, 2] = convert(T, 1 // 100)
    A[7, 3] = convert(T,
        big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509")
    A[7, 4] = convert(T,
        big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331")
    A[7, 5] = convert(T,
        big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677")
    A[7, 6] = convert(T,
        big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841")
    α[1] = convert(T,
        big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1")
    α[2] = convert(T, 1 // 100)
    α[3] = convert(T,
        big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509")
    α[4] = convert(T,
        big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331")
    α[5] = convert(T,
        big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677")
    α[6] = convert(T,
        big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841")
    αEEst[1] = convert(T,
        big".9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1")
    αEEst[2] = convert(T,
        big".9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2")
    αEEst[3] = convert(T,
        big".4877705284247615707855642599631228241516691959761363774365216240304071651579571959813")
    αEEst[4] = convert(T,
        big"1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761")
    αEEst[5] = convert(T,
        big"-2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702")
    αEEst[6] = convert(T,
        big"1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089")
    αEEst[7] = convert(T, 1 // 66)
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4, fsal = true,
        stability_size = 3.5068469938049547))
end

"""

An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28
"""
function constructBogakiShampine5(T::Type = Float64)
    A = zeros(T, 8, 8)
    c = zeros(T, 8)
    α = zeros(T, 8)
    αEEst = zeros(T, 8)
    αEEst2 = zeros(T, 8)
    c[2] = 1 // 6
    c[3] = 2 // 9
    c[4] = 3 // 7
    c[5] = 2 // 3
    c[6] = 3 // 4
    c[7] = 1
    c[8] = 1
    A[2, 1] = 1 // 6
    A[3, 1] = 2 // 27
    A[3, 2] = 4 // 27
    A[4, 1] = 183 // 1372
    A[4, 2] = -162 // 343
    A[4, 3] = 1053 // 1372
    A[5, 1] = 68 // 297
    A[5, 2] = -4 // 11
    A[5, 3] = 42 // 143
    A[5, 4] = 1960 // 3861
    A[6, 1] = 597 // 22528
    A[6, 2] = 81 // 352
    A[6, 3] = 63099 // 585728
    A[6, 4] = 58653 // 366080
    A[6, 5] = 4617 // 20480
    A[7, 1] = 174197 // 959244
    A[7, 2] = -30942 // 79937
    A[7, 3] = 8152137 // 19744439
    A[7, 4] = 666106 // 1039181
    A[7, 5] = -29421 // 29068
    A[7, 6] = 482048 // 414219
    A[8, 1] = 587 // 8064
    A[8, 2] = 0
    A[8, 3] = 4440339 // 15491840
    A[8, 4] = 24353 // 124800
    A[8, 5] = 387 // 44800
    A[8, 6] = 2152 // 5985
    A[8, 7] = 7267 // 94080
    α[1] = 587 // 8064
    α[2] = 0
    α[3] = 4440339 // 15491840
    α[4] = 24353 // 124800
    α[5] = 387 // 44800
    α[6] = 2152 // 5985
    α[7] = 7267 // 94080
    α[8] = 0
    αEEst[1] = 6059 // 80640
    αEEst[2] = 0
    αEEst[3] = 8559189 // 30983680
    αEEst[4] = 26411 // 124800
    αEEst[5] = -927 // 89600
    αEEst[6] = 443 // 1197
    αEEst[7] = 7267 // 94080
    αEEst2[1] = 2479 // 34992
    αEEst2[2] = 0
    αEEst2[3] = 123 // 416
    αEEst2[4] = 612941 // 3411720
    αEEst2[5] = 43 // 1440
    αEEst2[6] = 2272 // 6561
    αEEst2[7] = 79937 // 1113912
    αEEst2[8] = 3293 // 556956
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
        stability_size = 3.9879271987261333))
end

"""
Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum, by P.W.Sharp and E.Smart,
 Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.

"""
function constructSharpSmart5(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)
    αEEst2 = zeros(T, 7)
    c[2] = 16 // 105
    c[3] = 8 // 35
    c[4] = 9 // 20
    c[5] = 2 // 3
    c[6] = 7 // 9
    c[7] = 1
    A[2, 1] = 16 // 105
    A[3, 1] = 2 // 35
    A[3, 2] = 6 // 35
    A[4, 1] = 8793 // 40960
    A[4, 2] = -5103 // 8192
    A[4, 3] = 17577 // 20480
    A[5, 1] = 347 // 1458
    A[5, 2] = -7 // 20
    A[5, 3] = 3395 // 10044
    A[5, 4] = 49792 // 112995
    A[6, 1] = -1223224109959 // 9199771214400
    A[6, 2] = 1234787701 // 2523942720
    A[6, 3] = 568994101921 // 3168810084960
    A[6, 4] = -105209683888 // 891227836395
    A[6, 5] = 9 // 25
    A[7, 1] = 2462504862877 // 8306031988800
    A[7, 2] = -123991 // 287040
    A[7, 3] = 106522578491 // 408709510560
    A[7, 4] = 590616498832 // 804646848915
    A[7, 5] = -319138726 // 534081275
    A[7, 6] = 52758 // 71449
    α[1] = 1093 // 15120
    α[2] = 0
    α[3] = 60025 // 190992
    α[4] = 3200 // 20709
    α[5] = 1611 // 11960
    α[6] = 712233 // 2857960
    α[7] = 3 // 40
    αEEst[1] = 84018211 // 991368000
    αEEst[2] = 0
    αEEst[3] = 92098979 // 357791680
    αEEst[4] = 17606944 // 67891005
    αEEst[5] = 3142101 // 235253200
    αEEst[6] = 22004596809 // 70270091500
    αEEst[7] = 9 // 125
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
        stability_size = 3.9156746135081772))
end

"""
constructBogakiShampine3()

Constructs the tableau object for the Bogakai-Shampine Order 2/3 method.
"""
function constructBogakiShampine3(T::Type = Float64)
    A = [0 0 0 0
         1//2 0 0 0
         0 3//4 0 0
         2//9 1//3 4//9 0]
    c = [0; 1 // 2; 3 // 4; 1]
    α = [2 // 9; 1 // 3; 4 // 9; 0]
    αEEst = [7 // 24; 1 // 4; 1 // 3; 1 // 8]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 3, αEEst = αEEst, adaptiveorder = 2,
        stability_size = 2.5127453266183286))
end

"""
constructCashKarp()

Constructs the tableau object for the Cash-Karp Order 4/5 method.
"""
function constructCashKarp(T::Type = Float64)
    A = [0 0 0 0 0 0
         1//5 0 0 0 0 0
         3//40 9//40 0 0 0 0
         3//10 -9//10 6//5 0 0 0
         -11//54 5//2 -70//27 35//27 0 0
         1631//55296 175//512 575//13824 44275//110592 253//4096 0]
    c = [0; 1 // 5; 3 // 10; 3 // 5; 1; 7 // 8]
    α = [37 // 378; 0; 250 // 621; 125 // 594; 0; 512 // 1771]
    αEEst = [2825 // 27648; 0; 18575 // 48384; 13525 // 55296; 277 // 14336; 1 // 4]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
        stability_size = 3.7343596072347225))
end

"""
Runge-Kutta-Fehberg Order 4/3
"""
function constructRKF4(T::Type = Float64)
    c = [0; 1 // 4; 4 // 9; 6 // 7; 1]
    A = [0 0 0 0 0
         1//4 0 0 0 0
         4//81 32//81 0 0 0
         57//98 -432//343 1053//686 0 0
         1//6 0 27//52 49//156 0]
    α = [43 // 288; 0; 243 // 416; 343 // 1872; 1 // 12]
    αEEst = [1 // 6; 0; 27 // 52; 49 // 156; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 4, αEEst = αEEst, adaptiveorder = 3,
        stability_size = 4.109222736949077))
end

"""
Butcher's Third Order 6

On Runge-Kutta Processes of High Order, by J. C. Butcher,
 Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194

"""
function constructButcher63(T::Type = Float64)
    c = [0; 1 / 2; 2 / 3; 1 / 3; 5 / 6; 1 / 6; 1]
    A = [0 0 0 0 0 0 0
         1//2 0 0 0 0 0 0
         2//9 4//9 0 0 0 0 0
         7//36 2//9 -1//12 0 0 0 0
         -35//144 -55//36 35//48 15//8 0 0 0
         -1//360 -11//36 -1//8 1//2 1//10 0 0
         -41//260 22//13 43//156 -118//39 32//195 80//39 0]
    α = [13 // 200; 0; 11 // 40; 11 // 40; 4 // 25; 4 // 25; 13 // 200]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 2.856108978668386))
end

"""

Butcher's First Order 6 method

On Runge-Kutta Processes of High Order, by J. C. Butcher,
 Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194
"""
function constructButcher6(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)

    c[2] = 1 // 2 - 1 // 10 * 5^(1 // 2)
    c[3] = 1 // 2 + 1 // 10 * 5^(1 // 2)
    c[4] = 1 // 2 - 1 // 10 * 5^(1 // 2)
    c[5] = 1 // 2 + 1 // 10 * 5^(1 // 2)
    c[6] = 1 // 2 - 1 // 10 * 5^(1 // 2)
    c[7] = 1
    A[2, 1] = 1 // 2 - 1 // 10 * 5^(1 // 2)
    A[3, 1] = -1 // 10 * 5^(1 // 2)
    A[3, 2] = 1 // 2 + 1 // 5 * 5^(1 // 2)
    A[4, 1] = 7 // 20 * 5^(1 // 2) - 3 // 4
    A[4, 2] = 1 // 4 * 5^(1 // 2) - 1 // 4
    A[4, 3] = 3 // 2 - 7 // 10 * 5^(1 // 2)
    A[5, 1] = 1 // 12 - 1 // 60 * 5^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = 1 // 6
    A[5, 4] = 7 // 60 * 5^(1 // 2) + 1 // 4
    A[6, 1] = 1 // 60 * 5^(1 // 2) + 1 // 12
    A[6, 2] = 0
    A[6, 3] = 3 // 4 - 5 // 12 * 5^(1 // 2)
    A[6, 4] = 1 // 6
    A[6, 5] = -1 // 2 + 3 // 10 * 5^(1 // 2)
    A[7, 1] = 1 // 6
    A[7, 2] = 0
    A[7, 3] = -55 // 12 + 25 // 12 * 5^(1 // 2)
    A[7, 4] = -7 // 12 * 5^(1 // 2) - 25 // 12
    A[7, 5] = 5 - 2 * 5^(1 // 2)
    A[7, 6] = 5 // 2 + 1 // 2 * 5^(1 // 2)
    α[1] = 1 // 12
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 5 // 12
    α[6] = 5 // 12
    α[7] = 1 // 12

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 2.865620537276994))
end

"""

Butcher's Second Order 6 method

On Runge-Kutta Processes of High Order, by J. C. Butcher,
 Journal of the Australian Mathematical Society, Vol. 4, (1964), pages 179 to 194
"""
function constructButcher62(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)

    c[2] = 1 // 3
    c[3] = 2 // 3
    c[4] = 1 // 3
    c[5] = 1 // 2
    c[6] = 1 // 2
    c[7] = 1
    A[2, 1] = 1 // 3
    A[3, 1] = 0
    A[3, 2] = 2 // 3
    A[4, 1] = 1 // 12
    A[4, 2] = 1 // 3
    A[4, 3] = -1 // 12
    A[5, 1] = -1 // 16
    A[5, 2] = 9 // 8
    A[5, 3] = -3 // 16
    A[5, 4] = -3 // 8
    A[6, 1] = 0
    A[6, 2] = 9 // 8
    A[6, 3] = -3 // 8
    A[6, 4] = -3 // 4
    A[6, 5] = 1 // 2
    A[7, 1] = 9 // 44
    A[7, 2] = -9 // 11
    A[7, 3] = 63 // 44
    A[7, 4] = 18 // 11
    A[7, 5] = 0
    A[7, 6] = -16 // 11
    α[1] = 11 // 120
    α[2] = 0
    α[3] = 27 // 40
    α[4] = 27 // 40
    α[5] = -4 // 15
    α[6] = -4 // 15
    α[7] = 11 // 120

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 2.8561089786683866))
end

"""

Verner Order 5/6 method

A Contrast of a New RK56 pair with DP56, by Jim Verner,
 Department of Mathematics. Simon Fraser University, Burnaby, Canada, 2006.

"""
function constructVerner6(T::Type = Float64)
    A = zeros(T, 8, 8)
    c = zeros(T, 8)
    α = zeros(T, 8)
    αEEst = zeros(T, 8)

    c[2] = 1 // 7
    c[3] = 2 // 9
    c[4] = 3 // 7
    c[5] = 2 // 3
    c[6] = 3 // 4
    c[7] = 1
    c[8] = 1
    A[2, 1] = 1 // 7
    A[3, 1] = 4 // 81
    A[3, 2] = 14 // 81
    A[4, 1] = 291 // 1372
    A[4, 2] = -27 // 49
    A[4, 3] = 1053 // 1372
    A[5, 1] = 86 // 297
    A[5, 2] = -14 // 33
    A[5, 3] = 42 // 143
    A[5, 4] = 1960 // 3861
    A[6, 1] = -267 // 22528
    A[6, 2] = 189 // 704
    A[6, 3] = 63099 // 585728
    A[6, 4] = 58653 // 366080
    A[6, 5] = 4617 // 20480
    A[7, 1] = 10949 // 6912
    A[7, 2] = -69 // 32
    A[7, 3] = -90891 // 68096
    A[7, 4] = 112931 // 25920
    A[7, 5] = -69861 // 17920
    A[7, 6] = 26378 // 10773
    A[8, 1] = 1501 // 19008
    A[8, 2] = -21 // 88
    A[8, 3] = 219519 // 347776
    A[8, 4] = 163807 // 926640
    A[8, 5] = -417 // 640
    A[8, 6] = 1544 // 1539
    A[8, 7] = 0
    α[1] = 79 // 1080
    α[2] = 0
    α[3] = 19683 // 69160
    α[4] = 16807 // 84240
    α[5] = 0
    α[6] = 2816 // 7695
    α[7] = 1 // 100
    α[8] = 187 // 2800
    αEEst[1] = 763 // 10800
    αEEst[2] = 0
    αEEst[3] = 59049 // 197600
    αEEst[4] = 88837 // 526500
    αEEst[5] = 243 // 4000
    αEEst[6] = 12352 // 38475
    αEEst[7] = 0
    αEEst[8] = 2 // 25

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst,
        stability_size = 3.9541297306311853))
end

"""

Dormand-Prince Order 5//6 method

P.J. Prince and J. R. Dormand, High order embedded Runge-Kutta formulae,
Journal of Computational and Applied Mathematics . 7 (1981), pp. 67-75.
"""
function constructDormandPrince6(T::Type = Float64)
    A = zeros(T, 8, 8)
    c = zeros(T, 8)
    α = zeros(T, 8)
    αEEst = zeros(T, 8)

    c[2] = 1 // 10
    c[3] = 2 // 9
    c[4] = 3 // 7
    c[5] = 3 // 5
    c[6] = 4 // 5
    c[7] = 1
    c[8] = 1
    A[2, 1] = 1 // 10
    A[3, 1] = -2 // 81
    A[3, 2] = 20 // 81
    A[4, 1] = 615 // 1372
    A[4, 2] = -270 // 343
    A[4, 3] = 1053 // 1372
    A[5, 1] = 3243 // 5500
    A[5, 2] = -54 // 55
    A[5, 3] = 50949 // 71500
    A[5, 4] = 4998 // 17875
    A[6, 1] = -26492 // 37125
    A[6, 2] = 72 // 55
    A[6, 3] = 2808 // 23375
    A[6, 4] = -24206 // 37125
    A[6, 5] = 338 // 459
    A[7, 1] = 5561 // 2376
    A[7, 2] = -35 // 11
    A[7, 3] = -24117 // 31603
    A[7, 4] = 899983 // 200772
    A[7, 5] = -5225 // 1836
    A[7, 6] = 3925 // 4056
    A[8, 1] = 465467 // 266112
    A[8, 2] = -2945 // 1232
    A[8, 3] = -5610201 // 14158144
    A[8, 4] = 10513573 // 3212352
    A[8, 5] = -424325 // 205632
    A[8, 6] = 376225 // 454272
    A[8, 7] = 0
    α[1] = 61 // 864
    α[2] = 0
    α[3] = 98415 // 321776
    α[4] = 16807 // 146016
    α[5] = 1375 // 7344
    α[6] = 1375 // 5408
    α[7] = -37 // 1120
    α[8] = 1 // 10
    αEEst[1] = 821 // 10800
    αEEst[2] = 0
    αEEst[3] = 19683 // 71825
    αEEst[4] = 175273 // 912600
    αEEst[5] = 395 // 3672
    αEEst[6] = 785 // 2704
    αEEst[7] = 3 // 50
    αEEst[8] = 0

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst,
        stability_size = 3.954129730631186))
end

"""

Sharp-Verner Order 5/6 method

Completely Imbedded Runge-Kutta Pairs, by P. W. Sharp and J. H. Verner,
 SIAM Journal on Numerical Analysis, Vol. 31, No. 4. (Aug., 1994), pages. 1169 to 1190.
"""
function constructSharpVerner6(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 12
    c[3] = 2 // 15
    c[4] = 1 // 5
    c[5] = 8 // 15
    c[6] = 2 // 3
    c[7] = 19 // 20
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 12
    A[3, 1] = 2 // 75
    A[3, 2] = 8 // 75
    A[4, 1] = 1 // 20
    A[4, 2] = 0
    A[4, 3] = 3 // 20
    A[5, 1] = 88 // 135
    A[5, 2] = 0
    A[5, 3] = -112 // 45
    A[5, 4] = 64 // 27
    A[6, 1] = -10891 // 11556
    A[6, 2] = 0
    A[6, 3] = 3880 // 963
    A[6, 4] = -8456 // 2889
    A[6, 5] = 217 // 428
    A[7, 1] = 1718911 // 4382720
    A[7, 2] = 0
    A[7, 3] = -1000749 // 547840
    A[7, 4] = 819261 // 383488
    A[7, 5] = -671175 // 876544
    A[7, 6] = 14535 // 14336
    A[8, 1] = 85153 // 203300
    A[8, 2] = 0
    A[8, 3] = -6783 // 2140
    A[8, 4] = 10956 // 2675
    A[8, 5] = -38493 // 13375
    A[8, 6] = 1152 // 425
    A[8, 7] = -7168 // 40375
    A[9, 1] = 53 // 912
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 5 // 16
    A[9, 5] = 27 // 112
    A[9, 6] = 27 // 136
    A[9, 7] = 256 // 969
    A[9, 8] = -25 // 336
    α[1] = 53 // 912
    α[2] = 0
    α[3] = 0
    α[4] = 5 // 16
    α[5] = 27 // 112
    α[6] = 27 // 136
    α[7] = 256 // 969
    α[8] = -25 // 336
    αEEst[1] = 617 // 10944
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 241 // 756
    αEEst[5] = 69 // 320
    αEEst[6] = 435 // 1904
    αEEst[7] = 10304 // 43605
    αEEst[8] = 0
    αEEst[9] = -1 // 18

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst, fsal = true,
        stability_size = 4.470828492062036))
end

"""

Verner 1991 Second Order 5/6 method

Some Ruge-Kutta Formula Pairs, by J.H.Verner,
 SIAM Journal on Numerical Analysis, Vol. 28, No. 2 (April 1991), pages 496 to 511.
"""
function constructVerner9162(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 8
    c[3] = 1 // 6
    c[4] = 1 // 4
    c[5] = 1 // 2
    c[6] = 3 // 5
    c[7] = 4 // 5
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 8
    A[3, 1] = 1 // 18
    A[3, 2] = 1 // 9
    A[4, 1] = 1 // 16
    A[4, 2] = 0
    A[4, 3] = 3 // 16
    A[5, 1] = 1 // 4
    A[5, 2] = 0
    A[5, 3] = -3 // 4
    A[5, 4] = 1
    A[6, 1] = 134 // 625
    A[6, 2] = 0
    A[6, 3] = -333 // 625
    A[6, 4] = 476 // 625
    A[6, 5] = 98 // 625
    A[7, 1] = -98 // 1875
    A[7, 2] = 0
    A[7, 3] = 12 // 625
    A[7, 4] = 10736 // 13125
    A[7, 5] = -1936 // 1875
    A[7, 6] = 22 // 21
    A[8, 1] = 9 // 50
    A[8, 2] = 0
    A[8, 3] = 21 // 25
    A[8, 4] = -2924 // 1925
    A[8, 5] = 74 // 25
    A[8, 6] = -15 // 7
    A[8, 7] = 15 // 22
    A[9, 1] = 11 // 144
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 256 // 693
    A[9, 5] = 0
    A[9, 6] = 125 // 504
    A[9, 7] = 125 // 528
    A[9, 8] = 5 // 72
    α[1] = 11 // 144
    α[2] = 0
    α[3] = 0
    α[4] = 256 // 693
    α[5] = 0
    α[6] = 125 // 504
    α[7] = 125 // 528
    α[8] = 5 // 72
    αEEst[1] = 1 // 18
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 32 // 63
    αEEst[5] = -2 // 3
    αEEst[6] = 125 // 126
    αEEst[7] = 0
    αEEst[8] = -5 // 63
    αEEst[9] = 4 // 21

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst,
        stability_size = 4.239991338165166))
end

"""

Verner 1991 First Order 5/6 method

Some Ruge-Kutta Formula Pairs, by J.H.Verner,
 SIAM Journal on Numerical Analysis, Vol. 28, No. 2 (April 1991), pages 496 to 511.
"""
function constructVerner916(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 8
    c[3] = 4 // 9 - 4 // 45 * 10^(1 // 2)
    c[4] = 2 // 3 - 2 // 15 * 10^(1 // 2)
    c[5] = 9 // 16
    c[6] = 1 // 2
    c[7] = 9 // 10
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 8
    A[3, 1] = -268 // 405 + 92 // 405 * 10^(1 // 2)
    A[3, 2] = 448 // 405 - 128 // 405 * 10^(1 // 2)
    A[4, 1] = 1 // 6 - 1 // 30 * 10^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = 1 // 2 - 1 // 10 * 10^(1 // 2)
    A[5, 1] = 11547 // 32768 + 405 // 16384 * 10^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -18225 // 32768 - 5103 // 16384 * 10^(1 // 2)
    A[5, 4] = 12555 // 16384 + 2349 // 8192 * 10^(1 // 2)
    A[6, 1] = 19662371 // 51149376 + 441281 // 12787344 * 10^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = -3786045 // 5683264 - 252663 // 710408 * 10^(1 // 2)
    A[6, 4] = 1570556745 // 1821486112 + 290041461 // 910743056 * 10^(1 // 2)
    A[6, 5] = -41227072 // 512292969 + 1374464 // 512292969 * 10^(1 // 2)
    A[7, 1] = -154207593 // 369412160 - 1829424339 // 11544130000 * 10^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 2659895739 // 1847060800 + 653855409 // 1154413000 * 10^(1 // 2)
    A[7, 4] = -349492176711 // 591982986400 - 359784638379 // 1479957466000 * 10^(1 // 2)
    A[7, 5] = 153920585664 // 92497341625 + 311066673408 // 462486708125 * 10^(1 // 2)
    A[7, 6] = -1944 // 1625 - 6804 // 8125 * 10^(1 // 2)
    A[8, 1] = 70594945601 // 21406013856 + 21473424323 // 21406013856 * 10^(1 // 2)
    A[8, 2] = 0
    A[8, 3] = -794525145 // 88090592 - 249156075 // 88090592 * 10^(1 // 2)
    A[8, 4] = 866290968775 // 254097312624 + 256998959765 // 254097312624 * 10^(1 // 2)
    A[8, 5] = -15964196472448 // 1286367645159 -
              5039429245312 // 1286367645159 * 10^(1 // 2)
    A[8, 6] = 17017 // 1116 + 5075 // 1116 * 10^(1 // 2)
    A[8, 7] = 42875 // 90396 + 16625 // 90396 * 10^(1 // 2)
    A[9, 1] = 31 // 324 - 37 // 4860 * 10^(1 // 2)
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 37435 // 69228 - 3235 // 69228 * 10^(1 // 2)
    A[9, 5] = -1245184 // 1090341 + 9699328 // 16355115 * 10^(1 // 2)
    A[9, 6] = 71 // 54 - 74 // 135 * 10^(1 // 2)
    A[9, 7] = 625 // 486 - 250 // 729 * 10^(1 // 2)
    A[9, 8] = -23 // 21 + 37 // 105 * 10^(1 // 2)
    α[1] = 31 // 324 - 37 // 4860 * 10^(1 // 2)
    α[2] = 0
    α[3] = 0
    α[4] = 37435 // 69228 - 3235 // 69228 * 10^(1 // 2)
    α[5] = -1245184 // 1090341 + 9699328 // 16355115 * 10^(1 // 2)
    α[6] = 71 // 54 - 74 // 135 * 10^(1 // 2)
    α[7] = 625 // 486 - 250 // 729 * 10^(1 // 2)
    α[8] = -23 // 21 + 37 // 105 * 10^(1 // 2)
    αEEst[1] = 5 // 54 - 2 // 135 * 10^(1 // 2)
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 2390 // 17307 + 2290 // 17307 * 10^(1 // 2)
    αEEst[5] = 40960 // 121149 + 262144 // 605745 * 10^(1 // 2)
    αEEst[6] = 2 // 27 - 64 // 135 * 10^(1 // 2)
    αEEst[7] = 0
    αEEst[8] = 150029 // 443709 - 236267 // 2218545 * 10^(1 // 2)
    αEEst[9] = 2411 // 126774 + 1921 // 63387 * 10^(1 // 2)

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst,
        stability_size = 4.250607671658528))
end

"""

From Verner's Website
"""
function constructVernerRobust6(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 9 // 50
    c[3] = 1 // 6
    c[4] = 1 // 4
    c[5] = 53 // 100
    c[6] = 3 // 5
    c[7] = 4 // 5
    c[8] = 1
    c[9] = 1
    A[2, 1] = 9 // 50
    A[3, 1] = 29 // 324
    A[3, 2] = 25 // 324
    A[4, 1] = 1 // 16
    A[4, 2] = 0
    A[4, 3] = 3 // 16
    A[5, 1] = 79129 // 250000
    A[5, 2] = 0
    A[5, 3] = -261237 // 250000
    A[5, 4] = 19663 // 15625
    A[6, 1] = 1336883 // 4909125
    A[6, 2] = 0
    A[6, 3] = -25476 // 30875
    A[6, 4] = 194159 // 185250
    A[6, 5] = 8225 // 78546
    A[7, 1] = -2459386 // 14727375
    A[7, 2] = 0
    A[7, 3] = 19504 // 30875
    A[7, 4] = 2377474 // 13615875
    A[7, 5] = -6157250 // 5773131
    A[7, 6] = 902 // 735
    A[8, 1] = 2699 // 7410
    A[8, 2] = 0
    A[8, 3] = -252 // 1235
    A[8, 4] = -1393253 // 3993990
    A[8, 5] = 236875 // 72618
    A[8, 6] = -135 // 49
    A[8, 7] = 15 // 22
    A[9, 1] = 11 // 144
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 256 // 693
    A[9, 5] = 0
    A[9, 6] = 125 // 504
    A[9, 7] = 125 // 528
    A[9, 8] = 5 // 72
    α[1] = 11 // 144
    α[2] = 0
    α[3] = 0
    α[4] = 256 // 693
    α[5] = 0
    α[6] = 125 // 504
    α[7] = 125 // 528
    α[8] = 5 // 72
    αEEst[1] = 28 // 477
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 212 // 441
    αEEst[5] = -312500 // 366177
    αEEst[6] = 2125 // 1764
    αEEst[7] = 0
    αEEst[8] = -2105 // 35532
    αEEst[9] = 2995 // 17766

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst,
        stability_size = 4.3240258134292215))
end

"""

From Verner's Website
"""
function constructVernerEfficient6(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 3 // 50
    c[3] = 1439 // 15000
    c[4] = 1439 // 10000
    c[5] = 4973 // 10000
    c[6] = 389 // 400
    c[7] = 1999 // 2000
    c[8] = 1
    c[9] = 1
    A[2, 1] = 3 // 50
    A[3, 1] = 519479 // 27000000
    A[3, 2] = 2070721 // 27000000
    A[4, 1] = 1439 // 40000
    A[4, 2] = 0
    A[4, 3] = 4317 // 40000
    A[5, 1] = 109225017611 // 82828840000
    A[5, 2] = 0
    A[5, 3] = -417627820623 // 82828840000
    A[5, 4] = 43699198143 // 10353605000
    A[6, 1] = -8036815292643907349452552172369 // 191934985946683241245914401600
    A[6, 2] = 0
    A[6, 3] = 246134619571490020064824665 // 1543816496655405117602368
    A[6, 4] = -13880495956885686234074067279 // 113663489566254201783474344
    A[6, 5] = 755005057777788994734129 // 136485922925633667082436
    A[7, 1] = -BigInt(1663299841566102097180506666498880934230261) //
              BigInt(30558424506156170307020957791311384232000)
    A[7, 2] = 0
    A[7, 3] = 130838124195285491799043628811093033 // 631862949514135618861563657970240
    A[7, 4] = -BigInt(3287100453856023634160618787153901962873) //
              BigInt(20724314915376755629135711026851409200)
    A[7, 5] = 2771826790140332140865242520369241 // 396438716042723436917079980147600
    A[7, 6] = -1799166916139193 // 96743806114007800
    A[8, 1] = -BigInt(832144750039369683895428386437986853923637763) //
              BigInt(15222974550069600748763651844667619945204887)
    A[8, 2] = 0
    A[8, 3] = 818622075710363565982285196611368750 // 3936576237903728151856072395343129
    A[8, 4] = -BigInt(9818985165491658464841194581385463434793741875) //
              BigInt(61642597962658994069869370923196463581866011)
    A[8, 5] = BigInt(31796692141848558720425711042548134769375) //
              BigInt(4530254033500045975557858016006308628092)
    A[8, 6] = -14064542118843830075 // 766928748264306853644
    A[8, 7] = -1424670304836288125 // 2782839104764768088217
    A[9, 1] = 382735282417 // 11129397249634
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 5535620703125000 // 21434089949505429
    A[9, 5] = 13867056347656250 // 32943296570459319
    A[9, 6] = 626271188750 // 142160006043
    A[9, 7] = -51160788125000 // 289890548217
    A[9, 8] = 163193540017 // 946795234
    α[1] = 382735282417 // 11129397249634
    α[2] = 0
    α[3] = 0
    α[4] = 5535620703125000 // 21434089949505429
    α[5] = 13867056347656250 // 32943296570459319
    α[6] = 626271188750 // 142160006043
    α[7] = -51160788125000 // 289890548217
    α[8] = 163193540017 // 946795234
    αEEst[1] = 124310637869885675646798613 // 2890072468789466426596827670
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 265863151737164990361330921875 // 1113197271463372303940319369579
    αEEst[5] = 3075493557174030806536302953125 // 6843749922042323876546949699876
    αEEst[6] = 67798000008733879813263055 // 29532792147666737550036372
    αEEst[7] = -1099436585155390846238326375 // 15055706496446408859196167
    αEEst[8] = 26171252653086373181571802 // 368794478890732346033505
    αEEst[9] = 1 // 30

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst, fsal = true,
        stability_size = 4.8552743145651425))
end

function constructVerner7(T::Type = Float64)
    A = zeros(T, 10, 10)
    c = zeros(T, 10)
    α = zeros(T, 10)
    αEEst = zeros(T, 10)
    c[2] = 1 // 200
    c[3] = 49 // 450
    c[4] = 49 // 300
    c[5] = 911 // 2000
    c[6] = 3480084980 // 5709648941
    c[7] = 221 // 250
    c[8] = 37 // 40
    c[9] = 1
    c[10] = 1
    A[2, 1] = 1 // 200
    A[3, 1] = -4361 // 4050
    A[3, 2] = 2401 // 2025
    A[4, 1] = 49 // 1200
    A[4, 3] = 49 // 400
    A[5, 1] = 2454451729 // 3841600000
    A[5, 3] = -9433712007 // 3841600000
    A[5, 4] = 4364554539 // 1920800000
    A[6, 1] = -BigInt(6187101755456742839167388910402379177523537620) //
              BigInt(2324599620333464857202963610201679332423082271)
    A[6, 3] = BigInt(27569888999279458303270493567994248533230000) //
              BigInt(2551701010245296220859455115479340650299761)
    A[6, 4] = -BigInt(37368161901278864592027018689858091583238040000) //
              BigInt(4473131870960004275166624817435284159975481033)
    A[6, 5] = BigInt(1392547243220807196190880383038194667840000000) //
              BigInt(1697219131380493083996999253929006193143549863)
    A[7, 1] = 11272026205260557297236918526339 // 1857697188743815510261537500000
    A[7, 3] = -48265918242888069 // 1953194276993750
    A[7, 4] = 26726983360888651136155661781228 // 1308381343805114800955157615625
    A[7, 5] = -2090453318815827627666994432 // 1096684189897834170412307919
    A[7, 6] = BigInt(1148577938985388929671582486744843844943428041509) //
              BigInt(1141532118233823914568777901158338927629837500000)
    A[8, 1] = BigInt(1304457204588839386329181466225966641) //
              BigInt(108211771565488329642169667802016000)
    A[8, 3] = -1990261989751005 // 40001418792832
    A[8, 4] = BigInt(2392691599894847687194643439066780106875) //
              BigInt(58155654089143548047476915856270826016)
    A[8, 5] = -BigInt(1870932273351008733802814881998561250) //
              BigInt(419326053051486744762255151208232123)
    A[8, 6] = BigInt(1043329047173803328972823866240311074041739158858792987034783181) //
              BigInt(510851127745017966999893975119259285040213723744255237522144000)
    A[8, 7] = -311918858557595100410788125 // 3171569057622789618800376448
    A[9, 1] = BigInt(17579784273699839132265404100877911157) //
              BigInt(1734023495717116205617154737841023480)
    A[9, 3] = -18539365951217471064750 // 434776548575709731377
    A[9, 4] = BigInt(447448655912568142291911830292656995992000) //
              BigInt(12511202807447096607487664209063950964109)
    A[9, 5] = -BigInt(65907597316483030274308429593905808000000) //
              BigInt(15158061430635748897861852383197382130691)
    A[9, 6] = BigInt(273847823027445129865693702689010278588244606493753883568739168819449761) //
              BigInt(136252034448398939768371761610231099586032870552034688235302796640584360)
    A[9, 7] = BigInt(694664732797172504668206847646718750) //
              BigInt(1991875650119463976442052358853258111)
    A[9, 8] = -19705319055289176355560129234220800 // 72595753317320295604316217197876507
    A[10, 1] = -511858190895337044664743508805671 // 11367030248263048398341724647960
    A[10, 3] = 2822037469238841750 // 15064746656776439
    A[10, 4] = -BigInt(23523744880286194122061074624512868000) //
               BigInt(152723005449262599342117017051789699)
    A[10, 5] = BigInt(10685036369693854448650967542704000000) //
               BigInt(575558095977344459903303055137999707)
    A[10, 6] = -BigInt(6259648732772142303029374363607629515525848829303541906422993) //
               BigInt(876479353814142962817551241844706205620792843316435566420120)
    A[10, 7] = 17380896627486168667542032602031250 // 13279937889697320236613879977356033
    α[1] = 96762636172307789 // 2051985304794103980
    α[4] = 312188947591288252500000 // 1212357694274963646019729
    α[5] = 13550580884964304000000000000 // 51686919683339547115937980629
    α[6] = BigInt(72367769693133178898676076432831566019684378142853445230956642801) //
           BigInt(475600216991873963561768100160364792981629064220601844848928537580)
    α[7] = 1619421054120605468750 // 3278200730370057108183
    α[8] = -66898316144057728000 // 227310933007074849597
    α[9] = 181081444637946577 // 2226845467039736466
    αEEst[1] = 117807213929927 // 2640907728177740
    αEEst[4] = 4758744518816629500000 // 17812069906509312711137
    αEEst[5] = 1730775233574080000000000 // 7863520414322158392809673
    αEEst[6] = BigInt(2682653613028767167314032381891560552585218935572349997) //
               BigInt(12258338284789875762081637252125169126464880985167722660)
    αEEst[7] = 40977117022675781250 // 178949401077111131341
    αEEst[10] = 2152106665253777 // 106040260335225546
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, adaptiveorder = 6, αEEst = αEEst,
        stability_size = 4.640024335337557))
end

"""
Verner Efficient 8
"""
function constructVerner8(T::Type = Float64)
    A = zeros(T, 13, 13)
    c = zeros(T, 13)
    α = zeros(T, 13)
    αEEst = zeros(T, 13)

    c[2] = 1 // 20
    c[3] = 341 // 3200
    c[4] = 1023 // 6400
    c[5] = 39 // 100
    c[6] = 93 // 200
    c[7] = 31 // 200
    c[8] = 943 // 1000
    c[9] = 7067558016280 // 7837150160667
    c[10] = 909 // 1000
    c[11] = 47 // 50
    c[12] = 1
    c[13] = 1
    A[2, 1] = 1 // 20
    A[3, 1] = -7161 // 1024000
    A[3, 2] = 116281 // 1024000
    A[4, 1] = 1023 // 25600
    A[4, 3] = 3069 // 25600
    A[5, 1] = 4202367 // 11628100
    A[5, 3] = -3899844 // 2907025
    A[5, 4] = 3982992 // 2907025
    A[6, 1] = 5611 // 114400
    A[6, 4] = 31744 // 135025
    A[6, 5] = 923521 // 5106400
    A[7, 1] = 21173 // 343200
    A[7, 4] = 8602624 // 76559175
    A[7, 5] = -26782109 // 689364000
    A[7, 6] = 5611 // 283500
    A[8, 1] = -1221101821869329 // 690812928000000
    A[8, 4] = -125 // 2
    A[8, 5] = -1024030607959889 // 168929280000000
    A[8, 6] = 1501408353528689 // 265697280000000
    A[8, 7] = 6070139212132283 // 92502016000000
    A[9, 1] = -BigInt(1472514264486215803881384708877264246346044433307094207829051978044531801133057155) //
              BigInt(1246894801620032001157059621643986024803301558393487900440453636168046069686436608)
    A[9, 4] = -BigInt(5172294311085668458375175655246981230039025336933699114138315270772319372469280000) //
              BigInt(124619381004809145897278630571215298365257079410236252921850936749076487132995191)
    A[9, 5] = -BigInt(12070679258469254807978936441733187949484571516120469966534514296406891652614970375) //
              BigInt(2722031154761657221710478184531100699497284085048389015085076961673446140398628096)
    A[9, 6] = BigInt(780125155843893641323090552530431036567795592568497182701460674803126770111481625) //
              BigInt(183110425412731972197889874507158786859226102980861859505241443073629143100805376)
    A[9, 7] = BigInt(664113122959911642134782135839106469928140328160577035357155340392950009492511875) //
              BigInt(15178465598586248136333023107295349175279765150089078301139943253016877823170816)
    A[9, 8] = BigInt(10332848184452015604056836767286656859124007796970668046446015775000000) //
              BigInt(1312703550036033648073834248740727914537972028638950165249582733679393783)
    A[10, 1] = -BigInt(29055573360337415088538618442231036441314060511) //
               BigInt(22674759891089577691327962602370597632000000000)
    A[10, 4] = -20462749524591049105403365239069 // 454251913499893469596231268750
    A[10, 5] = -180269259803172281163724663224981097 // 38100922558256871086579832832000000
    A[10, 6] = BigInt(21127670214172802870128286992003940810655221489) //
               BigInt(4679473877997892906145822697976708633673728000)
    A[10, 7] = BigInt(318607235173649312405151265849660869927653414425413) //
               BigInt(6714716715558965303132938072935465423910912000000)
    A[10, 8] = 212083202434519082281842245535894 // 20022426044775672563822865371173879
    A[10, 9] = -BigInt(2698404929400842518721166485087129798562269848229517793703413951226714583) //
               BigInt(469545674913934315077000442080871141884676035902717550325616728175875000000)
    A[11, 1] = -BigInt(2342659845814086836951207140065609179073838476242943917) //
               BigInt(1358480961351056777022231400139158760857532162795520000)
    A[11, 4] = -996286030132538159613930889652 // 16353068885996164905464325675
    A[11, 5] = -26053085959256534152588089363841 // 4377552804565683061011299942400
    A[11, 6] = BigInt(20980822345096760292224086794978105312644533925634933539) //
               BigInt(3775889992007550803878727839115494641972212962174156800)
    A[11, 7] = BigInt(890722993756379186418929622095833835264322635782294899) //
               BigInt(13921242001395112657501941955594013822830119803764736)
    A[11, 8] = BigInt(161021426143124178389075121929246710833125) //
               BigInt(10997207722131034650667041364346422894371443)
    A[11, 9] = BigInt(300760669768102517834232497565452434946672266195876496371874262392684852243925359864884962513) //
               BigInt(4655443337501346455585065336604505603760824779615521285751892810315680492364106674524398280000)
    A[11, 10] = -31155237437111730665923206875 // 392862141594230515010338956291
    A[12, 1] = -BigInt(2866556991825663971778295329101033887534912787724034363) //
               BigInt(868226711619262703011213925016143612030669233795338240)
    A[12, 4] = -BigInt(16957088714171468676387054358954754000) //
               BigInt(143690415119654683326368228101570221)
    A[12, 5] = -BigInt(4583493974484572912949314673356033540575) //
               BigInt(451957703655250747157313034270335135744)
    A[12, 6] = BigInt(2346305388553404258656258473446184419154740172519949575) //
               BigInt(256726716407895402892744978301151486254183185289662464)
    A[12, 7] = BigInt(1657121559319846802171283690913610698586256573484808662625) //
               BigInt(13431480411255146477259155104956093505361644432088109056)
    A[12, 8] = BigInt(345685379554677052215495825476969226377187500) //
               BigInt(74771167436930077221667203179551347546362089)
    A[12, 9] = -BigInt(3205890962717072542791434312152727534008102774023210240571361570757249056167015230160352087048674542196011) //
               BigInt(947569549683965814783015124451273604984657747127257615372449205973192657306017239103491074738324033259120)
    A[12, 10] = BigInt(40279545832706233433100438588458933210937500) //
                BigInt(8896460842799482846916972126377338947215101)
    A[12, 11] = -BigInt(6122933601070769591613093993993358877250) //
                BigInt(1050517001510235513198246721302027675953)
    A[13, 1] = -BigInt(618675905535482500672800859344538410358660153899637) //
               BigInt(203544282118214047100119475340667684874292102389760)
    A[13, 4] = -BigInt(4411194916804718600478400319122931000) //
               BigInt(40373053902469967450761491269633019)
    A[13, 5] = -BigInt(16734711409449292534539422531728520225) //
               BigInt(1801243715290088669307203927210237952)
    A[13, 6] = BigInt(135137519757054679098042184152749677761254751865630525) //
               BigInt(16029587794486289597771326361911895112703716593983488)
    A[13, 7] = BigInt(38937568367409876012548551903492196137929710431584875) //
               BigInt(340956454090191606099548798001469306974758443147264)
    A[13, 8] = -BigInt(6748865855011993037732355335815350667265625) //
               BigInt(7002880395717424621213565406715087764770357)
    A[13, 9] = -BigInt(1756005520307450928195422767042525091954178296002788308926563193523662404739779789732685671) //
               BigInt(348767814578469983605688098046186480904607278021030540735333862087061574934154942830062320)
    A[13, 10] = BigInt(53381024589235611084013897674181629296875) //
                BigInt(8959357584795694524874969598508592944141)
    α[1] = 44901867737754616851973 // 1014046409980231013380680
    α[6] = 791638675191615279648100000 // 2235604725089973126411512319
    α[7] = 3847749490868980348119500000 // 15517045062138271618141237517
    α[8] = -13734512432397741476562500000 // 875132892924995907746928783
    α[9] = BigInt(12274765470313196878428812037740635050319234276006986398294443554969616342274215316330684448207141) //
           BigInt(489345147493715517650385834143510934888829280686609654482896526796523353052166757299452852166040)
    α[10] = -9798363684577739445312500000 // 308722986341456031822630699
    α[11] = 282035543183190840068750 // 12295407629873040425991
    α[12] = -306814272936976936753 // 1299331183183744997286

    αEEst[1] = 10835401739407019406577 // 244521829356935137978320
    αEEst[6] = 13908189778321895491375000 // 39221135527894265375640567
    αEEst[7] = 73487947527027243487625000 // 296504045773342769773399443
    αEEst[8] = 68293140641257649609375000 // 15353208647806945749946119
    αEEst[9] = BigInt(22060647948996678611017711379974578860522018208949721559448560203338437626022142776381) //
               BigInt(1111542009262325874512959185795727215759010577565736079641376621381577236680929558640)
    αEEst[10] = -547971229495642458203125000 // 23237214025700991642563601
    αEEst[13] = -28735456870978964189 // 79783493704265043693

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, adaptiveorder = 7, αEEst = αEEst,
        stability_size = 5.864113239353444))
end

"""

Papakostas's Order 6

On Phase-Fitted modified Runge-Kutta Pairs of order 6(5), by Ch. Tsitouras and I. Th. Famelis,
 International Conference of Numerical Analysis and Applied Mathematics, Crete, (2006)
"""
function constructPapakostas6(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 17 // 183
    c[3] = 12 // 83
    c[4] = 18 // 83
    c[5] = 71 // 125
    c[6] = 42 // 59
    c[7] = 199 // 200
    c[8] = 1
    c[9] = 1
    A[2, 1] = 17 // 183
    A[3, 1] = 3756 // 117113
    A[3, 2] = 13176 // 117113
    A[4, 1] = 9 // 166
    A[4, 2] = 0
    A[4, 3] = 27 // 166
    A[5, 1] = 207751751 // 316406250
    A[5, 2] = 0
    A[5, 3] = -526769377 // 210937500
    A[5, 4] = 1524242129 // 632812500
    A[6, 1] = -4970082682619223281 // 2887511529739311186
    A[6, 2] = 0
    A[6, 3] = 97919278033879057 // 13556392158400522
    A[6, 4] = -407131674007930877068 // 74078904949579652469
    A[6, 5] = 1237601855204268750000 // 1753200750473385108433
    A[7, 1] = 176597685527535385020980411 // 42773485015591331328000000
    A[7, 2] = 0
    A[7, 3] = -6793162515552646891859 // 401628967282547712000
    A[7, 4] = 12704926019361287204873446554247 // 886659402653054716778496000000
    A[7, 5] = -50728836334509259632278125 // 32657591718008685915971584
    A[7, 6] = 51536223982796190703 // 51293749413888000000
    A[8, 1] = 299033520572337573523 // 66918720793812357519
    A[8, 2] = 0
    A[8, 3] = -16550269823961899 // 902146153892364
    A[8, 4] = 49920346343238033627496282 // 3215735869387500624775563
    A[8, 5] = -1686432488955761721093750 // 978844996793357447730403
    A[8, 6] = 161901609084039 // 149698803705724
    A[8, 7] = -305146137600000 // 54760341991955873
    A[9, 1] = 24503 // 381483
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 1366847103121 // 4106349847584
    A[9, 5] = 20339599609375 // 75933913767768
    A[9, 6] = 35031290651 // 194765546144
    A[9, 7] = 16620160000000 // 11001207123543
    A[9, 8] = -14933 // 11016
    α[1] = 24503 // 381483
    α[2] = 0
    α[3] = 0
    α[4] = 1366847103121 // 4106349847584
    α[5] = 20339599609375 // 75933913767768
    α[6] = 35031290651 // 194765546144
    α[7] = 16620160000000 // 11001207123543
    α[8] = -14933 // 11016
    αEEst[1] = 61010485298317 // 979331468960880
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 320207313882553286621 // 941222813406992395200
    αEEst[5] = 6845867841119140625 // 29008216787127405534
    αEEst[6] = 124109197949158875473 // 562495660250110816320
    αEEst[7] = 19339714537078400000 // 16810691577722216811
    αEEst[8] = -211029377951 // 210416202900
    αEEst[9] = -1 // 150

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst, fsal = true,
        stability_size = 4.449984142037165))
end

"""

Lawson's Order 6

An Order 6 Runge-Kutta Process with an Extended Region of Stability, by J. D. Lawson,
 Siam Journal on Numerical Analysis, Vol. 4, No. 4 (Dec. 1967) pages 620-625.
"""
function constructLawson6(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)

    c[2] = 26 // 105 - 2 // 315 * 51^(1 // 2)
    c[3] = 13 // 35 - 1 // 105 * 51^(1 // 2)
    c[4] = 7 // 8
    c[5] = 1 // 2
    c[6] = 1 // 8
    c[7] = 1
    A[2, 1] = 26 // 105 - 2 // 315 * 51^(1 // 2)
    A[3, 1] = 13 // 140 - 1 // 420 * 51^(1 // 2)
    A[3, 2] = 39 // 140 - 1 // 140 * 51^(1 // 2)
    A[4, 1] = 917 // 1024 + 133 // 2048 * 51^(1 // 2)
    A[4, 2] = -3339 // 1024 - 567 // 2048 * 51^(1 // 2)
    A[4, 3] = 1659 // 512 + 217 // 1024 * 51^(1 // 2)
    A[5, 1] = -653 // 5684 - 1313 // 45472 * 51^(1 // 2)
    A[5, 2] = 477 // 464 + 81 // 928 * 51^(1 // 2)
    A[5, 3] = -2478015 // 5433904 - 21365 // 339619 * 51^(1 // 2)
    A[5, 4] = 14568 // 339619 + 1528 // 339619 * 51^(1 // 2)
    A[6, 1] = 9515 // 50176 + 1443 // 100352 * 51^(1 // 2)
    A[6, 2] = -477 // 1024 - 81 // 2048 * 51^(1 // 2)
    A[6, 3] = 25166643 // 29980160 - 2073611 // 59960320 * 51^(1 // 2)
    A[6, 4] = 111 // 46844 - 275 // 46844 * 51^(1 // 2)
    A[6, 5] = -141 // 320 + 21 // 320 * 51^(1 // 2)
    A[7, 1] = -790 // 49 - 2243 // 3528 * 51^(1 // 2)
    A[7, 2] = 159 // 4 + 27 // 8 * 51^(1 // 2)
    A[7, 3] = -44723003 // 702660 - 258368 // 175665 * 51^(1 // 2)
    A[7, 4] = 327392 // 316197 - 2528 // 105399 * 51^(1 // 2)
    A[7, 5] = 3158 // 135 - 56 // 45 * 51^(1 // 2)
    A[7, 6] = 448 // 27
    α[1] = 1 // 70
    α[2] = 0
    α[3] = 0
    α[4] = 256 // 945
    α[5] = 58 // 135
    α[6] = 256 // 945
    α[7] = 1 // 70

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 6.463163498924733))
end

"""

Tsitouras-Papakostas's Order 6

Cheap Error Estimation for Runge-Kutta methods, by Ch. Tsitouras and S.N. Papakostas,
Siam Journal on Scientific Computing, Vol. 20, Issue 6, Nov 1999.
"""
function constructTsitourasPapakostas6(T::Type = Float64)
    A = zeros(T, 8, 8)
    c = zeros(T, 8)
    α = zeros(T, 8)
    αEEst = zeros(T, 8)

    c[2] = 4 // 27
    c[3] = 2 // 9
    c[4] = 3 // 7
    c[5] = 11 // 16
    c[6] = 10 // 13
    c[7] = 1
    c[8] = 1
    A[2, 1] = 4 // 27
    A[3, 1] = 1 // 18
    A[3, 2] = 1 // 6
    A[4, 1] = 66 // 343
    A[4, 2] = -729 // 1372
    A[4, 3] = 1053 // 1372
    A[5, 1] = 13339 // 49152
    A[5, 2] = -4617 // 16384
    A[5, 3] = 5427 // 53248
    A[5, 4] = 95207 // 159744
    A[6, 1] = -6935 // 57122
    A[6, 2] = 23085 // 48334
    A[6, 3] = 33363360 // 273642941
    A[6, 4] = 972160 // 118442467
    A[6, 5] = 172687360 // 610434253
    A[7, 1] = 611 // 1891
    A[7, 2] = -4617 // 7564
    A[7, 3] = 6041007 // 13176488
    A[7, 4] = 12708836 // 22100117
    A[7, 5] = -35840000 // 62461621
    A[7, 6] = 6597591 // 7972456
    A[8, 1] = -11107621 // 14976720
    A[8, 2] = 193023 // 332816
    A[8, 3] = 1515288591 // 724706840
    A[8, 4] = -801736138 // 352888965
    A[8, 5] = 814013312 // 606245145
    A[8, 6] = 0
    A[8, 7] = 0
    α[1] = 131 // 1800
    α[2] = 0
    α[3] = 1121931 // 3902080
    α[4] = 319333 // 1682928
    α[5] = 262144 // 2477325
    α[6] = 4084223 // 15177600
    α[7] = 1891 // 25200
    α[8] = 0
    αEEst[1] = 4093 // 60720
    αEEst[2] = 0
    αEEst[3] = 16273467 // 51284480
    αEEst[4] = 693889 // 5376020
    αEEst[5] = 146456576 // 626763225
    αEEst[6] = 15737111 // 93089280
    αEEst[7] = 1 // 25
    αEEst[8] = 1 // 23

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst,
        stability_size = 3.954129730631186))
end

"""

DormandLockyerMcCorriganPrince Order 6 Global Error Estimation

Global Error estimation with Runge-Kutta triples, by J.R.Dormand, M.A.Lockyer, N.E.McCorrigan and P.J.Prince,
 Computers and Mathematics with Applications, 18 (1989) pages 835-846.
"""
function constructDormandLockyerMcCorriganPrince6(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 4 // 39
    c[3] = 2 // 13
    c[4] = 3 // 13
    c[5] = 13021 // 22659
    c[6] = 39 // 67
    c[7] = 86 // 87
    c[8] = 1
    c[9] = 1
    A[2, 1] = 4 // 39
    A[3, 1] = 1 // 26
    A[3, 2] = 3 // 26
    A[4, 1] = 3 // 52
    A[4, 2] = 0
    A[4, 3] = 9 // 52
    A[5, 1] = 1406640413621 // 2478209482476
    A[5, 2] = 0
    A[5, 3] = -1755653396555 // 826069827492
    A[5, 4] = 1321105868272 // 619552370619
    A[6, 1] = 1463083752990765495750771 // 2783511518115684243554356
    A[6, 2] = 0
    A[6, 3] = -417501847634533557363 // 213770948323146013636
    A[6, 4] = 414679390177938970209021 // 208212903666744217281464
    A[6, 5] = 48490458547529962724706855 // 2711140218644676453221942744
    A[7, 1] = -1146771707244809451668952985850 // 1178428995610817474751161698881
    A[7, 2] = 0
    A[7, 3] = 883524649813655720289029 // 257840992692019416968811
    A[7, 4] = -550972740958654450507278066587 // 325473716439106878117398000544
    A[7, 5] = -968282586950392419883943203143069455 // 32828460835559176127341032228568032
    A[7, 6] = 3230428272165469542719 // 108684219530393291772
    A[8, 1] = -368234904360842614256649493 // 316816365518493517722015828
    A[8, 2] = 0
    A[8, 3] = 19247613365107236707 // 4836252322132133628
    A[8, 4] = -3627331815384766429266963559 // 1852940251413900055822878936
    A[8, 5] = -BigInt(1113629335962635330712822690622675431) //
              BigInt(31397585101055611361934971740947112)
    A[8, 6] = 7515696221383336 // 210977455127283
    A[8, 7] = -2466239729887929744 // 169565652664150899277
    A[9, 1] = 265211783 // 3930519060
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 53198489747 // 147124055808
    A[9, 5] = -165257255035734106911939 // 60068315060067285425920
    A[9, 6] = 1687862952891371 // 536423949472320
    A[9, 7] = 397538864947251 // 474823057340620
    A[9, 8] = -7142267 // 10794560
    α[1] = 265211783 // 3930519060
    α[2] = 0
    α[3] = 0
    α[4] = 53198489747 // 147124055808
    α[5] = -165257255035734106911939 // 60068315060067285425920
    α[6] = 1687862952891371 // 536423949472320
    α[7] = 397538864947251 // 474823057340620
    α[8] = -7142267 // 10794560
    αEEst[1] = 61909109615135759874805 // 909927606074377568473224
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 191026171670373376608531013 // 532182573462887006324068800
    αEEst[5] = -16134045022121298883366692111903 // 6240618480477179290700503465280
    αEEst[6] = 123307988643675607258315033 // 41372548832224521099482580
    αEEst[7] = 1648598887654728061640361417 // 1901782716664438059984219160
    αEEst[8] = -17 // 25
    αEEst[9] = -259237562821839 // 28937895739220050

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, adaptiveorder = 5, αEEst = αEEst, fsal = true,
        stability_size = 4.371949555657634))
end

"""

constructTanakaKasugaYamashitaYazaki Order 6 D

On the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,
 by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,
 Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.
"""
function constructTanakaKasugaYamashitaYazaki6D(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 250
    c[3] = 61 // 500
    c[4] = 333 // 1000
    c[5] = 56 // 125
    c[6] = 221 // 250
    c[7] = 31 // 40
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 250
    A[3, 1] = -3477 // 2000
    A[3, 2] = 3721 // 2000
    A[4, 1] = 5893 // 960
    A[4, 2] = -155333 // 24000
    A[4, 3] = 2 // 3
    A[5, 1] = 24108461262911872104361127416 // 10623418096530694059554537475
    A[5, 2] = -5917919395479054157557671728 // 2612315925376400178578984625
    A[5, 3] = 232280630406488701857197728 // 1274810171583683287146544497
    A[5, 4] = 5471535080448445690971520 // 20898527403011201428631877
    A[6, 1] = BigInt(6437510603218083862126582571557692575111484748570978316432703829) //
              BigInt(462655382163643505993146289902892129798145865902401535270912000)
    A[6, 2] = -25083125032238929237211514053 // 1741543950250933452385989750
    A[6, 3] = BigInt(191471592827838703414373042295007003190385413366181230722958611) //
              BigInt(201998644533947923598793335502244867385083328952030670310246400)
    A[6, 4] = -BigInt(723497411149192638041504597211363403675892243966786723532767) //
              BigInt(584075332429951001052530262062711448991365353250616692288000)
    A[6, 5] = BigInt(20368395873479894700136894919415569811519) //
              BigInt(12245317457885350045378058014027636736000)
    A[7, 1] = -BigInt(1053235524590544505552847720336831796873415344431325525436387723231012409077771) //
              BigInt(151882847900737440255262782278159042659114250070600779935950570906179610869760)
    A[7, 2] = 110556263584594415870053017499108315 // 15711834237341169934278611772340752
    A[7, 3] = BigInt(27883202411953328597139497795565290131372048742640205053329023156822546337651) //
              BigInt(66313136270946971682878125476803367732416846682610519097035561761715919388672)
    A[7, 4] = -BigInt(134178686375307065766767719823759096451245561636180411142603830822435042239) //
              BigInt(191743203036279335802347994313396039656741862830639889775698722784806666240)
    A[7, 5] = BigInt(1192365386732594028658711096219658222748871123673792181) //
              BigInt(1339984907884349782764828286663427050229770700990709760)
    A[7, 6] = 1 // 16
    A[8, 1] = BigInt(139472860369418680318405579593202554799958275392817771935632674919979204348213061279) //
              BigInt(37743426679794053215512832350173350810071958850079214921681182582585691688393293824)
    A[8, 2] = -764277739442303538051509575419150250 // 241488900344982122797870353007980909
    A[8, 3] = -BigInt(7371553444464925806572282862981502319484039203525663191183225508630714910425342675) //
              BigInt(4846779318912629207559289658412491740053883371241789941360845137312080681624033792)
    A[8, 4] = BigInt(8801132642676112554927251590711440081546416618547191842516081677278335900180085) //
              BigInt(2071689842640920660980285733506294870283516476694921779219912685665930091210112)
    A[8, 5] = -BigInt(298368164571373898500032163301238477118159212135073380804555) //
              BigInt(92828694411466377636989261345454987067336292167034303922176)
    A[8, 6] = -240786613447142598518650 // 1482036690339563123748203
    A[8, 7] = 52860332724128242560 // 47245264125077723987
    A[9, 1] = -BigInt(70335440697678472884389500304621917516220835017456192585613672496070708987720514364185979355621) //
              BigInt(5219173052726261361794285880851433852375067076864459834761850509124293962243718819226202095616)
    A[9, 2] = BigInt(55057602528492150957512473843952756334951250) //
              BigInt(3825542672812430870769605520946554801303339)
    A[9, 3] = -BigInt(77980782561050861485751338404689118536953083418300092316918550174770045337324006382574405615) //
              BigInt(39424363920896283457695961995433007195554371610645094834381872512385376722954438065398121984)
    A[9, 4] = BigInt(2177533338045948502564188609357006384381697842038103984795932632439956259386338338574157919) //
              BigInt(387582376874948376545286159615164080595940845312336117877987331598679530792396290947668352)
    A[9, 5] = -BigInt(1429465993644112145770349737782620061868013755109155817982734692969685809) //
              BigInt(295236772314145005435682097407317274403116518334589230312060489793912832)
    A[9, 6] = -BigInt(162882286796273095967517045376244881875035) //
              BigInt(1226056358355165076488310595889719561316528)
    A[9, 7] = 61 // 43
    A[9, 8] = 0
    α[1] = 4783097999 // 163657290888
    α[2] = 0
    α[3] = 1043945712500000 // 4811958427553991
    α[4] = 12731480000000 // 212031396478881
    α[5] = 3420427578125 // 10842916254408
    α[6] = 22394286718750 // 335690566308279
    α[7] = 601452300800 // 2318247773037
    α[8] = 11555307017 // 221477311314
    α[9] = 0
    αEEst[1] = 401359714631498171030059910 // 15267054426382006481740339287
    αEEst[2] = 0
    αEEst[3] = 1872374247198383241346435120000 // 8180262298260999309783933205623
    αEEst[4] = 215533585653762668869792000 // 12486209626895034583833718599
    αEEst[5] = 10643952265127448732528470125 // 29318848171314267144678228486
    αEEst[6] = 49992752929730967217476426250 // 454671173972183983205823101943
    αEEst[7] = 41077106106719584492498868480 // 192232681858829974979332840131
    αEEst[8] = 0
    αEEst[9] = 13 // 318

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 7.723403386984321))
end

"""

constructTanakaKasugaYamashitaYazaki Order 6 C

On the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,
 by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,
 Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.
"""
function constructTanakaKasugaYamashitaYazaki6C(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 100
    c[3] = 11 // 100
    c[4] = 33 // 100
    c[5] = 43 // 100
    c[6] = 177 // 200
    c[7] = 39 // 50
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 100
    A[3, 1] = -99 // 200
    A[3, 2] = 121 // 200
    A[4, 1] = 931 // 600
    A[4, 2] = -1133 // 600
    A[4, 3] = 2 // 3
    A[5, 1] = 91958765679448742111 // 109106230010372040600
    A[5, 2] = -2894503641421476941 // 3306249394253698200
    A[5, 3] = 42536885251557644 // 181843716683953401
    A[5, 4] = 124767011204926940 // 545531150051860203
    A[6, 1] = BigInt(2432582077139652291040741961658278931048241123) //
              BigInt(1042517831971586211298945274530293397770600000)
    A[6, 2] = -5086533956627126527 // 2404545014002689600
    A[6, 3] = -BigInt(28233375209799471901757749652001556425049843) //
              BigInt(193956805948202085823059585959124353073600000)
    A[6, 4] = -BigInt(64143717120482360099075816554051641763250417) //
              BigInt(80815335811750869092941494149635147114000000)
    A[6, 5] = 368355687225024998241998861697 // 229321716687851233749224000000
    A[7, 1] = -BigInt(218146655544801488428299362829574333746140908535614198091) //
              BigInt(129671440168534479253493745415914920591854154832474480000)
    A[7, 2] = 3225815231796122073236882 // 1899026582543293621912725
    A[7, 3] = BigInt(90474621818961969641326264897298400503244755498598833923) //
              BigInt(192999352808981550516827900153919881811131765332055040000)
    A[7, 4] = -BigInt(1469845402576002022391874744301978646600063888157276277) //
              BigInt(1827645386448688925348749054487877668666020505038400000)
    A[7, 5] = BigInt(236537380548621166783390004339783981292483) //
              BigInt(228189687334919622052120277661453593600000)
    A[7, 6] = 1 // 16
    A[8, 1] = BigInt(20745045233020258374773011834908499926526821480199726607317) //
              BigInt(6703021612773885250214340203318686755448684469761576087200)
    A[8, 2] = -10062163770239468036471200 // 3356073465031554356885103
    A[8, 3] = -BigInt(6600722659737487671039676673294259975246135762925896417) //
              BigInt(9474444736357442876691573864537554989039417965607475200)
    A[8, 4] = BigInt(80225477702329113082476652949959330669319741271791060069) //
              BigInt(21801989308095252074205042131464260060005478841312656000)
    A[8, 5] = -BigInt(572101553192975384084889663939632620168841743) //
              BigInt(192662470105149826757328803363450656933632000)
    A[8, 6] = -223351084874296955 // 945158103640587564
    A[8, 7] = 3931692514354 // 3491636633667
    A[9, 1] = BigInt(1042961260483134998083167734826788831236189026066652724982180836527) //
              BigInt(166222076505494822625791597828801737913859767695104851064439232400)
    A[9, 2] = -4754359005740226898713768141325 // 838855890397374389475867956178
    A[9, 3] = -BigInt(42750676060439966787211838095564987064848978722439603628501619959) //
              BigInt(27488922212795009426487499382540855882395413779035427566212172800)
    A[9, 4] = BigInt(205686747816311239126907254218859168711623965851131881271189228377) //
              BigInt(77312593723485964011996092013396157169237101253537140029971736000)
    A[9, 5] = -BigInt(75946282944430327533335124997772458772050053465158929) //
              BigInt(62109575359095533697158937720595685916206114892672000)
    A[9, 6] = -264375495891732597330199815785 // 1220640854393687869579464225198
    A[9, 7] = 8 // 11
    A[9, 8] = 0
    α[1] = 92725517 // 4525454934
    α[2] = 0
    α[3] = 203217125 // 939422946
    α[4] = 16794200 // 510230259
    α[5] = 7552075 // 21858018
    α[6] = 22967296000 // 562167528741
    α[7] = 229749950 // 798566769
    α[8] = 67671311 // 1203893922
    α[9] = 0
    αEEst[1] = 1773713828809813620047 // 76546041056100425133324
    αEEst[2] = 0
    αEEst[3] = 576936961199857272475 // 2856614455209710336064
    αEEst[4] = 28907772101756401885 // 257621374646831423844
    αEEst[5] = 26154985271235615605 // 103780652247492444224
    αEEst[6] = -12285231659070809284400 // 206713297399385556110331
    αEEst[7] = 237440851439576891930 // 613972674142959478347
    αEEst[8] = 0
    αEEst[9] = 1 // 12

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 9.727143096983585))
end

"""

constructTanakaKasugaYamashitaYazaki Order 6 B

On the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,
 by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,
 Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.
"""
function constructTanakaKasugaYamashitaYazaki6B(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 5
    c[3] = 3 // 20
    c[4] = 2 // 5
    c[5] = 1 // 2
    c[6] = 3 // 4
    c[7] = 4 // 5
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 5
    A[3, 1] = 3 // 32
    A[3, 2] = 9 // 160
    A[4, 1] = -71 // 400
    A[4, 2] = -53 // 400
    A[4, 3] = 71 // 100
    A[5, 1] = 20861387 // 246441600
    A[5, 2] = -2674671 // 27382400
    A[5, 3] = 14952973 // 61610400
    A[5, 4] = 555163 // 2053680
    A[6, 1] = 197611764884981873 // 575894510914896000
    A[6, 2] = 16376643 // 136912000
    A[6, 3] = -281745778473052331 // 1007815394101068000
    A[6, 4] = -761016212927233 // 4799120924290800
    A[6, 5] = 7119424000 // 9814726677
    A[7, 1] = 2141650356947829667 // 15684127070697870750
    A[7, 2] = 124460793 // 710231000
    A[7, 3] = 729393184893351907 // 54894444747442547625
    A[7, 4] = -857270676309077671 // 4182433885519432200
    A[7, 5] = 496013770812769 // 855353429900550
    A[7, 6] = 1 // 10
    A[8, 1] = 59331555925548335981 // 3757538915366421931200
    A[8, 2] = -807212549 // 1636098400
    A[8, 3] = 3149950445245624075843 // 6575693101891238379600
    A[8, 4] = 213748931907436689007 // 250502594357761462080
    A[8, 5] = -229203323900597 // 10246103543907504
    A[8, 6] = -1306535 // 1043952
    A[8, 7] = 35275 // 24856
    A[9, 1] = -82982860685611836065717 // 78861151462250930572800
    A[9, 2] = -4258977307 // 34337529600
    A[9, 3] = 257897426390213500944749 // 138007015058939128502400
    A[9, 4] = 7684590445667080762721 // 5257410097483395371520
    A[9, 5] = -418846778208596923 // 215039562243682176
    A[9, 6] = 14813735 // 606173568
    A[9, 7] = 10 // 13
    A[9, 8] = 0
    α[1] = 253 // 6048
    α[2] = 0
    α[3] = 71200 // 292383
    α[4] = 725 // 7056
    α[5] = 124 // 441
    α[6] = -160 // 1323
    α[7] = 10375 // 26208
    α[8] = 239 // 4284
    α[9] = 0
    αEEst[1] = 1061069 // 26544672
    αEEst[2] = 0
    αEEst[3] = 17566690 // 75486411
    αEEst[4] = 3286405 // 10322928
    αEEst[5] = -148157 // 1935549
    αEEst[6] = 3649510 // 5806647
    αEEst[7] = -27404525 // 115026912
    αEEst[8] = 0
    αEEst[9] = 2 // 21

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 6.007860841958881))
end

"""

TanakaKasugaYamashitaYazaki Order 6 A

On the Optimization of Some Eight-stage Sixth-order Explicit Runge-Kutta Method,
 by M. Tanaka, K. Kasuga, S. Yamashita and H. Yazaki,
 Journal of the Information Processing Society of Japan, Vol. 34, No. 1 (1993), pages 62 to 74.
"""
function constructTanakaKasugaYamashitaYazaki6A(T::Type = Float64)
    A = zeros(T, 9, 9)
    c = zeros(T, 9)
    α = zeros(T, 9)
    αEEst = zeros(T, 9)

    c[2] = 1 // 100
    c[3] = 63 // 500
    c[4] = 63 // 200
    c[5] = 91 // 200
    c[6] = 21 // 25
    c[7] = 121 // 200
    c[8] = 1
    c[9] = 1
    A[2, 1] = 1 // 100
    A[3, 1] = -3339 // 5000
    A[3, 2] = 3969 // 5000
    A[4, 1] = 7409 // 2400
    A[4, 2] = -2751 // 800
    A[4, 3] = 2 // 3
    A[5, 1] = 879267900536971128212809 // 854114648586066093506400
    A[5, 2] = -32049516354411396813791 // 31633875873558003463200
    A[5, 3] = 415928927703664756394 // 3202929932197747850649
    A[5, 4] = 4945815476584113561248 // 16014649660988739253245
    A[6, 1] = BigInt(264364422150553063350374034500534568610052754809099) //
              BigInt(76594645373513572971811613560122656201571540854400)
    A[6, 2] = -4942738251715985006237 // 988558621048687608225
    A[6, 3] = BigInt(857070668647981630049588712158864150151776268553) //
              BigInt(235675831918503301451728041723454326774066279552)
    A[6, 4] = -BigInt(575165070977427883238172569267537197276341339157) //
              BigInt(157117221279002200967818694482302884516044186368)
    A[6, 5] = 41539505110086641686706761875 // 17218030314852652247321445632
    A[7, 1] = -BigInt(7384332843467384251748880311988170834105043085902030057004642579) //
              BigInt(1315398593371357506652159044711447383311024734178303038737049600)
    A[7, 2] = 5656750329476775131209405217353 // 831860122007390417396028410400
    A[7, 3] = -BigInt(333346623364947547771242996577246247659950898456372585977245839) //
              BigInt(190226873502934777885081461850578544663440500019631516371204096)
    A[7, 4] = BigInt(7100883704616862706920092769411499000150383998766128086388523) //
              BigInt(4497089208107205150947552289611785925849657210865993294827520)
    A[7, 5] = -BigInt(2179578925723796720158217363641743790096295) //
              BigInt(4632538471751007430155293544377007681933312)
    A[7, 6] = 1 // 16
    A[8, 1] = -BigInt(2032945134330261992294105795025746391376544519999264654886089003925834569) //
              BigInt(370204608795106233716533306294462150635128163745469969879178202980349056)
    A[8, 2] = 104400636802444703223540033778325 // 10372562449415222120508569935887
    A[8, 3] = -BigInt(1282736937871158931724790946085686543914309745805108288599439251128389375) //
              BigInt(121675640652937013878860597173704343215741424447811808281967661119415424)
    A[8, 4] = BigInt(69591790067148620077849942355176681820069270105493261193720164528698455) //
              BigInt(5177686836295192079951514773349120987903890402034545033275219622102784)
    A[8, 5] = -BigInt(277247257204675144642737406657510117972977064771225085) //
              BigInt(32594437558796549389422402191407749206480247247281408)
    A[8, 6] = 206057847729988755370340 // 362497215363232474360143
    A[8, 7] = 482015452915471040 // 328927231085842477
    A[9, 1] = BigInt(149325572050522001715919965533803802996087500406626797668220209959218571) //
              BigInt(32753013690352528622466162111781793695504165546622151035212132158026496)
    A[9, 2] = -BigInt(2802049259729837394660870725638222600) //
              BigInt(886389794903952023701108666481909127)
    A[9, 3] = -BigInt(247046029472637555138263024224503808349122551711211979490561742618517516481875) //
              BigInt(62386925133959095591686249139496495395591773904073928005189977758511970759424)
    A[9, 4] = BigInt(17088511387117619026306232328319322341432000611501090216270103426420066685195) //
              BigInt(2654762771657833854965372303808361506195394634215911830008084159936679606784)
    A[9, 5] = -BigInt(1905127968294488657830239165534292922200553202371049568415) //
              BigInt(428517757803501533934438346448102436732262910333700649472)
    A[9, 6] = BigInt(348519474109306281943020669241827485) //
              BigInt(1086125473870224085200713153977775988)
    A[9, 7] = 5 // 4
    A[9, 8] = 0
    α[1] = 703450322 // 19272872619
    α[2] = 0
    α[3] = 1197916718750000 // 6147559103622993
    α[4] = 20671732000 // 146035199457
    α[5] = 5917108000 // 36937869969
    α[6] = 23017645375 // 95720436504
    α[7] = 23374664000 // 131057876103
    α[8] = 4163735081 // 86609369112
    α[9] = 0
    αEEst[1] = 572720177479183548778421 // 9495052764959466387724734
    αEEst[2] = 0
    αEEst[3] = 4729369086248342935597812500 // 45049042428818825645161957701
    αEEst[4] = 1262015363899500143599090 // 3413510704758155271620049
    αEEst[5] = -106945493251106516170 // 843200958241401898187
    αEEst[6] = 711860694955815806929355 // 3483263158495231241410497
    αEEst[7] = 1763218383368811180347180 // 5312516873137942812595713
    αEEst[8] = 0
    αEEst[9] = 1 // 18

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 4.730235594027788))
end

"""

Mikkawy-Eisa Order 6

A general four-parameter non-FSAL embedded Runge–Kutta algorithm of orders 6 and 4 in seven stages,
 by M.E.A. El-Mikkawy and M.M.M. Eisa,
 Applied Mathematics and Computation, Vol. 143, No. 2, (2003) pages 259 to 267.

"""
function constructMikkawyEisa(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)
    αEEst = zeros(T, 7)

    c[2] = 3 // 19
    c[3] = 9 // 38
    c[4] = 342 // 683
    c[5] = 22 // 27
    c[6] = 9 // 11
    c[7] = 1
    A[2, 1] = 3 // 19
    A[3, 1] = 9 // 152
    A[3, 2] = 27 // 152
    A[4, 1] = 94474764 // 318611987
    A[4, 2] = -310753854 // 318611987
    A[4, 3] = 375818328 // 318611987
    A[5, 1] = -1073687692 // 12631821129
    A[5, 2] = 1421618 // 911979
    A[5, 3] = -797223265256 // 505937677851
    A[5, 4] = 8812260835312 // 9612815879169
    A[6, 1] = -140767614 // 1662675883
    A[6, 2] = 302157 // 185009
    A[6, 3] = -6054413848056 // 3590051357651
    A[6, 4] = 1284501562654185 // 1332472080292672
    A[6, 5] = -34222143195 // 4324664670656
    A[7, 1] = 666881867 // 2248862211
    A[7, 2] = -19133 // 162361
    A[7, 3] = 337642670416 // 1442156395089699
    A[7, 4] = 388258038731717539 // 1050707717099750592
    A[7, 5] = 11900400291315 // 1115308365632
    A[7, 6] = -1730758084946 // 169374832839
    α[1] = 456023 // 6094440
    α[2] = 0
    α[3] = 85078761640 // 257589787311
    α[4] = 5499272526534791 // 27584466984914400
    α[5] = 4653828837 // 944530400
    α[6] = -962601827 // 208639800
    α[7] = 162361 // 1977800
    αEEst[1] = 851969693 // 10815599520
    αEEst[2] = 0
    αEEst[3] = 6675679897963 // 21163704130305
    αEEst[4] = 14982723272271221 // 67174249160747520
    αEEst[5] = 86506623776547 // 20114719398400
    αEEst[6] = -4
    αEEst[7] = 162361 // 1977800

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 3.932606644665211))
end

"""

Chummund's First Order 6 method

A three-dimensional family of seven-step Runge-Kutta methods of order 6, by G. M. Chammud (Hammud),
Numerical Methods and programming, 2001, Vol.2, 2001, pages 159-166
(Advanced Computing Scientific journal published by the Research Computing Center of the Lomonosov Moscow State Univeristy)
"""
function constructChummund6(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)

    c[2] = 4 // 7
    c[3] = 5 // 7
    c[4] = 6 // 7
    c[5] = 1 // 2 - 1 // 10 * 5^(1 // 2)
    c[6] = 1 // 2 + 1 // 10 * 5^(1 // 2)
    c[7] = 1
    A[2, 1] = 4 // 7
    A[3, 1] = 115 // 112
    A[3, 2] = -5 // 16
    A[4, 1] = 589 // 630
    A[4, 2] = 5 // 18
    A[4, 3] = -16 // 45
    A[5, 1] = -29 // 6000 * 5^(1 // 2) + 229 // 1200
    A[5, 2] = -187 // 1200 * 5^(1 // 2) + 119 // 240
    A[5, 3] = 34 // 375 * 5^(1 // 2) - 14 // 75
    A[5, 4] = -3 // 100 * 5^(1 // 2)
    A[6, 1] = 71 // 2400 - 587 // 12000 * 5^(1 // 2)
    A[6, 2] = 187 // 480 - 391 // 2400 * 5^(1 // 2)
    A[6, 3] = -38 // 75 + 26 // 375 * 5^(1 // 2)
    A[6, 4] = 27 // 80 - 3 // 400 * 5^(1 // 2)
    A[6, 5] = 1 // 4 * 5^(1 // 2) + 1 // 4
    A[7, 1] = 43 // 160 * 5^(1 // 2) - 49 // 480
    A[7, 2] = 51 // 32 * 5^(1 // 2) - 425 // 96
    A[7, 3] = -4 // 5 * 5^(1 // 2) + 52 // 15
    A[7, 4] = 3 // 16 * 5^(1 // 2) - 27 // 16
    A[7, 5] = -3 // 4 * 5^(1 // 2) + 5 // 4
    A[7, 6] = 5 // 2 - 1 // 2 * 5^(1 // 2)
    α[1] = 1 // 12
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 5 // 12
    α[6] = 5 // 12
    α[7] = 1 // 12

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 2.9052272747717307))
end

"""

Chummund's Second Order 6 method

A three-dimensional family of seven-step Runge-Kutta methods of order 6, by G. M. Chammud (Hammud),
Numerical Methods and programming, 2001, Vol.2, 2001, pages 159-166
(Advanced Computing Scientific journal published by the Research Computing Center of the Lomonosov Moscow State Univeristy)
"""
function constructChummund62(T::Type = Float64)
    A = zeros(T, 7, 7)
    c = zeros(T, 7)
    α = zeros(T, 7)

    c[2] = 750557 // 18870600
    c[3] = 748997 // 1685240
    c[4] = 624713 // 833636
    c[5] = 1 // 2 - 1 // 10 * 5^(1 // 2)
    c[6] = 1 // 2 + 1 // 10 * 5^(1 // 2)
    c[7] = 1
    A[2, 1] = 750557 // 18870600
    A[3, 1] = -BigInt(2618524936181374161531835574563010037) //
              BigInt(1243120420984397996179713098690223560) +
              3789435780183611636743 // 41406596067467160927808051177 *
              826321815619^(1 // 2)
    A[3, 2] = 79275599164011507825658766043059867 // 31078010524609949904492827467255589 -
              3789435780183611636743 // 41406596067467160927808051177 *
              826321815619^(1 // 2)
    A[4, 1] = -BigInt(268026226098849050460940570252841432600176823648465011650077) //
              BigInt(5150236239747870879984989921782016384957547807596571966124) +
              BigInt(75921270833932035355544618340621749422299940322072982) //
              BigInt(1287559059936967719996247480445504096239386951899142991531) *
              826321815619^(1 // 2)
    A[4, 2] = BigInt(31518353641767830515869136922708314461169361966520350023632789875) //
              BigInt(553066166067338591675176701039997652878686435650554094276321374) -
              BigInt(35636316203172669679112815645100268396443783685272961209975) //
              BigInt(553066166067338591675176701039997652878686435650554094276321374) *
              826321815619^(1 // 2)
    A[4, 3] = -BigInt(2316654837809440895270062146979974449795559590002488266576771951) //
              BigInt(551916642155010749548578353874329493980040828668519599378437454) +
              BigInt(3018332469640261586051209426912089268215645947971953376787) //
              BigInt(551916642155010749548578353874329493980040828668519599378437454) *
              826321815619^(1 // 2)
    A[5, 1] = 7557928766176693537 // 21071504819547814620 -
              3808101371666611447 // 21071504819547814620 * 5^(1 // 2) +
              1830930386299 // 21071504819547814620 * 826321815619^(1 // 2) +
              2746507673963 // 21071504819547814620 * 826321815619^(1 // 2) * 5^(1 // 2)
    A[5, 2] = -306676166043910279990585069725 // 1347799564492996272083268634084 +
              281193211701441857849531780175 // 1347799564492996272083268634084 *
              5^(1 // 2) -
              135831974232490143067425 // 1347799564492996272083268634084 *
              826321815619^(1 // 2) -
              203756277349889111562225 // 1347799564492996272083268634084 *
              826321815619^(1 // 2) * 5^(1 // 2)
    A[5, 3] = 59599519641651017073742851703 // 154850028368562600922782726012 -
              19070452237413934132405792717 // 154850028368562600922782726012 * 5^(1 // 2) +
              3249940179998629473139 // 154850028368562600922782726012 *
              826321815619^(1 // 2) +
              4875109240133213731043 // 154850028368562600922782726012 *
              826321815619^(1 // 2) * 5^(1 // 2)
    A[5, 4] = -179555217756675442542605916547 // 11203233050287498632877256911020 -
              10652462513464648264279700899 // 2240646610057499726575451382204 *
              5^(1 // 2) -
              79525200332352343539019 // 11203233050287498632877256911020 *
              826321815619^(1 // 2) -
              119292669246564200876603 // 11203233050287498632877256911020 *
              826321815619^(1 // 2) * 5^(1 // 2)
    A[6, 1] = -6195448976102809771 // 10535752409773907310 -
              1420753767230413123 // 10535752409773907310 * 5^(1 // 2) +
              2689284093484 // 5267876204886953655 * 826321815619^(1 // 2) +
              57223580479 // 5267876204886953655 * 826321815619^(1 // 2) * 5^(1 // 2)
    A[6, 2] = 583482891804126952797345387075 // 673899782246498136041634317042 +
              55992621594804402807018397650 // 336949891123249068020817158521 * 5^(1 // 2) -
              254249689295131189339575 // 336949891123249068020817158521 *
              826321815619^(1 // 2) -
              49153620024967407924975 // 673899782246498136041634317042 *
              826321815619^(1 // 2) * 5^(1 // 2)
    A[6, 3] = 19364321623328820251751318839 // 77425014184281300461391363006 +
              1040656060288496357918934418 // 38712507092140650230695681503 * 5^(1 // 2) -
              11049041979964718221319 // 38712507092140650230695681503 *
              826321815619^(1 // 2) -
              19382660452189380569035 // 77425014184281300461391363006 *
              826321815619^(1 // 2) * 5^(1 // 2)
    A[6, 4] = BigInt(115027957037506741565069794992944458832461) //
              BigInt(237309179139566721503969545043331027555510) -
              BigInt(19371908299625437579435144814464589438224) //
              BigInt(118654589569783360751984772521665513777755) * 5^(1 // 2) +
              BigInt(64263508455127995016232148201271346) //
              BigInt(118654589569783360751984772521665513777755) * 826321815619^(1 // 2) -
              BigInt(41163600684913040546114995026601907) //
              BigInt(237309179139566721503969545043331027555510) * 826321815619^(1 // 2) *
              5^(1 // 2)
    A[6, 5] = 43435632261 // 211822050005 * 5^(1 // 2) - 43432990301 // 84728820002 +
              205839 // 423644100010 * 826321815619^(1 // 2) * 5^(1 // 2) -
              514 // 42364410001 * 826321815619^(1 // 2)
    A[7, 1] = 9047270149938488929 // 4214300963909562924 +
              2216536302042479231 // 1404766987969854308 * 5^(1 // 2) -
              991800665293 // 1404766987969854308 * 826321815619^(1 // 2) * 5^(1 // 2) -
              12588066760235 // 4214300963909562924 * 826321815619^(1 // 2)
    A[7, 2] = -4301448087821718128020528522125 // 1347799564492996272083268634084 -
              2525818490403297345388026853875 // 1347799564492996272083268634084 *
              5^(1 // 2) +
              1510317586999119637060875 // 1347799564492996272083268634084 *
              826321815619^(1 // 2) * 5^(1 // 2) +
              5764153657065074502128625 // 1347799564492996272083268634084 *
              826321815619^(1 // 2)
    A[7, 3] = -491640814441543287886227446905 // 154850028368562600922782726012 +
              752920605871714580844952275 // 1564141700692551524472552788 * 5^(1 // 2) +
              1711626851729573101365 // 1564141700692551524472552788 *
              826321815619^(1 // 2) * 5^(1 // 2) +
              204731138699301217060685 // 154850028368562600922782726012 *
              826321815619^(1 // 2)
    A[7, 4] = -BigInt(222449163212150849227576515026323419478375) //
              BigInt(94923671655826688601587818017332411022204) +
              BigInt(805495551976820677052288591080735143509) //
              BigInt(958824966220471602036240586033660717396) * 5^(1 // 2) -
              BigInt(253684995628220543918448388547756365) //
              BigInt(94923671655826688601587818017332411022204) * 826321815619^(1 // 2) +
              BigInt(882636009291931420795306294477883) //
              BigInt(958824966220471602036240586033660717396) * 826321815619^(1 // 2) *
              5^(1 // 2)
    A[7, 5] = -44506854521 // 84728820002 * 5^(1 // 2) + 214493500755 // 42364410001 -
              205839 // 84728820002 * 826321815619^(1 // 2) * 5^(1 // 2) +
              2570 // 42364410001 * 826321815619^(1 // 2)
    A[7, 6] = 5 // 2 - 1 // 2 * 5^(1 // 2)
    α[1] = 1 // 12
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 5 // 12
    α[6] = 5 // 12
    α[7] = 1 // 12

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 2.87113148977662))
end

"""

Anton Hutas Second Order 6 method

Une amélioration de la méthode de Runge-Kutta-Nyström pour la résolution numérique des équations différentielles du premièr ordre,
by Anton Huta,
Acta Fac. Nat. Univ. Comenian Math., Vol. 1, pages 201-224 (1956).

"""
function constructHuta62(T::Type = Float64)
    A = zeros(T, 8, 8)
    c = zeros(T, 8)
    α = zeros(T, 8)

    c[2] = 1 // 9
    c[3] = 1 // 6
    c[4] = 1 // 3
    c[5] = 1 // 2
    c[6] = 2 // 3
    c[7] = 5 // 6
    c[8] = 1
    A[2, 1] = 1 // 9
    A[3, 1] = 1 // 24
    A[3, 2] = 1 // 8
    A[4, 1] = 1 // 6
    A[4, 2] = -1 // 2
    A[4, 3] = 2 // 3
    A[5, 1] = 139 // 272
    A[5, 2] = -945 // 544
    A[5, 3] = 105 // 68
    A[5, 4] = 99 // 544
    A[6, 1] = -53 // 3
    A[6, 2] = 91 // 2
    A[6, 3] = -52 // 3
    A[6, 4] = -107 // 6
    A[6, 5] = 8
    A[7, 1] = 55487 // 22824
    A[7, 2] = -83 // 16
    A[7, 3] = 2849 // 1902
    A[7, 4] = 34601 // 15216
    A[7, 5] = -640 // 2853
    A[7, 6] = 107 // 2536
    A[8, 1] = -101195 // 25994
    A[8, 2] = 351 // 41
    A[8, 3] = -35994 // 12997
    A[8, 4] = -26109 // 25994
    A[8, 5] = -10000 // 12997
    A[8, 6] = -36 // 12997
    A[8, 7] = 36 // 41
    α[1] = 41 // 840
    α[2] = 0
    α[3] = 9 // 35
    α[4] = 9 // 280
    α[5] = 34 // 105
    α[6] = 9 // 280
    α[7] = 9 // 35
    α[8] = 41 // 840

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 4.042886866789528))
end

"""

Anton Hutas First Order 6 method

Une amélioration de la méthode de Runge-Kutta-Nyström pour la résolution numérique des équations différentielles du premièr ordre,
by Anton Huta,
Acta Fac. Nat. Univ. Comenian Math., Vol. 1, pages 201-224 (1956).

"""
function constructHuta6(T::Type = Float64)
    A = zeros(T, 8, 8)
    c = zeros(T, 8)
    α = zeros(T, 8)

    c[2] = 1 // 9
    c[3] = 1 // 6
    c[4] = 1 // 3
    c[5] = 1 // 2
    c[6] = 2 // 3
    c[7] = 5 // 6
    c[8] = 1
    A[2, 1] = 1 // 9
    A[3, 1] = 1 // 24
    A[3, 2] = 1 // 8
    A[4, 1] = 1 // 6
    A[4, 2] = -1 // 2
    A[4, 3] = 2 // 3
    A[5, 1] = -5 // 8
    A[5, 2] = 27 // 8
    A[5, 3] = -3
    A[5, 4] = 3 // 4
    A[6, 1] = 221 // 9
    A[6, 2] = -981 // 9
    A[6, 3] = 289 // 3
    A[6, 4] = -34 // 3
    A[6, 5] = 1 // 9
    A[7, 1] = -61 // 16
    A[7, 2] = 113 // 8
    A[7, 3] = -59 // 6
    A[7, 4] = -11 // 8
    A[7, 5] = 5 // 3
    A[7, 6] = 1 // 16
    A[8, 1] = 358 // 41
    A[8, 2] = -2079 // 82
    A[8, 3] = 501 // 41
    A[8, 4] = 417 // 41
    A[8, 5] = -227 // 41
    A[8, 6] = -9 // 82
    A[8, 7] = 36 // 41
    α[1] = 41 // 840
    α[2] = 0
    α[3] = 9 // 35
    α[4] = 9 // 280
    α[5] = 34 // 105
    α[6] = 9 // 280
    α[7] = 9 // 35
    α[8] = 41 // 840

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 6, stability_size = 3.8400244379056505))
end

"""
The Relative Efficiency of Alternative Defect Control Schemes for High-Order Continuous Runge-Kutta Formulas
 W. H. Enright SIAM Journal on Numerical Analysis, Vol. 30, No. 5. (Oct., 1993), pp. 1419-1445.

"""
function constructEnrightVerner7(T::Type = Float64)
    A = zeros(T, 10, 10)
    c = zeros(T, 10)
    α = zeros(T, 10)
    αEEst = zeros(T, 10)
    c[2] = 1 // 18
    c[3] = 1 // 9
    c[4] = 1 // 6
    c[5] = 4 // 9
    c[6] = 19 // 39
    c[7] = 7 // 9
    c[8] = 8 // 9
    c[9] = 1
    c[10] = 1
    A[2, 1] = 1 // 18
    A[3, 1] = 0
    A[3, 2] = 1 // 9
    A[4, 1] = 1 // 24
    A[4, 2] = 0
    A[4, 3] = 1 // 8
    A[5, 1] = 44 // 81
    A[5, 2] = 0
    A[5, 3] = -56 // 27
    A[5, 4] = 160 // 81
    A[6, 1] = 91561 // 685464
    A[6, 2] = 0
    A[6, 3] = -12008 // 28561
    A[6, 4] = 55100 // 85683
    A[6, 5] = 29925 // 228488
    A[7, 1] = -1873585 // 1317384
    A[7, 2] = 0
    A[7, 3] = 15680 // 2889
    A[7, 4] = -4003076 // 1083375
    A[7, 5] = -43813 // 21400
    A[7, 6] = 5751746 // 2287125
    A[8, 1] = 50383360 // 12679821
    A[8, 2] = 0
    A[8, 3] = -39440 // 2889
    A[8, 4] = 1258442432 // 131088375
    A[8, 5] = 222872 // 29425
    A[8, 6] = -9203268152 // 1283077125
    A[8, 7] = 24440 // 43197
    A[9, 1] = -22942833 // 6327608
    A[9, 2] = 0
    A[9, 3] = 71784 // 5947
    A[9, 4] = -572980 // 77311
    A[9, 5] = -444645 // 47576
    A[9, 6] = 846789710 // 90281407
    A[9, 7] = -240750 // 707693
    A[9, 8] = 3972375 // 14534468
    A[10, 1] = 3379947 // 720328
    A[10, 2] = 0
    A[10, 3] = -10656 // 677
    A[10, 4] = 78284 // 7447
    A[10, 5] = 71865 // 5416
    A[10, 6] = -2803372 // 218671
    A[10, 7] = 963000 // 886193
    A[10, 8] = 0
    A[10, 9] = 0
    α[1] = 28781 // 595840
    α[2] = 0
    α[3] = 0
    α[4] = 820752 // 3128125
    α[5] = 11259 // 280000
    α[6] = 188245551 // 625100000
    α[7] = 8667 // 43120
    α[8] = 286011 // 2737280
    α[9] = 5947 // 140000
    α[10] = 0
    αEEst[1] = 577 // 10640
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 8088 // 34375
    αEEst[5] = 3807 // 10000
    αEEst[6] = -1113879 // 16150000
    αEEst[7] = 8667 // 26180
    αEEst[8] = 0
    αEEst[9] = 0
    αEEst[10] = 677 // 10000

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 4.499873870264921))
end

"""
From Verner's website
"""
function constructVernerRobust7(T::Type = Float64)
    A = zeros(T, 10, 10)
    c = zeros(T, 10)
    α = zeros(T, 10)
    αEEst = zeros(T, 10)

    c[2] = 1 // 200
    c[3] = 49 // 450
    c[4] = 49 // 300
    c[5] = 91 // 200
    c[6] = 34704460 // 57271701
    c[7] = 167 // 200
    c[8] = 183 // 200
    c[9] = 1
    c[10] = 1
    A[2, 1] = 1 // 200
    A[3, 1] = -4361 // 4050
    A[3, 2] = 2401 // 2025
    A[4, 1] = 49 // 1200
    A[4, 2] = 0
    A[4, 3] = 49 // 400
    A[5, 1] = 1781 // 2800
    A[5, 2] = 0
    A[5, 3] = -13689 // 5600
    A[5, 4] = 507 // 224
    A[6, 1] = -275776923568321554889485313326460 // 108782544039075797415588764982099
    A[6, 2] = 0
    A[6, 3] = 9576001512158705648097438730000 // 929765333667314507825544999847
    A[6, 4] = -22178538465642954902290458689600 // 2789296001001943523476634999541
    A[6, 5] = 12323686750109283453854913779200 // 15540363434153685345084109283157
    A[7, 1] = 8221508014471276464993 // 8206108584945090400000
    A[7, 2] = 0
    A[7, 3] = -9852144759099 // 2364568872400
    A[7, 4] = 15322550932778398907299 // 3996134347464283007200
    A[7, 5] = -18338463121898520004 // 36506562121215938675
    A[7, 6] = BigInt(23340475544602125119307511373519383499) //
              BigInt(34957329425893779598660175543672800000)
    A[8, 1] = 81088643022740545033730780169 // 2975182110231937140152800000
    A[8, 2] = 0
    A[8, 3] = -2837794586103 // 67559110640
    A[8, 4] = -14167575606881316095038341141 // 1344719188593468553337836000
    A[8, 5] = 2395552232834276839307772 // 29760062864388068268875
    A[8, 6] = -BigInt(4076715891031001341580357765362043260356514682697) //
              BigInt(60535801523558513633981092635987721507186400000)
    A[8, 7] = 36551527355459957808 // 2801171464968864455
    A[9, 1] = -3347747115771808477876752833599 // 1101327591307901549464211073280
    A[9, 2] = 0
    A[9, 3] = 1214704878477625125 // 119815105452943264
    A[9, 4] = -65581118617864038124456452415 // 10200342297342072539428709672
    A[9, 5] = -133373082911575479273298406095 // 84070826916189821955373830704
    A[9, 6] = BigInt(622515683654039386383701463758952267447736841050281950137693) //
              BigInt(328994218860140584540186142455568117669077094653332432085760)
    A[9, 7] = 46169188671551441399013763125 // 2343692704156841275628930358208
    A[9, 8] = 18880867865877597493091875 // 3469664148196911643475533504
    A[10, 1] = -74309815528722196127526037 // 51427190037752065777334640
    A[10, 2] = 0
    A[10, 3] = 913722369253875 // 113761793498234
    A[10, 4] = -440658227159292602060396890 // 58109996881684093545238767
    A[10, 5] = 37290888293935763341814380 // 10411746696360295961914239
    A[10, 6] = -BigInt(645152888113581065780360392105179310452496326847) //
               BigInt(264735425121804814898131042131367320451487910960)
    A[10, 7] = 473757486190086838681892500 // 556321100802942639360459887
    A[10, 8] = 0
    A[10, 9] = 0
    α[1] = 9420080774669597 // 198627609019792680
    α[2] = 0
    α[3] = 0
    α[4] = 18658605936510000 // 72821569629535727
    α[5] = 296950875175030000 // 1101802245630054969
    α[6] = BigInt(18875276980274212686824835566143151189850553896330009) //
           BigInt(148780947139609706104394157596648357994575577036224440)
    α[7] = 18663850606812500 // 74993696164706319
    α[8] = 179884749312500 // 58508928482581269
    α[9] = 349315176247648 // 7273791403140339
    α[10] = 0
    αEEst[1] = 7362904929137 // 155056681514280
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 7505178129270000 // 29317774785916981
    αEEst[5] = 1851744839320000 // 6843492208882329
    αEEst[6] = BigInt(750882778189818437512810407839645051128089) //
               BigInt(6004363295715536270735789992319176322209240)
    αEEst[7] = 27902602073000000 // 110704980052661709
    αEEst[8] = 0
    αEEst[9] = 0
    αEEst[10] = 331667036438 // 6791588611709

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 4.635489309972309))
end

"""
From Verner's website
"""
function constructVernerEfficient7(T::Type = Float64)
    A = zeros(T, 10, 10)
    c = zeros(T, 10)
    α = zeros(T, 10)
    αEEst = zeros(T, 10)

    c[2] = 1 // 200
    c[3] = 1633 // 15000
    c[4] = 1633 // 10000
    c[5] = 911 // 2000
    c[6] = 3872020203200 // 6348224000949
    c[7] = 221 // 250
    c[8] = 37 // 40
    c[9] = 1
    c[10] = 1
    A[2, 1] = 1 // 200
    A[3, 1] = -2421739 // 2250000
    A[3, 2] = 2666689 // 2250000
    A[4, 1] = 1633 // 40000
    A[4, 2] = 0
    A[4, 3] = 4899 // 40000
    A[5, 1] = 13638791441 // 21333512000
    A[5, 2] = 0
    A[5, 3] = -10484391993 // 4266702400
    A[5, 4] = 1212514581 // 533337800
    A[6, 1] = -BigInt(10566420537453573046911093467384714791794598586757942301092800) //
              BigInt(3945478828013189718441489268080494867217306323934034936408479)
    A[6, 2] = 0
    A[6, 3] = BigInt(5231021510299954960251585060495720489874110055798400000000) //
              BigInt(481214639347870437668189933904195007588401795820714103721)
    A[6, 4] = -BigInt(53191822847901986858463766895922715840689391279600297600000000) //
              BigInt(6327491292785148384899029440906260154779895213246569749827429)
    A[6, 5] = BigInt(1781287527622331824616139968260471142273403953102720000000) //
              BigInt(2161611109404685052753814119556349840946763773078759856171)
    A[7, 1] = BigInt(40316732614812499600926954604381622079953) //
              BigInt(6622027015251398289270845012723400000000)
    A[7, 2] = 0
    A[7, 3] = -3102627096349187411472 // 125153701152335532175
    A[7, 4] = BigInt(24512593611811814197027847310011747907776) //
              BigInt(1196407956242848578947668299729517712775)
    A[7, 5] = -27992220377104300725776520574439264 // 14697934672024287924310069954740819
    A[7, 6] = BigInt(3894588020195704963647351313475389164102598021642890694420881) //
              BigInt(3885780076153082703909872726256212057602657456212198200000000)
    A[8, 1] = BigInt(219473830564834473126851658008599711669913643457) //
              BigInt(18148581918776897851104105262715070449854720000)
    A[8, 2] = 0
    A[8, 3] = -10861924990722808125 // 217658610699713969
    A[8, 4] = BigInt(119480318847905366486490498699126516657779132592500) //
              BigInt(2895978969529640099701689414747730700037272591519)
    A[8, 5] = -BigInt(1668548580652290279936994936449972411950403625) //
              BigInt(374242591559339035424931213674065468041890542)
    A[8, 6] = BigInt(38257721102788015064997204763622283132936825913246654962493266309964200473188255317) //
              BigInt(18800609083266437975749341539616748050410201552864110184725432140737391745319680000)
    A[8, 7] = -500096105036391897537714558203625 // 5080166748283963062825254164179392
    A[9, 1] = BigInt(73088758866605608581823594159257042017139758291) //
              BigInt(7204156842026682625849538837116971615241846400)
    A[9, 2] = 0
    A[9, 3] = -69262713771106609055901000000 // 1623478443813548051039428229
    A[9, 4] = BigInt(14675490755672909998080009537477454813814604865188000000) //
              BigInt(410204279941884058406992680456716961361878302419014157)
    A[9, 5] = -BigInt(56001310541527335660935463926718909867728400000) //
              BigInt(12914842534226649788973963163793850259828428307)
    A[9, 6] = BigInt(1658790111713212343550432782788025386861251060331884987447839679378079088262228191320420318059) //
              BigInt(829167683860659102329644660243394316974911850817625934703082071608371922272124378570833366400)
    A[9, 7] = BigInt(397383910096482616626611556102242919581343750) //
              BigInt(1139945435719262247312274838475128039778393907)
    A[9, 8] = -BigInt(15924946453851683996573317555310225254108800) //
              BigInt(58718594835565436263767442911746035434825067)
    A[10, 1] = -BigInt(14946539621537434020656235571945101642930646111) //
               BigInt(332010640625300658857148487489318781645129600)
    A[10, 2] = 0
    A[10, 3] = 171118989155232139053000000 // 913940284754207909241943
    A[10, 4] = -BigInt(250906345110237533448815794107434670644623116000000) //
               BigInt(1630126656666301353726306757388291345384354510761)
    A[10, 5] = BigInt(100053022832376562190601031855550739901200000) //
               BigInt(5412177019773639317674437354179369686847329)
    A[10, 6] = -BigInt(741124542370420850881498253695230801778840260257881176783159378749906762930019) //
               BigInt(104572466684527096537276822312914502674692290031200239724554350786783199401600)
    A[10, 7] = BigInt(61868781182570033929134892454042771640656250) //
               BigInt(47410089451637179941093334205322110905900061)
    A[10, 8] = 0
    A[10, 9] = 0
    α[1] = 23316791871424559928103 // 494567385514963690893600
    α[2] = 0
    α[3] = 0
    α[4] = 167985897310649194506250000000000000 // 652486332022351662028939973252423307
    α[5] = 1146393512631783066800000000000 // 4366210223630112740458758316827
    α[6] = BigInt(896878436896140749138361006566548336071636396490005664374818159563751598674344336941) //
           BigInt(5901910841358241153066354693443102222070394003395707585243672949102475401519736720800)
    α[7] = 7313950190577733068066406250 // 14824175638294700286557272779
    α[8] = -16063489150654390383296000 // 54687057461871573111333849
    α[9] = 608799317735794481861 // 7494797304744605111718
    α[10] = 0
    αEEst[1] = 1774600368527162619 // 39781803854163746050
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 3268266429605273275000000000000 // 12237407528692429753538888074653
    αEEst[5] = 1613021830993125656000000000 // 7306908939590938710641463789
    αEEst[6] = BigInt(791180295543850558970031000997303838661625778390906982427125664631928) //
               BigInt(3622667153393802697274794410484953168258089439622231002370760754183025)
    αEEst[7] = 3939837293445123160156250 // 17217393308123926000647239
    αEEst[8] = 0
    αEEst[9] = 0
    αEEst[10] = 79854863595916300271 // 3925846207247174106138

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 4.640792701835124))
end

"""
Completely Imbedded Runge-Kutta Pairs, by P.W.Sharp and J.H.Verner, Siam Journal on Numerical Analysis, Vol.31, No.4.
(August 1994) pages 1169-1190.
"""
function constructSharpVerner7(T::Type = Float64)
    A = zeros(T, 12, 12)
    c = zeros(T, 12)
    α = zeros(T, 12)
    αEEst = zeros(T, 12)

    c[2] = 1 // 12
    c[3] = 4 // 27
    c[4] = 2 // 9
    c[5] = 5 // 9
    c[6] = 2 // 3
    c[7] = 1 // 6
    c[8] = 4 // 9
    c[9] = 3 // 4
    c[10] = 11 // 12
    c[11] = 1
    c[12] = 1
    A[2, 1] = 1 // 12
    A[3, 1] = 4 // 243
    A[3, 2] = 32 // 243
    A[4, 1] = 1 // 18
    A[4, 2] = 0
    A[4, 3] = 1 // 6
    A[5, 1] = 5 // 9
    A[5, 2] = 0
    A[5, 3] = -25 // 12
    A[5, 4] = 25 // 12
    A[6, 1] = 1 // 15
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = 1 // 3
    A[6, 5] = 4 // 15
    A[7, 1] = 319 // 3840
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = 161 // 1536
    A[7, 5] = -41 // 960
    A[7, 6] = 11 // 512
    A[8, 1] = 245 // 5184
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 1627 // 10368
    A[8, 5] = 151 // 1296
    A[8, 6] = -445 // 10368
    A[8, 7] = 1 // 6
    A[9, 1] = -556349853 // 7539261440
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = -4356175383 // 3015704576
    A[9, 5] = -814787343 // 1884815360
    A[9, 6] = 831004641 // 3015704576
    A[9, 7] = 355452237 // 235601920
    A[9, 8] = 107943759 // 117800960
    A[10, 1] = -68998698967 // 1063035863040
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = -767387292485 // 425214345216
    A[10, 5] = -205995991597 // 265758965760
    A[10, 6] = -22181208863 // 141738115072
    A[10, 7] = 26226796959 // 15502606336
    A[10, 8] = 1614200643 // 1107329024
    A[10, 9] = 187 // 329
    A[11, 1] = 24511479161 // 17979371520
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 3889847115 // 217931776
    A[11, 5] = 22028391 // 3681280
    A[11, 6] = 614528179 // 217931776
    A[11, 7] = -148401247 // 10215552
    A[11, 8] = -3122234829 // 318384704
    A[11, 9] = -4160 // 1221
    A[11, 10] = 15040 // 20757
    A[12, 1] = 5519 // 110880
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 0
    A[12, 5] = 0
    A[12, 6] = 83 // 560
    A[12, 7] = 932 // 3675
    A[12, 8] = 282123 // 1047200
    A[12, 9] = 2624 // 24255
    A[12, 10] = 3008 // 19635
    A[12, 11] = 37 // 2100
    α[1] = 5519 // 110880
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 83 // 560
    α[7] = 932 // 3675
    α[8] = 282123 // 1047200
    α[9] = 2624 // 24255
    α[10] = 3008 // 19635
    α[11] = 37 // 2100
    αEEst[1] = 15509 // 341880
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 6827 // 15540
    αEEst[7] = 22138 // 81585
    αEEst[8] = 78003 // 387464
    αEEst[9] = -64144 // 299145
    αEEst[10] = 623408 // 2179485
    αEEst[11] = 0
    αEEst[12] = -15 // 518

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6, fsal = true,
        stability_size = 4.622084785418665))
end

"""
Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum, by P.W.Sharp and E.Smart,
 Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.
"""
function constructSharpSmart7(T::Type = Float64)
    A = zeros(T, 11, 11)
    c = zeros(T, 11)
    α = zeros(T, 11)
    αEEst = zeros(T, 11)

    c[2] = 1 // 50
    c[3] = 27 // 125
    c[4] = 41 // 100
    c[5] = 57 // 100
    c[6] = 43 // 50
    c[7] = 2272510 // 11977321
    c[8] = 18 // 25
    c[9] = 5 // 6
    c[10] = 1
    c[11] = 1
    A[2, 1] = 1 // 50
    A[3, 1] = -594 // 625
    A[3, 2] = 729 // 625
    A[4, 1] = 451 // 21600
    A[4, 2] = 0
    A[4, 3] = 1681 // 4320
    A[5, 1] = 19 // 160
    A[5, 2] = 0
    A[5, 3] = 361 // 3104
    A[5, 4] = 3249 // 9700
    A[6, 1] = -31 // 200
    A[6, 2] = 0
    A[6, 3] = 520921 // 412056
    A[6, 4] = -17371 // 11640
    A[6, 5] = 132023 // 106200
    A[7, 1] = 25959766877768976976598957736980 // 487594514129628295945513157189933
    A[7, 2] = 0
    A[7, 3] = BigInt(347890318302644246405985993187156250) //
              BigInt(1321817402067092875750818220388519949)
    A[7, 4] = -1717046972617147709491116450178750 // 7467894926932728111586543618014237
    A[7, 5] = 29780304732725103577764751746216250 // 258912687002832625147067486467854423
    A[7, 6] = -302662548054389051180423185000 // 25662869164717278733974376694207
    A[8, 1] = 42409705291266846 // 416462256407406875
    A[8, 2] = 0
    A[8, 3] = 3247095172038 // 883201854817
    A[8, 4] = -518509279926 // 374238074075
    A[8, 5] = 435669225629732566638 // 393965828849029186615
    A[8, 6] = -6468694559114760 // 61945939006089637
    A[8, 7] = -8593750881095206170491007194502 // 3213504543545558150903880585625
    A[9, 1] = -1401024812030113404025 // 19887564677841032175639
    A[9, 2] = 0
    A[9, 3] = 13281373111234375 // 5150833217292744
    A[9, 4] = -50491693720625 // 29100752640072
    A[9, 5] = 8909776468783164583973193125 // 6271093223575470807674793192
    A[9, 6] = -4792324941735635008750 // 159776107397443897190271
    A[9, 7] = -BigInt(1532806290465891141166096531902118541769245) //
              BigInt(1203242011387872547807852011647420329982736)
    A[9, 8] = -7500029126894375 // 132689679447323376
    A[10, 1] = 36393032615434450612 // 324390586094889663425
    A[10, 2] = 0
    A[10, 3] = -1462401427649331250 // 154787214582248211
    A[10, 4] = 4135780451822750 // 874504037187843
    A[10, 5] = -2349378733647002895234008950 // 1090914599757106529355865311
    A[10, 6] = -78686605908422443750 // 52446632451499515953
    A[10, 7] = BigInt(2315079813491204524435067899365885119542372444358703) //
               BigInt(316169042039527157595235231573788308031260760584200)
    A[10, 8] = -33473047374792524975 // 32907430028856870472
    A[10, 9] = 5594658687556280397846 // 1893189870520997940175
    A[11, 1] = 2508607706701842363083 // 197875357745688550590720
    A[11, 2] = 0
    A[11, 3] = -5122833329940625 // 508724268374592
    A[11, 4] = 13293920580875 // 2874148408896
    A[11, 5] = -599188464780493707137440161875 // 277270064173229869784600732736
    A[11, 6] = -3601465055348923762849875 // 2146128454918752594358208
    A[11, 7] = BigInt(606030238246181777051198920509497430523044409408159) //
               BigInt(74752050141640998967813674460513197348653288024576)
    A[11, 8] = -1922750201834125 // 1941504226023936
    A[11, 9] = 12539348439579 // 3975412795840
    A[11, 10] = 0
    α[1] = 771570009067 // 14036203465200
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 28304779228000000 // 53707434325074117
    α[6] = -296881060859375 // 515060733835389
    α[7] = BigInt(744858303758379680905615939985761920312207508379) //
           BigInt(2487223884477764590764433396524922145673887618400)
    α[8] = -5118512171875 // 11763620626464
    α[9] = 136801854099 // 127885521925
    α[10] = 103626500437 // 1717635089268
    α[11] = 0
    αEEst[1] = 448234490819 // 8120946290580
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 7786773134600000 // 14452831163890377
    αEEst[6] = -408698637296875 // 567617951573694
    αEEst[7] = BigInt(4426705150369152638325381078278067803359) //
               BigInt(14828075230102658203818343670586143438076)
    αEEst[8] = -5004542378125 // 10330679593521
    αEEst[9] = 154806770859 // 124231649870
    αEEst[10] = 0
    αEEst[11] = 16 // 243

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 3.8994534847171285))
end

"""
On the Optimization of Some Nine-Stage Seventh-order Runge-Kutta Method, by M. Tanaka, S. Muramatsu and S. Yamashita,
Information Processing Society of Japan, Vol. 33, No. 12 (1992) pages 1512-1526.
"""
function constructTanakaYamashitaEfficient7(T::Type = Float64)
    A = zeros(T, 10, 10)
    c = zeros(T, 10)
    α = zeros(T, 10)
    αEEst = zeros(T, 10)

    c[2] = 36259 // 463869
    c[3] = 36259 // 309246
    c[4] = 36259 // 206164
    c[5] = 76401 // 153188
    c[6] = 5164663400901569152 // 6688924124988083687
    c[7] = 3486 // 3517
    c[8] = 44151 // 44173
    c[9] = 1
    c[10] = 1
    A[2, 1] = 36259 // 463869
    A[3, 1] = 36259 // 1236984
    A[3, 2] = 36259 // 412328
    A[4, 1] = 36259 // 824656
    A[4, 2] = 0
    A[4, 3] = 108777 // 824656
    A[5, 1] = 17751533712975975593187 // 24112920357813127230992
    A[5, 2] = 0
    A[5, 3] = -68331192803887602162951 // 24112920357813127230992
    A[5, 4] = 7825717455900471140481 // 3014115044726640903874
    A[6, 1] = -BigInt(2425518501234340256175806929031336393991001205323654593685210322030691047097621496102266496) //
              BigInt(201073929944556265242953373967503382318096046559546854970564286270157897072030532387737241)
    A[6, 2] = 0
    A[6, 3] = BigInt(126875939114499086848675646731069753055308007638565564293214808307459627250976287910912) //
              BigInt(2631823273838775215546306644775636213113650954300949659959480717139276934490785884841)
    A[6, 4] = -BigInt(18238165682427587123600563411903599919711680222699338744428834349094610403849667513626245575680) //
              BigInt(479212348415302218688412787744011607018072627155280851781784515530195350673833210517995172867)
    A[6, 5] = BigInt(74777425357290689294313120787550134201356775453356604582280658347816977407509825814840320) //
              BigInt(27848089034948481594251542168496834020714916243735255636715495143105044359669539013058107)
    A[7, 1] = BigInt(42210784012026021620512889337138957588173072058924928398799062235) //
              BigInt(401168555464694196502745570125544252560955194769351196028554688)
    A[7, 2] = 0
    A[7, 3] = -BigInt(53537582181289418572806048482253962781541488) //
              BigInt(128102133978061070595749084326726258918069)
    A[7, 4] = BigInt(6373437319382536771018620806214785516542915567996760353063349991182871200304) //
              BigInt(19178871740288180724887392022898914045213833131528843480576173243533301485)
    A[7, 5] = -BigInt(836513109281956728811652083904588515347012294160401579661057793958992) //
              BigInt(42189346226535262916910956145917457264775063492307360825161811325023)
    A[7, 6] = BigInt(10038768138260655813133796321688310283082351149893792474426644227234755871856831386997923013888351) //
              BigInt(8279123943002224665888560854425725483235895533066047643118716510648226939201056966728652698557760)
    A[8, 1] = BigInt(1454976871505621321312348899226731229297985195430097820532172928754404221419640982320963761) //
              BigInt(12687546780768188413911065021432924447284583965992535848754097389537051103097048673168256)
    A[8, 2] = 0
    A[8, 3] = -BigInt(1452249436938195913836212549773886207822959770792) //
              BigInt(3187825000852340545619892931005470986913487349)
    A[8, 4] = BigInt(3193785703967379485471835519262043520640585789136428552340853315619929163223926155626278646291801931779256) //
              BigInt(8816743814108800069900425523882492176796603795861854625575345408990649746129323017714575203134405597571)
    A[8, 5] = -BigInt(314398569508916946629277462588835135011587938712337655816458752800894863689255534896547161759213480) //
              BigInt(14507196201560052990013371105817112064769849230048646555812475120383456376679192045076337148816813)
    A[8, 6] = BigInt(5021633516852870452803558794670341128133410978274753232000155240629688617274518068065484524425884625107263111090060721584249881611265924113) //
              BigInt(3807402575192378287101053794016079417728266285278436439472658972755893033722804748992796724254152818232996309281540415603729279478920107136)
    A[8, 7] = -BigInt(894451839895008223904010765658125850176064186717638397881061173697811879745) //
              BigInt(186244934020117483847289332768639722211239803963523669807238114327710091115676)
    A[9, 1] = BigInt(152015786770038627019906826956584678402371493198250158080970494807155603994339) //
              BigInt(1319428594672311986480108760138089275639618425553698631283119461253421932416)
    A[9, 2] = 0
    A[9, 3] = -BigInt(19887569115365707672105043997835466942389220328) //
              BigInt(43451712251082409470704235239058276887205131)
    A[9, 4] = BigInt(6298831527954572673520838478029639446424615570453903300371170696118960335541193275024146681623960) //
              BigInt(17307483347318198085207889427954666589398911583434527253470846782562794571553580157056644256313)
    A[9, 5] = -BigInt(16267621644623777942279856217571823792451732234540266142050307930357537283432611648312520) //
              BigInt(747020211145282116967827947968990352912884924402891384654470989583659988117513448655559)
    A[9, 6] = BigInt(491920517345271821393960134665582163547632868347911487496995665146055538579545277983570189994492481977206720065882583432234119698425636137169515) //
              BigInt(371241970695441505578374965290296000309261530083026613438333515399198575394818137422626328203755084156959422247928840402063855870066548878130304)
    A[9, 7] = -BigInt(17535891839112183607157943692398769696531153719141528498448224128785868799210475) //
              BigInt(3881175428498724209649715816699297677268154716152409333146177577349474565697791732)
    A[9, 8] = -BigInt(31140449219386755112730831706895080247696102690585728771850210691242594436100540310) //
              BigInt(58531715707220748822628340615174217489020037018063169180406742622693159384762890406389)
    A[10, 1] = BigInt(24861126512935523838485032295435745281790804119672244200744512677831357181363) //
               BigInt(215828469302253893975010055544246846578750854407392771457340001283636121600)
    A[10, 2] = 0
    A[10, 3] = -76626859319946149305867456329803 // 167454524692981091214376557800
    A[10, 4] = BigInt(257532657386915224604779230484778835596042580268896440943054087972106955277512448850995064336363) //
               BigInt(707777528357579864776572552477247532276956780876653359042572831013312547307465249178438602200)
    A[10, 5] = -BigInt(103092665221253777021612043042409780416654274677686197534469014507504059634284484983141143) //
               BigInt(4735075386204034224907103653335170134874540866215348781137359896717512695961598377363000)
    A[10, 6] = BigInt(1318945254307068672853031172410281620677291556423152759282406612372948205789241763483098989903852936890735513699395545618802215742952753372919) //
               BigInt(995520191927224509158660659519643916330847017611189618002256023928790665495276022949114110343406997764203331763292012060684160018593393766400)
    A[10, 7] = -BigInt(2175691361381933486174620849991740173349017185199505364607841) //
               BigInt(482872625303278742130341621563226511344221688759361797916327450)
    A[10, 8] = -BigInt(11327601987184122343710458559595782081610122892585097) //
               BigInt(21251874884678431935286330856983429378055579208005268000)
    A[10, 9] = 0
    α[1] = BigInt(677260699094873524061210073954310211) //
           BigInt(13212228177645157882237395248920447488)
    α[2] = 0
    α[3] = 0
    α[4] = BigInt(5627843976805934592544586970647029617399366281651959837492864) //
           BigInt(20448796992082885248862284273169726631726393791864145954479875)
    α[5] = BigInt(1359735671458057021603668186882234273947181034928034734244224) //
           BigInt(4035225037829041960922838374222759264846456609494840689395475)
    α[6] = BigInt(3575764371063841994042920363615768888383369782579963896064642431626191680598750790399139608006651160426580137040859330533720256407) //
           BigInt(18833618269956378326078572170759846509476617594300797062242096554507068838086062412372695473217373611870290738365243380652826304000)
    α[7] = BigInt(14322850798205614664394883796805489119964080948503151) //
           BigInt(1692788382425178679633337406927131793062126418747780)
    α[8] = -BigInt(16735096417960349589058935251250023138290806176584545269411) //
           BigInt(128573843052304513208482301684749747737236254208431871400)
    α[9] = 33050288141543277444692395096256051 // 271248590133163812341791503489000
    α[10] = 0
    αEEst[1] = BigInt(962650826879437817605721930727384851) //
               BigInt(18874611682350225546053421784172067840)
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = BigInt(99703652969826806275610089806158069716600653757297413344) //
               BigInt(361062893830367886445877712954352019629670588714825566425)
    αEEst[5] = BigInt(17540887447270394964911517553576959050951784592644178144) //
               BigInt(52550888012671962193116521992300249584518949945887203425)
    αEEst[6] = BigInt(101855668513773837712956593596043266148479443244790887636953159551191054134940671472736229702711787350735239179) //
               BigInt(504322587935299170723833764883183242017770187561624249681119708768991642691172146267201689787026963930014131200)
    αEEst[7] = BigInt(179578338747395946570172802104016572846366090083599) //
               BigInt(31203472487100067827342625012481692038011546889360)
    αEEst[8] = -BigInt(500374162579884236288722085953024481890963958534161489781) //
               BigInt(5844265593286568782203740985670443078965284282201448700)
    αEEst[9] = 0
    αEEst[10] = 80

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 4.912538831797948))
end

"""
On the Optimization of Some Nine-Stage Seventh-order Runge-Kutta Method, by M. Tanaka, S. Muramatsu and S. Yamashita,
Information Processing Society of Japan, Vol. 33, No. 12 (1992) pages 1512-1526.
"""
function constructTanakaYamashitaStable7(T::Type = Float64)
    A = zeros(T, 10, 10)
    c = zeros(T, 10)
    α = zeros(T, 10)
    αEEst = zeros(T, 10)

    c[2] = 1288 // 14535
    c[3] = 644 // 4845
    c[4] = 322 // 1615
    c[5] = 65 // 258
    c[6] = 627862490 // 27724306937
    c[7] = 78 // 115
    c[8] = 95 // 124
    c[9] = 1
    c[10] = 1
    A[2, 1] = 1288 // 14535
    A[3, 1] = 161 // 4845
    A[3, 2] = 161 // 1615
    A[4, 1] = 161 // 3230
    A[4, 2] = 0
    A[4, 3] = 483 // 3230
    A[5, 1] = 196347867755 // 3561236836416
    A[5, 2] = 0
    A[5, 3] = 134004261625 // 1187078945472
    A[5, 4] = 149425089125 // 1780618418208
    A[6, 1] = BigInt(874723327324627172137139789673935509613630495) //
              BigInt(56881344496107103495850556251109088355454531158)
    A[6, 2] = 0
    A[6, 3] = BigInt(140089490273660861720564275306545765967660125) //
              BigInt(4375488038162084884296196634700699104265733166)
    A[6, 4] = -BigInt(2038049847879400647989164901369906650290192935250) //
              BigInt(47909406273855748440601205051655304842157645301117)
    A[6, 5] = BigInt(2992403630086592541124850354372857004595944160) //
              BigInt(168193972876080132251638040959092347542006316207)
    A[7, 1] = -2001378790961964301303250341598299 // 131178829335937360185206084581250
    A[7, 2] = 0
    A[7, 3] = -1366679891168526950613 // 3342867750190010177170
    A[7, 4] = -BigInt(197077954039191584877658472075693196650) //
              BigInt(14482289235786224954374999272581521053)
    A[7, 5] = BigInt(2928205733652489758138852423071126752) //
              BigInt(289494157432907631631314224968221875)
    A[7, 6] = BigInt(32572843800597493853254181634376441943013874856495312) //
              BigInt(1642002429836009758962688168840520197029337863346875)
    A[8, 1] = BigInt(226949925367094612475083609619198193642397605) //
              BigInt(120406368918742115270494114142742317844627456)
    A[8, 2] = 0
    A[8, 3] = 83451940525721530822125 // 1129677771714575730562048
    A[8, 4] = BigInt(7035716180093388934005544535766324331669337496890597125) //
              BigInt(8744668672303692797949525789789288711635037496457428992)
    A[8, 5] = BigInt(1783910495800307104322539337559667105512922125) //
              BigInt(384741540181237650862158213355907627413029681664)
    A[8, 6] = -BigInt(84153602056538973791098303633128803165153465256807063611390103929570560775) //
              BigInt(37282901435832588263568764858094381410386083732352113498074163165795975168)
    A[8, 7] = BigInt(529655154424978769932790603243342890625) //
              BigInt(2074272966571578715335103162383459680256)
    A[9, 1] = -BigInt(2220302447236283385210081868020072818509) //
              BigInt(374126802552343922668161638021420098000)
    A[9, 2] = 0
    A[9, 3] = -580875348986851918117575 // 7422906155739208262352728
    A[9, 4] = -BigInt(10152884092399228192381460837845336141124812794348025) //
              BigInt(1298474693887469810743537803700532058113206856602944)
    A[9, 5] = BigInt(57628597675871150072147324302138021982593246488) //
              BigInt(8694208064927983865022808524707850651868593125)
    A[9, 6] = BigInt(25006323928409346448859146781582297955220041834414805003311931721685248934905197) //
              BigInt(3116797087097920659891402066923116953398124957660879445267880838924530153376000)
    A[9, 7] = -BigInt(991935992163983524020354479671037652370649875) //
              BigInt(1096924827756227471690652450958041839154828032)
    A[9, 8] = BigInt(11028636941861502413824025771962099757599945728) //
              BigInt(10166706010345110864067934052134974581343819375)
    A[10, 1] = BigInt(2847557162233802909802913419338134005277) //
               BigInt(175580852316165047798631596921711256000)
    A[10, 2] = 0
    A[10, 3] = 1552914837310075 // 7167358071597822
    A[10, 4] = BigInt(7813795507395804332400817811705117266280297151179075) //
               BigInt(609385085239477995902119898391644620983028518989568)
    A[10, 5] = -BigInt(306358654025510315315806741256227901425369583) //
               BigInt(37780230596405542839291492772111971777338125)
    A[10, 6] = -BigInt(30230616135053261889365940573714713926595600173797519397657905897494488117634591) //
               BigInt(1462739064177653822206630548583223553011887775071057363338006334613651993472000)
    A[10, 7] = -3556025825918703192187464108779875 // 11170666795578957880984290260063232
    A[10, 8] = 2140578935503723938488131712 // 2556174768949326564363043125
    A[10, 9] = 0
    α[1] = -28836965799708194669 // 40897924114041540000
    α[2] = 0
    α[3] = 0
    α[4] = -5319231056637407390089058139231875 // 4078513870347642725257280732562048
    α[5] = 2430832495624902882205404599808 // 1640677233246140147577399278125
    α[6] = BigInt(58846832125102891510730576086257275195560457005949780449038492433452787) //
           BigInt(55051022626529988904867724618756855605285697694381512902783594390080000)
    α[7] = 8765694250492187515737289375 // 142862058843931893355781359104
    α[8] = 282726763309436004945812396032 // 864839130161188820942829590625
    α[9] = 35795813026789129771 // 507885604115513709330
    α[10] = 0
    αEEst[1] = -4303806316703372599 // 11685121175440440000
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = -2188638181830974432849378205625 // 2703688346269567600435718085888
    αEEst[5] = 1313681506776569792299836438 // 1214416900996402773928496875
    αEEst[6] = BigInt(3205635250634133320066291736997892470430172563677735127041069) //
               BigInt(5224293273794925935053617011014728265925425567271645535360000)
    αEEst[7] = 153260086062341088187716875 // 1103181921574763655256998912
    αEEst[8] = 1144102534493369691260897984 // 4260291281582210940605071875
    αEEst[9] = 0
    αEEst[10] = 3 // 40

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 9.29902700964581))
end

"""
Some Explicit Runge-Kutta Methods of High Order, by G. J. Cooper and J. H. Verner,
 SIAM Journal on Numerical Analysis, Vol. 9, No. 3, (September 1972), pages 389 to 405
"""
function constructCooperVerner8(T::Type = Float64)
    A = zeros(T, 11, 11)
    c = zeros(T, 11)
    α = zeros(T, 11)

    c[2] = 1 // 2
    c[3] = 1 // 2
    c[4] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[5] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[6] = 1 // 2
    c[7] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[8] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[9] = 1 // 2
    c[10] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[11] = 1
    A[2, 1] = 1 // 2
    A[3, 1] = 1 // 4
    A[3, 2] = 1 // 4
    A[4, 1] = 1 // 7
    A[4, 2] = -1 // 14 + 3 // 98 * 21^(1 // 2)
    A[4, 3] = 3 // 7 - 5 // 49 * 21^(1 // 2)
    A[5, 1] = 11 // 84 - 1 // 84 * 21^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = 2 // 7 - 4 // 63 * 21^(1 // 2)
    A[5, 4] = 1 // 12 + 1 // 252 * 21^(1 // 2)
    A[6, 1] = 5 // 48 - 1 // 48 * 21^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 1 // 4 - 1 // 36 * 21^(1 // 2)
    A[6, 4] = -77 // 120 - 7 // 180 * 21^(1 // 2)
    A[6, 5] = 63 // 80 + 7 // 80 * 21^(1 // 2)
    A[7, 1] = 5 // 21 + 1 // 42 * 21^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = -48 // 35 - 92 // 315 * 21^(1 // 2)
    A[7, 4] = 211 // 30 + 29 // 18 * 21^(1 // 2)
    A[7, 5] = -36 // 5 - 23 // 14 * 21^(1 // 2)
    A[7, 6] = 9 // 5 + 13 // 35 * 21^(1 // 2)
    A[8, 1] = 1 // 14
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 1 // 9 + 1 // 42 * 21^(1 // 2)
    A[8, 6] = 13 // 63 + 1 // 21 * 21^(1 // 2)
    A[8, 7] = 1 // 9
    A[9, 1] = 1 // 32
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 91 // 576 + 7 // 192 * 21^(1 // 2)
    A[9, 6] = 11 // 72
    A[9, 7] = -385 // 1152 + 25 // 384 * 21^(1 // 2)
    A[9, 8] = 63 // 128 - 13 // 128 * 21^(1 // 2)
    A[10, 1] = 1 // 14
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 1 // 9
    A[10, 6] = -733 // 2205 + 1 // 15 * 21^(1 // 2)
    A[10, 7] = 515 // 504 - 37 // 168 * 21^(1 // 2)
    A[10, 8] = -51 // 56 + 11 // 56 * 21^(1 // 2)
    A[10, 9] = 132 // 245 - 4 // 35 * 21^(1 // 2)
    A[11, 1] = 0
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = -7 // 3 - 7 // 18 * 21^(1 // 2)
    A[11, 6] = -2 // 5 - 28 // 45 * 21^(1 // 2)
    A[11, 7] = -91 // 24 + 53 // 72 * 21^(1 // 2)
    A[11, 8] = 301 // 72 - 53 // 72 * 21^(1 // 2)
    A[11, 9] = 28 // 45 + 28 // 45 * 21^(1 // 2)
    A[11, 10] = 49 // 18 + 7 // 18 * 21^(1 // 2)
    α[1] = 1 // 20
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = 49 // 180
    α[9] = 16 // 45
    α[10] = 49 // 180
    α[11] = 1 // 20

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 8, stability_size = 4.142567934840205))
end

"""
Some Explicit Runge-Kutta Methods of High Order, by G. J. Cooper and J. H. Verner,
 SIAM Journal on Numerical Analysis, Vol. 9, No. 3, (September 1972), pages 389 to 405
"""
function constructCooperVerner82(T::Type = Float64)
    A = zeros(T, 11, 11)
    c = zeros(T, 11)
    α = zeros(T, 11)

    c[2] = 1 // 2
    c[3] = 1 // 2
    c[4] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[5] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[6] = 1 // 2
    c[7] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[8] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[9] = 1 // 2
    c[10] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[11] = 1
    A[2, 1] = 1 // 2
    A[3, 1] = 1 // 4
    A[3, 2] = 1 // 4
    A[4, 1] = 1 // 7
    A[4, 2] = -1 // 14 - 3 // 98 * 21^(1 // 2)
    A[4, 3] = 3 // 7 + 5 // 49 * 21^(1 // 2)
    A[5, 1] = 11 // 84 + 1 // 84 * 21^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = 2 // 7 + 4 // 63 * 21^(1 // 2)
    A[5, 4] = 1 // 12 - 1 // 252 * 21^(1 // 2)
    A[6, 1] = 5 // 48 + 1 // 48 * 21^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 1 // 4 + 1 // 36 * 21^(1 // 2)
    A[6, 4] = -77 // 120 + 7 // 180 * 21^(1 // 2)
    A[6, 5] = 63 // 80 - 7 // 80 * 21^(1 // 2)
    A[7, 1] = 5 // 21 - 1 // 42 * 21^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = -48 // 35 + 92 // 315 * 21^(1 // 2)
    A[7, 4] = 211 // 30 - 29 // 18 * 21^(1 // 2)
    A[7, 5] = -36 // 5 + 23 // 14 * 21^(1 // 2)
    A[7, 6] = 9 // 5 - 13 // 35 * 21^(1 // 2)
    A[8, 1] = 1 // 14
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 1 // 9 - 1 // 42 * 21^(1 // 2)
    A[8, 6] = 13 // 63 - 1 // 21 * 21^(1 // 2)
    A[8, 7] = 1 // 9
    A[9, 1] = 1 // 32
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 91 // 576 - 7 // 192 * 21^(1 // 2)
    A[9, 6] = 11 // 72
    A[9, 7] = -385 // 1152 - 25 // 384 * 21^(1 // 2)
    A[9, 8] = 63 // 128 + 13 // 128 * 21^(1 // 2)
    A[10, 1] = 1 // 14
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 1 // 9
    A[10, 6] = -733 // 2205 - 1 // 15 * 21^(1 // 2)
    A[10, 7] = 515 // 504 + 37 // 168 * 21^(1 // 2)
    A[10, 8] = -51 // 56 - 11 // 56 * 21^(1 // 2)
    A[10, 9] = 132 // 245 + 4 // 35 * 21^(1 // 2)
    A[11, 1] = 0
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = -7 // 3 + 7 // 18 * 21^(1 // 2)
    A[11, 6] = -2 // 5 + 28 // 45 * 21^(1 // 2)
    A[11, 7] = -91 // 24 - 53 // 72 * 21^(1 // 2)
    A[11, 8] = 301 // 72 + 53 // 72 * 21^(1 // 2)
    A[11, 9] = 28 // 45 - 28 // 45 * 21^(1 // 2)
    A[11, 10] = 49 // 18 - 7 // 18 * 21^(1 // 2)
    α[1] = 1 // 20
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = 49 // 180
    α[9] = 16 // 45
    α[10] = 49 // 180
    α[11] = 1 // 20

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 8, stability_size = 3.715376401690181))
end

"""
An Eighth Order Runge-Kutta process with Eleven Function Evaluations per Step, by A. R. Curtis,
 Numerische Mathematik, Vol. 16, No. 3 (1970), pages 268 to 277
"""
function constructCurtis8(T::Type = Float64)
    A = zeros(T, 11, 11)
    c = zeros(T, 11)
    α = zeros(T, 11)
    c[2] = 1 // 192
    c[3] = 1295 // 10302 - 267 // 24038 * 21^(1 // 2)
    c[4] = 1295 // 6868 - 801 // 48076 * 21^(1 // 2)
    c[5] = 13 // 19
    c[6] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[7] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[8] = 1 // 2
    c[9] = 1 // 2 - 1 // 14 * 21^(1 // 2)
    c[10] = 1 // 2 + 1 // 14 * 21^(1 // 2)
    c[11] = 1
    A[2, 1] = 1 // 192
    A[3, 1] = -203059043 // 123819738 + 10606041 // 41273246 * 21^(1 // 2)
    A[3, 2] = 109311824 // 61909869 - 790320 // 2948089 * 21^(1 // 2)
    A[4, 1] = 1295 // 27472 - 801 // 192304 * 21^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = 3885 // 27472 - 2403 // 192304 * 21^(1 // 2)
    A[5, 1] = 250635205301 // 56025436876 + 45433220625 // 56025436876 * 21^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -886475583357 // 56025436876 - 158352084189 // 56025436876 * 21^(1 // 2)
    A[5, 4] = 168543392927 // 14006359219 + 28229715891 // 14006359219 * 21^(1 // 2)
    A[6, 1] = 91309 // 1560468 - 2059 // 520156 * 21^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = 4826396191 // 17845079061 - 1256382753 // 41638517809 * 21^(1 // 2)
    A[6, 5] = 777405919 // 4545562476 - 8073043 // 216455356 * 21^(1 // 2)
    A[7, 1] = 69085 // 390117 + 331 // 18577 * 21^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = -4774533908353 // 9725568088245 - 416000886501 // 3241856029415 * 21^(1 // 2)
    A[7, 5] = 130436567291 // 840929058060 + 11894034143 // 280309686020 * 21^(1 // 2)
    A[7, 6] = 53133 // 80660 + 15753 // 112924 * 21^(1 // 2)
    A[8, 1] = 246163 // 3566784 + 32841 // 1188928 * 21^(1 // 2)
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = -40141172781731 // 77804544705960 -
              180700440759 // 3241856029415 * 21^(1 // 2)
    A[8, 5] = 15546835747 // 96106178064 - 510947487 // 32035392688 * 21^(1 // 2)
    A[8, 6] = 294903 // 645280 + 2464 // 20165 * 21^(1 // 2)
    A[8, 7] = 21 // 64 - 5 // 64 * 21^(1 // 2)
    A[9, 1] = 2738629 // 18725616 - 1862059 // 56176848 * 21^(1 // 2)
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 25953361719391 // 58353408529470 -
              20890037862259 // 1225421579118870 * 21^(1 // 2)
    A[9, 5] = -782725944593 // 2018229739344 + 64202866997 // 864955602576 * 21^(1 // 2)
    A[9, 6] = 750631 // 8711280 - 1748677 // 26133840 * 21^(1 // 2)
    A[9, 7] = -133 // 432 + 91 // 1296 * 21^(1 // 2)
    A[9, 8] = 14 // 27 - 8 // 81 * 21^(1 // 2)
    A[10, 1] = -15028217 // 20286084 - 555173 // 2898012 * 21^(1 // 2)
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 543944452315361 // 126432385147185 +
               115426942571351 // 126432385147185 * 21^(1 // 2)
    A[10, 5] = -67975591357 // 73865390235 - 10957712053 // 73865390235 * 21^(1 // 2)
    A[10, 6] = -355837 // 42510 - 324791 // 178542 * 21^(1 // 2)
    A[10, 7] = 35 // 156 + 35 // 468 * 21^(1 // 2)
    A[10, 8] = 56 // 39 + 32 // 117 * 21^(1 // 2)
    A[10, 9] = 297 // 65 + 63 // 65 * 21^(1 // 2)
    A[11, 1] = 134701 // 668772 + 10947 // 74308 * 21^(1 // 2)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = -80282345563462 // 29176704264735 -
               963735684048 // 3241856029415 * 21^(1 // 2)
    A[11, 5] = 15546835747 // 18019908387 - 170315829 // 2002212043 * 21^(1 // 2)
    A[11, 6] = 14338141 // 2177820 + 7985863 // 5081580 * 21^(1 // 2)
    A[11, 7] = 31 // 54 - 5 // 42 * 21^(1 // 2)
    A[11, 8] = -16 // 27
    A[11, 9] = -23 // 5 - 39 // 35 * 21^(1 // 2)
    A[11, 10] = 13 // 18 - 13 // 126 * 21^(1 // 2)
    α[1] = 1 // 20
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 13 // 180
    α[7] = 1 // 5
    α[8] = 16 // 45
    α[9] = 1 // 5
    α[10] = 13 // 180
    α[11] = 1 // 20
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 8, stability_size = 5.6599472616274795))
end

"""
The Relative Efficiency of Alternative Defect Control Schemes for High-Order Continuous Runge-Kutta Formulas
 W. H. Enright SIAM Journal on Numerical Analysis, Vol. 30, No. 5. (Oct., 1993), pp. 1419-1445.
"""
function constructEnrightVerner8(T::Type = Float64)
    A = zeros(T, 13, 13)
    c = zeros(T, 13)
    α = zeros(T, 13)
    αEEst = zeros(T, 13)

    c[2] = 139 // 2500
    c[3] = 473499 // 4616000
    c[4] = 1420497 // 9232000
    c[5] = 1923 // 5000
    c[6] = 923 // 2000
    c[7] = 769 // 5000
    c[8] = 8571 // 10000
    c[9] = 14280181565693441800 // 15023510624554266963
    c[10] = 3611 // 5000
    c[11] = 15 // 16
    c[12] = 1
    c[13] = 1
    A[2, 1] = 139 // 2500
    A[3, 1] = 94226774499 // 11846945536000
    A[3, 2] = 224201303001 // 2369389107200
    A[4, 1] = 1420497 // 36928000
    A[4, 2] = 0
    A[4, 3] = 4261491 // 36928000
    A[5, 1] = 3590880889914277 // 9341720958375000
    A[5, 2] = 0
    A[5, 3] = -1122077307242443 // 778476746531250
    A[5, 4] = 1683359084698258 // 1167715119796875
    A[6, 1] = 3279216863 // 71027928000
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = 1418466643672 // 6146756958375
    A[6, 5] = 725783021041 // 3932390759616
    A[7, 1] = 22628776819090891 // 378192197956950000
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = 116326878837670912 // 1047320101758291075
    A[7, 5] = -21022506263989 // 614436056190000
    A[7, 6] = 182198401 // 10649112500
    A[8, 1] = -150702609045465151 // 280142368857000000
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = -47048174572430533795112 // 6781573975400382234375
    A[8, 5] = -29670244019740727 // 6363654238520000
    A[8, 6] = 503350581600824913 // 125990342179812500
    A[8, 7] = 9
    A[9, 1] = -BigInt(16045918035544526531085666708200976799092954349059102554439418058855237034803138076235914485211856830024864166368323768200) //
              BigInt(9829483157728664844198544034921655975378257830080951893693520116322068439377446708493269493633244146421522730587295739613)
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = -BigInt(2358539147881877873039445811751103289867403691705393271280994482043263320271984518766662589229245453919623760032009543680000000) //
              BigInt(217835526181955565558022800341429464086385172122906608031349543583863022015887117168566338470163090032498196064101837445053273)
    A[9, 5] = -BigInt(991136638972626168678903371416456100093900405535164924683058122802429707354033382826947398158683765324439618282500000000) //
              BigInt(79848142008846002789298925227605775190331269194726743910364273272231784282184770467794155269096224513726772081370189773)
    A[9, 6] = BigInt(99411279821210413387149352497211785829547149358696952646781033905129048593757052549024957474512389085654050445280000000) //
              BigInt(10219750071154355529199360565479179548497227516418867455143546824991883616697162903707552409044338358355837306818391433)
    A[9, 7] = BigInt(194327599672380134095898291719912961363678073793023525007081328425098431574448809779310732532821200046895000000000) //
              BigInt(11996011488227722649656673931136490891256463292620473601841875115170259043531987881275965965525693289412424400177)
    A[9, 8] = -BigInt(4738143867684122189593816244199450540483384372163549951990525387550768038015218275414120082248510000000000) //
              BigInt(45627381556119209828916528385326434273376137158228892503158792567081801761628344925618992749059885819540261)
    A[10, 1] = BigInt(100509763879264306824096153463041174636629364248095333923106653001873) //
               BigInt(229490324007644628042361756217436155832461488260089524115475000000000)
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = BigInt(35076261889213578261995286390053983221937920015616) //
               BigInt(8903662403052468890234321680409039895089390971875)
    A[10, 5] = BigInt(29877053472248545227782869189767925950557009) //
               BigInt(10443709158645362958089905740134206606110000)
    A[10, 6] = -BigInt(72602025182798889442893966553844286012776770019588838776297031451) //
               BigInt(40918439380405007926071673834718276106178863851548244800289187500)
    A[10, 7] = -BigInt(15322912063864370512130145988492098605486502994107816190) //
               BigInt(3130069400645604876183669549605539076333579290524919889)
    A[10, 8] = BigInt(66085154677219418645471125072555541174985695924) //
               BigInt(310222736648235062495097951962638800384203417327)
    A[10, 9] = -BigInt(33475654618965607625490266678231230366345830527265525310030016230875755239420324728600957368877132012320553021) //
               BigInt(563607181486505082775419237614873016043169125130359330269370345097328575740250457475506596993780051575000000000)
    A[11, 1] = -BigInt(5101097760197841615137571256611109219669965728737004453654940399435749927779) //
               BigInt(3460254790065218498025394912113654440547987641021174968985010753183759952128)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = -BigInt(739282853758412257967453242147288028248514875000) //
               BigInt(67244182875654829098323943713374207709631980757)
    A[11, 5] = -BigInt(304301954438407952266341438991435702455078125) //
               BigInt(26817588459995310678515661957421952038871616)
    A[11, 6] = BigInt(7944893992399893116476631697520363855501654462308622798446100145004371780765625) //
               BigInt(887005161602928659252827604035768657298435884939489754462273147402132735951354)
    A[11, 7] = BigInt(271466662889835128614810796916069488062485834784939879774912795880859375) //
               BigInt(17080057868208807584911919572395871967414618088449127584294090356526117)
    A[11, 8] = -BigInt(547671639173493535187249719091874257724281347519101322819710937500) //
               BigInt(5545897345777551056363997101744668576293320106258184122878184564689)
    A[11, 9] = BigInt(278805054229456473051914785770130401056277521071386274820253852826023602554035252878888463451556710568650135086751734733498388408614837563) //
               BigInt(57032789758993842471673989706725819315996341861372243823982103984407642346136141015055914541200747733442561186185728215790313661884647470336)
    A[11, 10] = -BigInt(7641031089620713267943821483116886435546875) //
                BigInt(1865115549729059107236140767953375513727732867)
    A[12, 1] = -BigInt(2173296165244568434534168496725754283210370856714048295955495704392191998074219) //
               BigInt(826329709817197468253912814996276721297719271873047121496939265708086637102000)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = -BigInt(1176485889334472345397948024769479865991267657667189808793600000) //
               BigInt(128238274125745784857395023683055391720808636521055761250085619)
    A[12, 5] = -BigInt(901641526273088667293677983960434306145731566007759375000) //
               BigInt(47006051321531726477650870675115597779536141098588395119)
    A[12, 6] = BigInt(44619197709345843038609995810290222092713734357503573021289399188441075277893825000) //
               BigInt(3047226829867759667805678237083584916028456629069014496162891378892137605913688013)
    A[12, 7] = BigInt(948339742665210716931560767459705347432692833888702670755668800833857753125000) //
               BigInt(54100203068523816523920547203712158257383235486719126879290708372319555115953)
    A[12, 8] = -BigInt(469359986530784242367841163967997985280992052062479522694037718750000) //
               BigInt(1262000068797740029604893636717124646814536106030908513651797310900389)
    A[12, 9] = -BigInt(664127481006260866245417353167879111725438608160134129788715555378465573117813001757238260693051901812594341608922016461968723754067215934799054852630446938063) //
               BigInt(947405313675779327886462044861834846879473861539583479995101168291116168844040438751584392785913224239894029843821464363224405440877235651675827150051218506000)
    A[12, 10] = BigInt(124594688949431236718278718078688571245967961641584151713065625000) //
                BigInt(2442266119995123914599926487344630873565880304003136401896150897137)
    A[12, 11] = BigInt(561064742310110985699639624014768981793940948761680464947352) //
                BigInt(671379391526329055316725343694451428305592235340194400714475)
    A[13, 1] = BigInt(2112537365439605317842214400294939600845864356126794907105288476157299972777) //
               BigInt(9791129992576348954381966935712572673760150144277222568673000574989723295600)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = BigInt(63401732568650174882469655273037557867024014649473659852800000) //
               BigInt(7597437179564699900056826682790095876501867794005624909986891)
    A[13, 5] = BigInt(6086761737797404637560639819223305802258863489403125000) //
               BigInt(2784859079002788542258846760345051095481316496802932391)
    A[13, 6] = -BigInt(2559665601805637409198888850698428172565517218973577146475607855715703433900000) //
               BigInt(1517075781450111169555356439342115045941218379260955323225394329122945220575603)
    A[13, 7] = -BigInt(3562962724414944757891560164749973139486260772587248231213362868690625000) //
               BigInt(408976639178968899516222045518619849445955714595234573207901740578867341)
    A[13, 8] = BigInt(27274774093948890572327172846505411328624378015249298079900000000) //
               BigInt(1115922486918625141794316350997246426030217325603407120793811703663)
    A[13, 9] = BigInt(1517649980273287794810711480928689762969941359446589355244683243470868873980537462717945530891242892188019850453060256773054273287757818339) //
               BigInt(17931096351403549517324822951619492690590259109554169692111771494162805088297821242848492102923447532215742446890559805801281057084758250800)
    A[13, 10] = BigInt(36524642706569597276845612368244490741747970924877885921875000) //
                BigInt(67204508345515429491838411146746729660269255184823759889942481)
    A[13, 11] = 0
    A[13, 12] = 0
    α[1] = 55284231195707975647126708111723 // 1258814250475014554517603612114000
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = BigInt(101652048214282205518610445783893750000000) //
           BigInt(289586658278060310247144250081091360673509)
    α[7] = BigInt(6219731958882093270433753490048828125000000) //
           BigInt(25268792314562899182186819800834272764054341)
    α[8] = BigInt(471525882014932587321673707929687500000000) //
           BigInt(523728817391198795728810839649314044495553)
    α[9] = BigInt(63998419659074502960979467027044380533513499562179716145788153981258354882170183557294261050806789914954901252698438375102943237148428019253408419) //
           BigInt(14067383878224159980676935999500324267727723947516309536464061683130063437616452951326819824531063499969576146304645311203502485409773697867038000)
    α[10] = BigInt(594562755257530592552224345703125000000) //
            BigInt(123802720910327682301431417435953442122031)
    α[11] = -1886691133979705639959153870454656 // 397947585885835383951563312487675
    α[12] = -50061468875139778913910254637881 // 141186641936036184819986782313781
    α[13] = 0
    αEEst[1] = 3635543992294021475202312550589 // 83920950031667636967840240807600
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = BigInt(91472553308336221233020750122000000000) //
               BigInt(270389036674192633283981559366098375979)
    αEEst[7] = BigInt(31245710879106859500854236329453125000000) //
               BigInt(125747467177230198813995913261775760852127)
    αEEst[8] = BigInt(7580382785455138868782239796250000000000) //
               BigInt(33873008089978031564549954803189465564389)
    αEEst[9] = -BigInt(1264668207534772000389416909088295806823076185957298091511366356638145411319432807599860924173256812108555319505561297378438513) //
               BigInt(31458133682363879376359837328611310990099259233936346807079836410100132902013789877583819562913435723874810037875576750158649200)
    αEEst[10] = BigInt(2363202400377740499662213515625000000) //
                BigInt(19167474982246118950523520271861502109)
    αEEst[11] = 0
    αEEst[12] = 0
    αEEst[13] = 2965876353908674111604434604609 // 47062213978678728273328927437927

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.642563947370989))
end

"""
Jim Verner's "Maple" (dverk78)
"""
function constructdverk78(T::Type = Float64)
    A = zeros(T, 13, 13)
    c = zeros(T, 13)
    α = zeros(T, 13)
    αEEst = zeros(T, 13)

    c[2] = 1 // 16
    c[3] = 112 // 1065
    c[4] = 56 // 355
    c[5] = 39 // 100
    c[6] = 7 // 15
    c[7] = 39 // 250
    c[8] = 24 // 25
    c[9] = 14435868 // 16178861
    c[10] = 11 // 12
    c[11] = 19 // 20
    c[12] = 1
    c[13] = 1
    A[2, 1] = 1 // 16
    A[3, 1] = 18928 // 1134225
    A[3, 2] = 100352 // 1134225
    A[4, 1] = 14 // 355
    A[4, 2] = 0
    A[4, 3] = 42 // 355
    A[5, 1] = 94495479 // 250880000
    A[5, 2] = 0
    A[5, 3] = -352806597 // 250880000
    A[5, 4] = 178077159 // 125440000
    A[6, 1] = 12089 // 252720
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = 2505377 // 10685520
    A[6, 5] = 960400 // 5209191
    A[7, 1] = 21400899 // 350000000
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = 3064329829899 // 27126050000000
    A[7, 5] = -21643947 // 592609375
    A[7, 6] = 124391943 // 6756250000
    A[8, 1] = -15365458811 // 13609565775
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = -7 // 5
    A[8, 5] = -8339128164608 // 939060038475
    A[8, 6] = 341936800488 // 47951126225
    A[8, 7] = 1993321838240 // 380523459069
    A[9, 1] = -BigInt(1840911252282376584438157336464708426954728061551) //
              BigInt(2991923615171151921596253813483118262195533733898)
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = -BigInt(14764960804048657303638372252908780219281424435) //
              BigInt(2981692102565021975611711269209606363661854518)
    A[9, 5] = -BigInt(875325048502130441118613421785266742862694404520560000) //
              BigInt(170212030428894418395571677575961339495435011888324169)
    A[9, 6] = BigInt(7632051964154290925661849798370645637589377834346780) //
              BigInt(1734087257418811583049800347581865260479233950396659)
    A[9, 7] = BigInt(7519834791971137517048532179652347729899303513750000) //
              BigInt(1045677303502317596597890707812349832637339039997351)
    A[9, 8] = BigInt(1366042683489166351293315549358278750) //
              BigInt(144631418224267718165055326464180836641)
    A[10, 1] = -63077736705254280154824845013881 // 78369357853786633855112190394368
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = -31948346510820970247215 // 6956009216960026632192
    A[10, 5] = -3378604805394255292453489375 // 517042670569824692230499952
    A[10, 6] = 1001587844183325981198091450220795 // 184232684207722503701669953872896
    A[10, 7] = 187023075231349900768014890274453125 // 25224698849808178010752575653374848
    A[10, 8] = 1908158550070998850625 // 117087067039189929394176
    A[10, 9] = -BigInt(52956818288156668227044990077324877908565) //
               BigInt(2912779959477433986349822224412353951940608)
    A[11, 1] = -BigInt(10116106591826909534781157993685116703) //
               BigInt(9562819945036894030442231411871744000)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = -9623541317323077848129 // 3864449564977792573440
    A[11, 5] = -4823348333146829406881375 // 576413233634141239944816
    A[11, 6] = BigInt(6566119246514996884067001154977284529) //
               BigInt(970305487021846325473990863582315520)
    A[11, 7] = BigInt(2226455130519213549256016892506730559375) //
               BigInt(364880443159675255577435648380047355776)
    A[11, 8] = 39747262782380466933662225 // 1756032802431424164410720256
    A[11, 9] = BigInt(48175771419260955335244683805171548038966866545122229) //
               BigInt(1989786420513815146528880165952064118903852843612160000)
    A[11, 10] = -2378292068163246 // 47768728487211875
    A[12, 1] = -3218022174758599831659045535578571 // 1453396753634469525663775847094384
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 26290092604284231996745 // 5760876126062860430544
    A[12, 5] = -697069297560926452045586710000 // 41107967755245430594036502319
    A[12, 6] = BigInt(1827357820434213461438077550902273440) //
               BigInt(139381013914245317709567680839641697)
    A[12, 7] = BigInt(643504802814241550941949227194107500000) //
               BigInt(242124609118836550860494007545333945331)
    A[12, 8] = 162259938151380266113750 // 59091082835244183497007
    A[12, 9] = -BigInt(23028251632873523818545414856857015616678575554130463402) //
               BigInt(20013169183191444503443905240405603349978424504151629055)
    A[12, 10] = 7958341351371843889152 // 3284467988443203581305
    A[12, 11] = -507974327957860843878400 // 121555654819179042718967
    A[13, 1] = -549080624436801105208519835138333 // 353987109028707139687100885600400
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = 29116675312186033956481 // 5331818957833865866320
    A[13, 5] = -91153092961177216058210567600 // 7609267653017028089793994539
    A[13, 6] = 1540775569495234383390307262972464 // 164751248733597841985831445059895
    A[13, 7] = -BigInt(2027488254536386321212021357622300000) //
               BigInt(7563800784313191306927092359781550321)
    A[13, 8] = -707592954577756600025 // 2430671607007101253926
    A[13, 9] = -BigInt(315869406877370103440389763510384832076674110731) //
               BigInt(550845901016714794114516357859531754926300500125)
    A[13, 10] = 4305840920849725632512 // 5066413598663346018375
    A[13, 11] = 0
    A[13, 12] = 0
    α[1] = 4631674879841 // 103782082379976
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 14327219974204125 // 40489566827933216
    α[7] = 2720762324010009765625000 // 10917367480696813922225349
    α[8] = -498533005859375 // 95352091037424
    α[9] = BigInt(405932030463777247926705030596175437402459637909765779) //
           BigInt(78803919436321841083201886041201537229769115088303952)
    α[10] = -10290327637248 // 1082076946951
    α[11] = 863264105888000 // 85814662253313
    α[12] = -29746300739 // 247142463456
    α[13] = 0
    αEEst[1] = 325503096889 // 7282953149472
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 495437430316125 // 1396191959583904
    αEEst[7] = 760741072216796875000 // 3055518466469861159313
    αEEst[8] = 1186143278515625 // 444976424841312
    αEEst[9] = BigInt(13791579353894559147128282108092181066885426809) //
               BigInt(3515333032435874239702147176443897942760320040)
    αEEst[10] = -30733299644928 // 5410384734755
    αEEst[11] = 0
    αEEst[12] = 0
    αEEst[13] = -45884771325 // 82380821152

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.785500208507322))
end

"""
Cheap Error Estimation for Runge-Kutta methods, by Ch. Tsitouras and S.N. Papakostas,
 Siam Journal on Scientific Computing, Vol. 20, Issue 6, Nov 1999.
"""
function constructTsitourasPapakostas8(T::Type = Float64)
    A = zeros(T, 13, 13)
    c = zeros(T, 13)
    α = zeros(T, 13)
    αEEst = zeros(T, 13)

    c[2] = 9 // 142
    c[3] = 24514 // 238491
    c[4] = 12257 // 79497
    c[5] = 50 // 129
    c[6] = 34 // 73
    c[7] = 23 // 148
    c[8] = 142 // 141
    c[9] = 29104120198 // 33218531707
    c[10] = 83 // 91
    c[11] = 143 // 149
    c[12] = 1
    c[13] = 1
    A[2, 1] = 9 // 142
    A[3, 1] = 9950845450 // 511901613729
    A[3, 2] = 42666469916 // 511901613729
    A[4, 1] = 12257 // 317988
    A[4, 2] = 0
    A[4, 3] = 12257 // 105996
    A[5, 1] = 14131686489425 // 35833975601529
    A[5, 2] = 0
    A[5, 3] = -17700454220625 // 11944658533843
    A[5, 4] = 17619604667500 // 11944658533843
    A[6, 1] = 4418011 // 96055225
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = 9757757988 // 41995818067
    A[6, 5] = 19921512441 // 106300094275
    A[7, 1] = 480285555889619 // 7997789253816320
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = 52162106556696460538643 // 464887033284472899365888
    A[7, 5] = -16450769949384489 // 490009784398315520
    A[7, 6] = 864813381633887 // 51718297665732608
    A[8, 1] = -7001048352088587137143 // 4449830351030385148800
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = -2522765548294599044197935933 // 1915963195493758447677901792
    A[8, 5] = -318556182235222634647116091 // 27172411532484037214054400
    A[8, 6] = 1205563885850790193966807 // 132365727505912221222912
    A[8, 7] = 254 // 39
    A[9, 1] = -BigInt(20629396399716689122801179264428539394462855874226604463554767070845753369) //
              BigInt(42881759662770155956657513470012114488076017981128105307375594963152353200)
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = -BigInt(3315443074343659404779422149387397712986453181141168247590906370819301077749322753) //
              BigInt(498517112641608872807821838566847987729514312398278633788185835348997643218553476)
    A[9, 5] = -BigInt(273749409411654060286948141164828452109898379203526945684314474186724062841643) //
              BigInt(60427583951377552503967840653825589117167443021628367691773162913607984889600)
    A[9, 6] = BigInt(16656372518874512738268060504309924900437672263609245028809229865738327731797537) //
              BigInt(4276990138533930522782076385771016009097627930550149628282203007913538394549504)
    A[9, 7] = BigInt(42008080033354305590804322944084805264441066760038302359736803632) //
              BigInt(4865302423216534910074823287811605599629170295030631799935804001)
    A[9, 8] = BigInt(668459780930716338000066627236927417191947396177093524377824) //
              BigInt(71100452948884643779799087713322002499799983434685642170251833)
    A[10, 1] = -BigInt(1793603946322260900828212460706877142477132870159527) //
               BigInt(2313097568511990753781649719556084665131024900300800)
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = -BigInt(14776874123722838192315406145167687512425345723701) //
               BigInt(1847893530366076701102014146927696206329050105856)
    A[10, 5] = -BigInt(19587020919884661714856757105130246995757906603) //
               BigInt(2911893296942532038566182868170393629789388800)
    A[10, 6] = BigInt(6364380863259071677112259236455506477417699780613300364807) //
               BigInt(1150428174584942133406579091549443438814988349188394057728)
    A[10, 7] = BigInt(27725164402569748756040320433848245155581006369) //
               BigInt(2544159473655547770881695354241256106302348256)
    A[10, 8] = BigInt(10744247163960019876833255044784609639) //
               BigInt(534761804739901348825491947768304503296)
    A[10, 9] = -BigInt(50977737930792808232204417497248979399878217280011103197862899) //
               BigInt(1300915694564913675613280314081837358644964393191337994183389184)
    A[11, 1] = -BigInt(3587625717068952487214493441966897048737050755812600710793) //
               BigInt(3015733164033229624772006086429467685046639983706974412800)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = -BigInt(5453011711267804731211501837262816944661201619903) //
               BigInt(764973320899716397072448644710733101329290196992)
    A[11, 5] = -BigInt(884348836774584715070440485633026464653487653) //
               BigInt(92725983515963799584249219499787659014963200)
    A[11, 6] = BigInt(26823469063654084387375587616552322383082061411417182757389742951) //
               BigInt(3541299744763681675473620647087123057228744296123642301654630400)
    A[11, 7] = BigInt(142363419491686507162007051071007722765323162710521029) //
               BigInt(12634887202368261565807449771335728082273587905775840)
    A[11, 8] = BigInt(64747617454909275289531520412519442831235890581) //
               BigInt(1269317188117975960996670628974453443777854830080)
    A[11, 9] = BigInt(112633808253272720979874303367503891597499261046700689572459050065039333987335667) //
               BigInt(1404514291245034532812181377119034501014039830518218465064579291611563374608650240)
    A[11, 10] = -BigInt(10612202518573994431153697720606405883) //
                BigInt(67082546658259846778754594976831647575)
    A[12, 1] = -BigInt(7534081165544982478296202335922049210803875045423) //
               BigInt(19219575665440848756074598051658416387002479820800)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = BigInt(237696087452786717802270375283034262859273455) //
               BigInt(60688480889936839261131818793519097177034752)
    A[12, 5] = -BigInt(20610578209826329263318986584876108069323) //
               BigInt(7356333776438834884293014575840932659200)
    A[12, 6] = BigInt(51260471529841028040709654458903254781320136131844164563) //
               BigInt(20998023776318546907382302106788168158765514131585630208)
    A[12, 7] = -BigInt(3077214437173472971196810795615384000211457151011) //
               BigInt(1272435592582280820597059432060116893200684680384)
    A[12, 8] = -BigInt(1539218116260541896259682954580256454049) //
               BigInt(4534670830750360983422999946409310393344)
    A[12, 9] = BigInt(241886539350268429372116296787271276553970618941104594460614948326132797451456131) //
               BigInt(1240669632662465892528916120041123051054026283188324780101555345343662316090597376)
    A[12, 10] = -80556486832245966191717452425924975 // 414445409518676597565032008051106461
    A[12, 11] = BigInt(2944781680874500347594142792814463350) //
                BigInt(4965161383073676983610218096030654529)
    A[13, 1] = -BigInt(7757739937862944832927743694336116203639371542761) //
               BigInt(5225100678421794325654850845451340473577260032000)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = -BigInt(433889546009521405913741133329446636837810749) //
               BigInt(181488796115643001621310691169322233346938880)
    A[13, 5] = -BigInt(246044720162308748108107126829066792329071) //
               BigInt(21999103311417807170100413255475150848000)
    A[13, 6] = BigInt(2140331235425829844389060818616719848637810765257179167) //
               BigInt(245428182036011897638028252815369258646052549878743040)
    A[13, 7] = BigInt(1573990926219809229258666534611598771063240529) //
               BigInt(214535514317495538101650508596881057760818080)
    A[13, 8] = BigInt(62408280667309375445301959066100433563) //
               BigInt(4838320046251983849512355559327451645440)
    A[13, 9] = BigInt(1145609822249493677618113725359506642998153205603226883141207089968379) //
               BigInt(26902802166822768476367427840835743139559294105398677975568538466119680)
    A[13, 10] = -408950356875874683139089678053832 // 7674292714443204070455739595109785
    A[13, 11] = 0
    A[13, 12] = 0
    α[1] = 55038446513529253801 // 1239280570055853383520
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 2335496795323782464411846394611 // 6598368783294020109895379936256
    α[7] = 7636073376527143565375240869888 // 30725949199261296642046754748645
    α[8] = -4237087214169934312729607487 // 12735791394214625116604076160
    α[9] = BigInt(408505291291133241760995514121984335914363927884426780078325258228227984174126699) //
           BigInt(212624874612193697466655159405635202123821166475348052734199905546455860773575680)
    α[10] = -1108225296327029096435947 // 405679075893979729103310
    α[11] = 2460988291206213825688467985 // 1756342789520947764222671739
    α[12] = 4808707937311 // 50545545388065
    α[13] = 0
    αEEst[1] = 715953338020208413 // 16094552857868225760
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 62284335162928966987066121 // 175437206755843240272669696
    αEEst[7] = 184309146777472302831695872 // 742417767522160213324368465
    αEEst[8] = -509771598215811385123057257 // 210282269969085698264101760
    αEEst[9] = BigInt(2701602489646143640362891402924962500766379885231830470264478480621389) //
               BigInt(1688575269294196295712420798364925687791337062814161130291144634516480)
    αEEst[10] = -7218534073012286740367561 // 3986456306153224954098780
    αEEst[11] = 0
    αEEst[12] = 0
    αEEst[13] = 5752173075461 // 1925544586212

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.906883483660496))
end

"""
From Verner's Webiste
"""
function constructVernerRobust9(T::Type = Float64)
    A = zeros(T, 16, 16)
    c = zeros(T, 16)
    α = zeros(T, 16)
    αEEst = zeros(T, 16)

    c[2] = 1 // 25
    c[3] = BigInt(48) // BigInt(335) - BigInt(32) // BigInt(1675) * 6^(1 // 2)
    c[4] = BigInt(72) // BigInt(335) - BigInt(48) // BigInt(1675) * 6^(1 // 2)
    c[5] = 72 // 125
    c[6] = BigInt(48) // BigInt(125) - BigInt(8) // BigInt(125) * 6^(1 // 2)
    c[7] = BigInt(48) // BigInt(125) + BigInt(8) // BigInt(125) * 6^(1 // 2)
    c[8] = 16 // 25
    c[9] = 12 // 25
    c[10] = 3377 // 50000
    c[11] = 1 // 4
    c[12] = 57617028499878 // 85094827871699
    c[13] = 1623 // 2000
    c[14] = 453 // 500
    c[15] = 1
    c[16] = 1
    A[2, 1] = 1 // 25
    A[3, 1] = -BigInt(15792) // BigInt(112225) + BigInt(5536) // BigInt(112225) * 6^(1 // 2)
    A[3, 2] = BigInt(31872) // BigInt(112225) - BigInt(1536) // BigInt(22445) * 6^(1 // 2)
    A[4, 1] = BigInt(18) // BigInt(335) - BigInt(12) // BigInt(1675) * 6^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = BigInt(54) // BigInt(335) - BigInt(36) // BigInt(1675) * 6^(1 // 2)
    A[5, 1] = BigInt(4014) // BigInt(3125) + BigInt(252) // BigInt(625) * 6^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -BigInt(14742) // BigInt(3125) - BigInt(972) // BigInt(625) * 6^(1 // 2)
    A[5, 4] = BigInt(12528) // BigInt(3125) + BigInt(144) // BigInt(125) * 6^(1 // 2)
    A[6, 1] = BigInt(1232) // BigInt(16875) - BigInt(152) // BigInt(16875) * 6^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = BigInt(29684) // BigInt(106875) - BigInt(13372) // BigInt(320625) * 6^(1 // 2)
    A[6, 5] = BigInt(2132) // BigInt(64125) - BigInt(284) // BigInt(21375) * 6^(1 // 2)
    A[7, 1] = BigInt(2032) // BigInt(16875) + BigInt(152) // BigInt(16875) * 6^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = -BigInt(7348) // BigInt(98325) - BigInt(33652) // BigInt(294975) * 6^(1 // 2)
    A[7, 5] = BigInt(10132) // BigInt(64125) - BigInt(716) // BigInt(21375) * 6^(1 // 2)
    A[7, 6] = BigInt(2592) // BigInt(14375) + BigInt(2912) // BigInt(14375) * 6^(1 // 2)
    A[8, 1] = 16 // 225
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 0
    A[8, 6] = BigInt(64) // BigInt(225) + BigInt(4) // BigInt(225) * 6^(1 // 2)
    A[8, 7] = BigInt(64) // BigInt(225) - BigInt(4) // BigInt(225) * 6^(1 // 2)
    A[9, 1] = 57 // 800
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 0
    A[9, 6] = BigInt(177) // BigInt(800) + BigInt(69) // BigInt(1600) * 6^(1 // 2)
    A[9, 7] = BigInt(177) // BigInt(800) - BigInt(69) // BigInt(1600) * 6^(1 // 2)
    A[9, 8] = -27 // 800
    A[10, 1] = 2844530829046074022657 // 58982400000000000000000
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 0
    A[10, 6] = BigInt(4287156859652598464203) // BigInt(58982400000000000000000) -
               BigInt(1598864762333658025459) // BigInt(117964800000000000000000) *
               6^(1 // 2)
    A[10, 7] = BigInt(4287156859652598464203) // BigInt(58982400000000000000000) +
               BigInt(1598864762333658025459) // BigInt(117964800000000000000000) *
               6^(1 // 2)
    A[10, 8] = -141033886218604337343 // 6553600000000000000000
    A[10, 9] = -21409264848554971927 // 204800000000000000000
    A[11, 1] = -BigInt(72189389771) // BigInt(9959178240000) -
               BigInt(459663572789) // BigInt(59755069440000) * 6^(1 // 2)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = 0
    A[11, 6] = 1 // 30
    A[11, 7] = -BigInt(14201240926266911) // BigInt(557169364500480000) -
               BigInt(31790792357660029) // BigInt(557169364500480000) * 6^(1 // 2)
    A[11, 8] = BigInt(22414436941) // BigInt(1563197440000) +
               BigInt(459663572789) // BigInt(56275107840000) * 6^(1 // 2)
    A[11, 9] = BigInt(154180604903) // BigInt(2534154240000) +
               BigInt(459663572789) // BigInt(11403694080000) * 6^(1 // 2)
    A[11, 10] = BigInt(21871487332435000000) // BigInt(125536952879579583419) +
                BigInt(18386542911560000000) // BigInt(1129832575916216250771) * 6^(1 // 2)
    A[12, 1] = -BigInt(178144571353393080183496267158614821877982611914666395752937745405391408707734804982447062502773) //
               BigInt(1247718010112994054746145516410425353598134947600568397156373491324203879120405304413829778240000) +
               BigInt(352591194575569317115651180991223097026568880384478484650944931412648046608828531034562538513357) //
               BigInt(8110167065734461355849945856667764798387877159403694581516427693607325214282634478689893558560000) *
               6^(1 // 2)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 0
    A[12, 5] = 0
    A[12, 6] = -BigInt(39115022545645779688585831988975140882502245161831122703201435591036520210945751850583137867) //
               BigInt(38425428798268102367071108589483043166777031255688218331140907046999467997785653437677908480) -
               BigInt(490826700396287454540598331961129757400186839154291202145347161815801053957874498846640625) //
               BigInt(15370171519307240946828443435793217266710812502275287332456362818799787199114261375071163392) *
               6^(1 // 2)
    A[12, 7] = -BigInt(622064296680932516193525865473060264652463486635595667447855699920023546211607797417635623004817316459) //
               BigInt(907451703567059050763665262359998604536360756710836393258911136378360869155327677312519401429106240000) +
               BigInt(643210041535328932923955834959360270277930780334485030265105750796567439186095378077346484617022694073) //
               BigInt(1814903407134118101527330524719997209072721513421672786517822272756721738310655354625038802858212480000) *
               6^(1 // 2)
    A[12, 8] = BigInt(1203943546728385294644268186854769106596033156989737459970836017592684437096940490927391123081521) //
               BigInt(6546732431504927940789740125933173479539636700187880198168132038132534360122730704444373657280000) -
               BigInt(151110511960958278763850506139095613011386663021919350564690684891134877118069370443383945077153) //
               BigInt(3273366215752463970394870062966586739769818350093940099084066019066267180061365352222186828640000) *
               6^(1 // 2)
    A[12, 9] = BigInt(63681466156701378449525903160928961790100434377497936201427140303316098275608362465676059840449) //
               BigInt(119057634931893490852780569777780769114850167606078444507680127107913165342297705956314829715000) -
               BigInt(352591194575569317115651180991223097026568880384478484650944931412648046608828531034562538513357) //
               BigInt(1547749254114615381086147407111149998493052178879019778599841652402871149449870177432092786295000) *
               6^(1 // 2)
    A[12, 10] = BigInt(10278934048763239705668635930656276438279802537812076691952406255660602976082835598179622579836784640000000) //
                BigInt(20706607333558650004563853916074751077994608347553087587301432465002027701979512351053114927376343490460197) -
                BigInt(24757946279683404392669266925829425235785590869511266396518921812563538267024485653444025561440747520000000) //
                BigInt(269185895336262450059330100908971764013929908518190138634918622045026360125733660563690494055892465375982561) *
                6^(1 // 2)
    A[12, 11] = BigInt(19635000096509466380843956455829094932847113883632439882926745444809822658280742880787456) //
                BigInt(15009933124323477487137151792766813737022277834253210285601916815234167186635020874092933)
    A[13, 1] = -BigInt(25488511950948766602163761842966272037005387677568247) //
               BigInt(25343364340644945281003771622773961523200000000000000) +
               BigInt(2525608241949563386308964624438617443527) //
               BigInt(12416017897756354830807859200000000000000) * 6^(1 // 2)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = 0
    A[13, 5] = 0
    A[13, 6] = -BigInt(2617546081675469247418718340204655213) //
               BigInt(431392587109884206933606400000000000) -
               552470350996365859393640759989 // 2400793528263703412500070400000 *
               6^(1 // 2)
    A[13, 7] = -BigInt(45969294618407232267578352626581642155421231201) //
               BigInt(10187731154133781589406192186163200000000000000) +
               BigInt(1540493536582818303738906021531546510663696889) //
               BigInt(885889665576850572991842798796800000000000000) * 6^(1 // 2)
    A[13, 8] = BigInt(203013873418014401588800777489550153682185921397719) //
               BigInt(6003783517480134311881093736968119910400000000000000) -
               BigInt(3968812951635028178485515838403541696971) //
               BigInt(18374628007211630439078297600000000000000) * 6^(1 // 2)
    A[13, 9] = BigInt(57340072791914637839492204156400815228449519011392389) //
               BigInt(19744238839693126186283635091213610188800000000000000) -
               BigInt(677602211254760908521917338264019314117) //
               BigInt(635714271434389910412902400000000000000) * 6^(1 // 2)
    A[13, 10] = BigInt(24100788715039192225758197856261740786803059875758180717) //
                BigInt(9662717001210818114953908042643919942065136939286292000) -
                BigInt(8800028717594297513271653743688562521) //
                BigInt(20449362746897197076153911415852248500) * 6^(1 // 2)
    A[13, 11] = BigInt(124358916033523439225154730110040589545064737) //
                BigInt(19935912449945669485002510178704400000000000)
    A[13, 12] = BigInt(594716139297486674475082356103592029886330739685357945223223985387172161751824241784555387316314439) //
                BigInt(818019847119678359631725769101238313386585208743795596694052943582031156000597320370585600000000000)
    A[14, 1] = BigInt(184592679361470753674594239208189459317466722387929369596771) //
               BigInt(127891217091830656233364014205731655064571321600000000000000) -
               BigInt(293880952873230803935166416296767553377) //
               BigInt(866938749696854853908947200000000000000) * 6^(1 // 2)
    A[14, 2] = 0
    A[14, 3] = 0
    A[14, 4] = 0
    A[14, 5] = 0
    A[14, 6] = BigInt(5966995986367380181027852718263477) //
               BigInt(684582962942931480729600000000000) +
               BigInt(2897291244884828193281565089) //
               BigInt(19049265055803310768128000000) * 6^(1 // 2)
    A[14, 7] = BigInt(99009240614446727611594225452277068906413179) //
               BigInt(16167053833464253010532287404800000000000000) -
               BigInt(3748619435624523040662096435375495632409181) //
               BigInt(1405830768127326348741938035200000000000000) * 6^(1 // 2)
    A[14, 8] = BigInt(4632645823113234486269357844283603936661991916850927460433) //
               BigInt(5096189922215990577140952099757378096638035200000000000000) +
               BigInt(41982993267604400562166630899538221911) //
               BigInt(116635822311401951029305600000000000000) * 6^(1 // 2)
    A[14, 9] = -BigInt(12135475677184688070559819501757154968443992507349800863475059) //
               BigInt(3337534088557268857775113382843672473526974754800000000000000) +
               BigInt(293880952873230803935166416296767553377) //
               BigInt(165447122399672764770545400000000000000) * 6^(1 // 2)
    A[14, 10] = -BigInt(8031428603597180964147750321204281053557209400385686643427859804352) //
                BigInt(2513586440178391306900178164909088113749933981949109368563389778875) +
                BigInt(21495292553013453087829315020563569618432) //
                BigInt(29973653446264566614372595657785433238875) * 6^(1 // 2)
    A[14, 11] = -BigInt(44866632128961269158599595825347949153517725172080973) //
                BigInt(5429607175821308141232283339099249860726338769531250)
    A[14, 12] = -BigInt(359988533934752543813368225945058137961504440675167771570563740579175484411998973149563531328527853237660249760109643) //
                BigInt(234128137445668240736067711662442379892483470931740696615716602345729933343308008758732590194294958617150200000000000)
    A[14, 13] = 1011136807359189181222571916914688 // 2927578889427871359661497484969729
    A[15, 1] = -BigInt(43117938494612449223384237106139955697243762787772847383) //
               BigInt(14030104862294690648579566521727488445420965260989440000) +
               BigInt(4220334940168717525563719481764309232553) //
               BigInt(5553268895030690235789996003491120640000) * 6^(1 // 2)
    A[15, 2] = 0
    A[15, 3] = 0
    A[15, 4] = 0
    A[15, 5] = 0
    A[15, 6] = -BigInt(492170504431044248358505817476691843) //
               BigInt(26311016381549020957251979880325120) -
               BigInt(65372649360291914372644744384375) //
               BigInt(457582893592156886213077910962176) * 6^(1 // 2)
    A[15, 7] = -BigInt(60242025722206083647640080390153021797730167) //
               BigInt(4671870346091448695506240355969515576320000) +
               BigInt(312292306982904450161582808948983992497619279) //
               BigInt(54031196176535884913246084116864832317440000) * 6^(1 // 2)
    A[15, 8] = -BigInt(10128486019754425336362855829794634260766582925227304979) //
               BigInt(2051440458733049774272501566471676463052668060718080000) -
               BigInt(4220334940168717525563719481764309232553) //
               BigInt(5229862652007483519857107084128791040000) * 6^(1 // 2)
    A[15, 9] = BigInt(90164396775224825317874747969895939573106537263617688583) //
               BigInt(9932194582195572820411492715910530662743409752296080000) -
               BigInt(4220334940168717525563719481764309232553) //
               BigInt(1059789239915401287502749181781142480000) * 6^(1 // 2)
    A[15, 10] = BigInt(822220421784353922416106346105112190310210366765138782108892460117760000000) //
                BigInt(112539961436930436130682455392076582934435837857314380718907401253008522853) -
                BigInt(43216229787327667461772487493266526541342720000000) //
                BigInt(26879922067906624005297106621487120800647611685231) * 6^(1 // 2)
    A[15, 11] = BigInt(209814213871916216679569640811090009729840605120) //
                BigInt(11651974660960738096987230622554824075229159217)
    A[15, 12] = BigInt(3941550892863952281184122715485813150214519858148244578543766661999079625453432169420046111603168071904019820197053706993310801913) //
                BigInt(576348806945747835561749398160966945390916610765210182493857593996732757363976976471505530337022074730630638475548489174498085584)
    A[15, 13] = -BigInt(4759443892077050695292772766912275498897121280000000000) //
                BigInt(4614026158486489304426198097242323507384603506591241137)
    A[15, 14] = BigInt(38733614144315448443538220691259283298308750000000) //
                BigInt(93806291160324220065540465392781147766752569954001)
    A[16, 1] = BigInt(18092408213832965447840389945226814727838921316097789) //
               BigInt(4568260908663202892681034746459555501706277524480000) -
               BigInt(156317653995567850184921063401246768027) //
               BigInt(212359240683067701291661898043594240000) * 6^(1 // 2)
    A[16, 2] = 0
    A[16, 3] = 0
    A[16, 4] = 0
    A[16, 5] = 0
    A[16, 6] = BigInt(65946837846509057750066287355594411) //
               BigInt(3018431611722608724311451319541760) +
               BigInt(41984484090942118562371132909375) //
               BigInt(52494462812567108248894805557248) * 6^(1 // 2)
    A[16, 7] = BigInt(8692965898134239206513518248119006769975359) //
               BigInt(535962614823230784958890847664050959360000) -
               BigInt(38800907634871851557429961190168666395776983) //
               BigInt(6198524154042582121698476759940763269120000) * 6^(1 // 2)
    A[16, 8] = BigInt(2641630262222426284468500071736583566105608571199325883) //
               BigInt(235343729805469979809271071241826755571130067491840000) +
               BigInt(156317653995567850184921063401246768027) //
               BigInt(199992055607259790546219146109360640000) * 6^(1 // 2)
    A[16, 9] = -BigInt(193422103474923353694542105476136625982637168991506421) //
               BigInt(16048357343710211432072833749116213173856151273040000) +
               BigInt(156317653995567850184921063401246768027) //
               BigInt(40526767661819895652001992553977680000) * 6^(1 // 2)
    A[16, 10] = -BigInt(19419587503350096681694941672242361995243037384025897671049255680000000) //
                BigInt(2155738929233095407456284002740551925607350553898010429805538450410571) +
                BigInt(1600692776914614785893591689228766904596480000000) //
                BigInt(1027899053306897212338272667278244239283200723671) * 6^(1 // 2)
    A[16, 11] = -BigInt(27285760059113968837908932750404506400064658240) //
                BigInt(1336728621411169422812695849204871581103756891)
    A[16, 12] = -BigInt(194853016061874333335484555436551917329750324944158872786342511110263341294405932910123138589192221578993200454377) //
                BigInt(15840605894584140753571150900860666508000461134947537749588341642048479030548925808450142813874420689113167632464)
    A[16, 13] = BigInt(89488735256675489275877211680603518179732480000000000) //
                BigInt(58814073298029244653337926984054935524805348324518339)
    A[16, 14] = 0
    A[16, 15] = 0
    α[1] = 100976787617015984669475787 // 6921502952403262310437464576
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = 8877148253451235588984375 // 4385514038208782482187638272
    α[9] = 961916572949681511747758515625 // 4416417715504587515036809762176
    α[10] = BigInt(1967337516701564001434375000000000000000000000000) //
            BigInt(15431364863119854943071851131903908575429289017877)
    α[11] = 2323713252076974806855457536 // 10352378514220126928031114081
    α[12] = BigInt(28308600293241456954311939117138937391141264054325158121088859798090056999504636993802634122671806566061266984631226363619) //
            BigInt(158391546540670080968065757007173890885877631132509858232987982213169380898029758612356596464561682203552186100842929847264)
    α[13] = BigInt(296136352341197653422080000000000000) //
            BigInt(3899432561650270968394778037550931439)
    α[14] = 117048651891177050452812500000000 // 903958175807874008864483503817601
    α[15] = 145778296653275182685983 // 4945417885871057962703934
    α[16] = 0
    αEEst[1] = 34542436255316150799031 // 1697695107285568386175488
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 0
    αEEst[7] = 0
    αEEst[8] = 11089416604965799654140625 // 10367645480399012960254464
    αEEst[9] = 19661377435148646805703125 // 255979697183364488207083392
    αEEst[10] = BigInt(688790393936688343750000000000000000000000) //
                BigInt(6091295374506475008386869675533556427693991)
    αEEst[11] = 2387113868151976968347648 // 9351742108599933990994683
    αEEst[12] = -BigInt(12428699133337309103655890990808226652772000718643045082918121983471100560973332724469319414156511978021429) //
                BigInt(12648919237095643063555422745116719205173888728920704277722434371356124846214898200346632677347950735102176)
    αEEst[13] = 185286960316915979617280000000000000 // 465364379785423284294027008596613217
    αEEst[14] = 0
    αEEst[15] = 0
    αEEst[16] = 16723862451391031122709 // 339047561915509883967882

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 9, αEEst = αEEst, adaptiveorder = 8,
        stability_size = 4.5214477514362255))
end

"""
From Verner's Webiste
"""
function constructVernerEfficient9(T::Type = Float64)
    A = zeros(T, 16, 16)
    c = zeros(T, 16)
    α = zeros(T, 16)
    αEEst = zeros(T, 16)

    c[2] = 1731 // 50000
    c[3] = BigInt(7630049) // BigInt(53810000) -
           BigInt(983539) // BigInt(53810000) * 6^(1 // 2)
    c[4] = BigInt(22890147) // BigInt(107620000) -
           BigInt(2950617) // BigInt(107620000) * 6^(1 // 2)
    c[5] = 561 // 1000
    c[6] = BigInt(387) // BigInt(1000) - BigInt(129) // BigInt(2000) * 6^(1 // 2)
    c[7] = BigInt(387) // BigInt(1000) + BigInt(129) // BigInt(2000) * 6^(1 // 2)
    c[8] = 129 // 200
    c[9] = 387 // 800
    c[10] = 6757 // 100000
    c[11] = 1 // 4
    c[12] = 1427971650951258372 // 2166662646162554701
    c[13] = 4103 // 5000
    c[14] = 2253 // 2500
    c[15] = 1
    c[16] = 1
    A[2, 1] = 1731 // 50000
    A[3, 1] = -BigInt(177968356965557) // BigInt(1002427673820000) +
              BigInt(14180534491313) // BigInt(250606918455000) * 6^(1 // 2)
    A[3, 2] = BigInt(64021741529527) // BigInt(200485534764000) -
              BigInt(7504450763411) // BigInt(100242767382000) * 6^(1 // 2)
    A[4, 1] = BigInt(22890147) // BigInt(430480000) -
              BigInt(2950617) // BigInt(430480000) * 6^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = BigInt(68670441) // BigInt(430480000) -
              BigInt(8851851) // BigInt(430480000) * 6^(1 // 2)
    A[5, 1] = BigInt(592203994261020339) // BigInt(513126355505556250) +
              BigInt(730386990293623641) // BigInt(2052505422022225000) * 6^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -BigInt(8712153884182794903) // BigInt(2052505422022225000) -
              BigInt(2843421359195851533) // BigInt(2052505422022225000) * 6^(1 // 2)
    A[5, 4] = BigInt(1873698362223295443) // BigInt(513126355505556250) +
              BigInt(528258592225556973) // BigInt(513126355505556250) * 6^(1 // 2)
    A[6, 1] = BigInt(11380823631) // BigInt(157617812000) -
              BigInt(339148869) // BigInt(39404453000) * 6^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = BigInt(16193232887091831) // BigInt(58864341808507450) -
              BigInt(2355345717024309) // BigInt(58864341808507450) * 6^(1 // 2)
    A[6, 5] = BigInt(165912282616977) // BigInt(4179075230308000) -
              BigInt(33181894472511) // BigInt(2089537615154000) * 6^(1 // 2)
    A[7, 1] = BigInt(26523528363) // BigInt(231790900000) +
              BigInt(863255358) // BigInt(123138915625) * 6^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = -BigInt(38208748178016484817787) // BigInt(842517966262441068418750) -
              BigInt(86118788556282369822807) // BigInt(842517966262441068418750) *
              6^(1 // 2)
    A[7, 5] = BigInt(92362336407446913) // BigInt(290322814529044000) -
              BigInt(232039320950012997) // BigInt(2467743923496874000) * 6^(1 // 2)
    A[7, 6] = -BigInt(362925891) // BigInt(1690350537500) +
              BigInt(857800423623) // BigInt(3380701075000) * 6^(1 // 2)
    A[8, 1] = 43 // 600
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 0
    A[8, 6] = BigInt(43) // BigInt(150) + BigInt(43) // BigInt(2400) * 6^(1 // 2)
    A[8, 7] = BigInt(43) // BigInt(150) - BigInt(43) // BigInt(2400) * 6^(1 // 2)
    A[9, 1] = 7353 // 102400
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 0
    A[9, 6] = BigInt(22833) // BigInt(102400) + BigInt(8901) // BigInt(204800) * 6^(1 // 2)
    A[9, 7] = BigInt(22833) // BigInt(102400) - BigInt(8901) // BigInt(204800) * 6^(1 // 2)
    A[9, 8] = -3483 // 102400
    A[10, 1] = 376708742472214988700853 // 7788456028125000000000000
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 0
    A[10, 6] = BigInt(187914666753956840195279) // BigInt(2596152009375000000000000) -
               BigInt(210440846556290693268911) // BigInt(15576912056250000000000000) *
               6^(1 // 2)
    A[10, 7] = BigInt(187914666753956840195279) // BigInt(2596152009375000000000000) +
               BigInt(210440846556290693268911) // BigInt(15576912056250000000000000) *
               6^(1 // 2)
    A[10, 8] = -18552667221896744226647 // 865384003125000000000000
    A[10, 9] = -3167799860072183913409 // 30423656359863281250000
    A[11, 1] = -BigInt(426968570497) // BigInt(54394415898750) -
               BigInt(92754382349) // BigInt(12087647977500) * 6^(1 // 2)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = 0
    A[11, 6] = 1 // 30
    A[11, 7] = -BigInt(2865012129681958) // BigInt(114898584332330625) -
               BigInt(12962517687655099) // BigInt(229797168664661250) * 6^(1 // 2)
    A[11, 8] = BigInt(4389715333607) // BigInt(309890657317500) +
               BigInt(92754382349) // BigInt(11477431752500) * 6^(1 // 2)
    A[11, 9] = BigInt(4990058173976) // BigInt(83757096376875) +
               BigInt(371017529396) // BigInt(9306344041875) * 6^(1 // 2)
    A[11, 10] = BigInt(1099523524595993125000) // BigInt(6257667909869756018891) +
                BigInt(100957348037989687500) // BigInt(6257667909869756018891) * 6^(1 // 2)
    A[12, 1] = BigInt(18382031104798403869938539009154656587521498573595595063164077882800315372787284683238439478955141517997198007108623761931447163756) //
               BigInt(13974256944499724344918960993890933614161025322970450047932688998095008528620821239604734608111291769444706187497807869179550841329375) +
               BigInt(407885778185158609210793892517582595305896470756467612636796259611491408260896413446883450891351622914818800693274034252252905536) //
               BigInt(28084926388601226073624096169175002956970191576455110633226765141161372294098693275117181239385312198137508846535933127837167926875) *
               6^(1 // 2)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 0
    A[12, 5] = 0
    A[12, 6] = -BigInt(333881311789849411971573472868128281438202210721723123251742145367734582887577395547778228760174068758086134389952015563403904) //
               BigInt(2270872004608103037127689848604039623086639035441372934050180593816493796129405349914148981460714202232988727738778494557727635) +
               BigInt(4819272892477768171373308666720689121421091953625792970278044071549950640195056472955523769829034800621890424847009130000000) //
               BigInt(23162894447002650978702436455761204155483718161502003927311842056928236720519934569124319610899284862776485022935540644488821877) *
               6^(1 // 2)
    A[12, 7] = -BigInt(136666607496463622270135608863772076443625468798139480390426740993024803946981763209348364716108721312822619845726151693667598437699964416) //
               BigInt(3719286465342404274788585327254180828195282427342057650194855634917821113563432870681372043512520401887141437067106105683944802332422369375) +
               BigInt(169845085565361336805556009296394374527636952379388961026066628725155521832762086875632366996477567928657535912191396155566765457826139904) //
               BigInt(1593979913718173260623679425966077497797978183146596135797795272107637620098614087434873732933937315094489187314474045293119200999609586875) *
               6^(1 // 2)
    A[12, 8] = BigInt(5610987899273278525411960528081442902198567594809764379756195673673265700551076812883925583370253765702553235594764427173637673766208) //
               BigInt(92881598198144033018278804740626334135423356791639598109358867770361609232846012626732332450844264293840456574956036349633197336361875) -
               BigInt(5587476413495323413846491678323049250765705078855720721052003556321800113162964567765526724539063327600257543743479921263738432) //
               BigInt(365303089362201664516413596925286161494473575337115296250511752859728108868696929614024803255122785403232359817965288739565550625) *
               6^(1 // 2)
    A[12, 9] = BigInt(54598539818083615233566148602203244896696958910734339754065270985433507945162707737759469214674480807272210648148477499238783276259328) //
               BigInt(301247919092298852634886875129959310794662932014184499827145075851637298698312074030567479239502011693447423026416040794479934024058125) -
               BigInt(6526172450962537747372702280281321524894343532103481802188740153783862532174342615150135214261625966637100811092384548036046488576) //
               BigInt(86490932843037281836028387921320502668579653176624892284566487468170341285762869374265713247057712228954184044334206372230816544375) *
               6^(1 // 2)
    A[12, 10] = BigInt(9391667348404584010955422210328707125006120661611061908889750805619418785820948002455890360939221912190524731087070645107486913457760000000) //
                BigInt(58157266968773020612419028503738708303515285854970725662326801531295387265784849843172223645193277229358434488742203091272981931739152584783) -
                BigInt(8108825145085088104344721048166325225173729495689364696426720161112012414227752328969720658987315654179873760357725235734000399440000000) //
                BigInt(265558296661064021061274102756797754810572081529546692522040189640618206693081506133206500662983001047298787619827411375675716583283801757) *
                6^(1 // 2)
    A[12, 11] = BigInt(123461712659887915177271339396606860810479028777869348014870450606260914019560285661288212498128400476015695960341952) //
                BigInt(281629106670320674754245209358840703704235147307838896741075511220826056829047205614324978253226176275078922716132461)
    A[13, 1] = -BigInt(56042772675322042139227629978042586330633622706053363946766144416933631) //
               BigInt(58808540772323190525590122613223430507352118534557342666015625000000000) +
               BigInt(281404579734699232141455524604487724159024972527) //
               BigInt(1478009944832743180452316204077188415527343750000) * 6^(1 // 2)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = 0
    A[13, 5] = 0
    A[13, 6] = -BigInt(1027163900229750356561238237947225332675621517) //
               BigInt(179261894431132664078747698292867431640625000) -
               BigInt(2745292391641202525373103979336813513372321) //
               BigInt(11702216468464340311060649744558385937500000) * 6^(1 // 2)
    A[13, 7] = -BigInt(157229999853748227305165773364426925282378072238332930121) //
               BigInt(36699907367985458573273204094330716033963413238525390625) +
               BigInt(5757606442802795095318986067317837904184278650664590252101) //
               BigInt(3523191107326604023034227593055748739260487670898437500000) *
               6^(1 // 2)
    A[13, 8] = -BigInt(9311448168593934146015965019904013602133802943325818346622781285907057) //
               BigInt(4255970849010124217193135449668739985401313363005576159362792968750000) -
               BigInt(844213739204097696424366573813463172477074917581) //
               BigInt(4210188359946578336976868164966163024902343750000) * 6^(1 // 2)
    A[13, 9] = BigInt(885774233856672590222951867695327816457340130391639153070521335485617578) //
               BigInt(301098541380295011015469248465465290112505656143757799934635162353515625) -
               BigInt(281404579734699232141455524604487724159024972527) //
               BigInt(284481916364737983221402322504830303192138671875) * 6^(1 // 2)
    A[13, 10] = BigInt(315479116729780153956412124052199685097744239386639023787359107959254802182) //
                BigInt(134481850506505848012587842215515574380212543200894932329128471154748828125) -
                BigInt(2940396453647872276646068776592292229737651937934623) //
                BigInt(7345465058781983710795837429530784777245286520703125) * 6^(1 // 2)
    A[13, 11] = BigInt(2250996163406545378616532039018846586217631599453822541) //
                BigInt(382491303797095993563304148204275636433504028320312500)
    A[13, 12] = BigInt(2689340957307691853294902388334454003959378146957529866233529251986359392336044151708949720958809747970514366293458424272174024493) //
                BigInt(959516386019578808500569114780871708466894752280482835105408027815194895319055443842782227102120493960805649575561796875000000000)
    A[14, 1] = BigInt(47342003848024391498707976847688893013083074441159779465719863625051668939887702630319) //
               BigInt(44802546873926050730401222636656855760802419993852060264615320801485392456054687500000) -
               BigInt(866369530987077991125562402829092187100493209601) //
               BigInt(3325522375873672156017711459173673934936523437500) * 6^(1 // 2)
    A[14, 2] = 0
    A[14, 3] = 0
    A[14, 4] = 0
    A[14, 5] = 0
    A[14, 6] = BigInt(871779321807802447463310035318238762878527157) //
               BigInt(134446420823349498059060773719650573730468750) +
               BigInt(107641268480999396081848975271849857994818) //
               BigInt(1097082793918531904161935913552348681640625) * 6^(1 // 2)
    A[14, 7] = BigInt(496103786351862292800034805114190705484800743513354117014) //
               BigInt(110099722103956375719819612282992148101890239715576171875) -
               BigInt(1329938412606197485769312599390307351191540891599374831099) //
               BigInt(660598332623738254318917673697952888611341438293457031250) *
               6^(1 // 2)
    A[14, 8] = BigInt(40774077277747636354598451708891165494123131383777235229538611989392175193285994266471) //
               BigInt(15264290546248162101058985941588079518256741255377031736357946125713524703979492187500) +
               BigInt(123767075855296855875080343261298883871499029943) //
               BigInt(451091609994276250390378731960660324096679687500) * 6^(1 // 2)
    A[14, 9] = -BigInt(10522038608500556459828649038302068473735749030796372764961618751973793724796364606986664) //
               BigInt(3899417425005422254034574000397382862235892829653375835197340918271556055507659912109375) +
               BigInt(3465478123948311964502249611316368748401972838404) //
               BigInt(2560337247282641848992620902543472728729248046875) * 6^(1 // 2)
    A[14, 10] = -BigInt(27843764471262693189365201135620670490328475323282820219474851621693895769527094334687108984) //
                BigInt(12257041066285164222002594300605593929434139193022166317802121412999357024704596261133984375) +
                BigInt(574774300271998598683873114105472016699241495055292) //
                BigInt(1049352151254569101542262489932969253892183788671875) * 6^(1 // 2)
    A[14, 11] = -BigInt(34241134351848245624232809437676889009431930503529853032576417589898516) //
                BigInt(5613347824358651981100985009024281007603230062439942682713165283203125)
    A[14, 12] = -BigInt(3432044375893932378102368568052286501033850910516999202088532705211633432793920547702800961532438008401883737341854688972639605334600163938610268855705742764072609) //
                BigInt(1143174106341682260971647690410567292143926198650927778920823267461111371275907599801714870165813394147519068210931766844494994616580258435518181434575195312500000)
    A[14, 13] = BigInt(4746930876023919335079451612726717649218264199984) //
                BigInt(18592065538407049755200144388134089346432755594877)
    A[15, 1] = -BigInt(25188329249258825443748527038142409879923012133738985313265430932280250855708601) //
               BigInt(11370641325574469312056961874077298550827642308774647316995717036347558064286250) +
               BigInt(1234273058981860170179592598535508631343082535549881956) //
               BigInt(2105633771469628744518390642968552144069898845895808125) * 6^(1 // 2)
    A[15, 2] = 0
    A[15, 3] = 0
    A[15, 4] = 0
    A[15, 5] = 0
    A[15, 6] = -BigInt(54821142119685055562477216205428613949905430396088) //
               BigInt(3959439837009461289085587746748097947393101278095) -
               BigInt(1511276753825982856072891469504471256664975925000) //
               BigInt(40386286337496505148672995016830599063409633036569) * 6^(1 // 2)
    A[15, 7] = -BigInt(60922424274061599918603524049390657305431262635197540405697952) //
               BigInt(6484861747489032169774584624759953148531564032417461909516875) +
               BigInt(84558575751635978733109961893984238786929550462615375699341616) //
               BigInt(19454585242467096509323753874279859445594692097252385728550625) *
               6^(1 // 2)
    A[15, 8] = -BigInt(116118147575045169733222875835719955334334798191459879782123534889390467935109772) //
               BigInt(8810626901954835245672275131295870892503713957512170681453300814988417642493125) -
               BigInt(176324722711694310025656085505072661620440362221411708) //
               BigInt(285619406719829107485771207042040133465420149964555625) * 6^(1 // 2)
    A[15, 9] = BigInt(17769448722513898342276837490665097286927607247073335618566987143467294900183033216) //
               BigInt(2551217008137889615056342146084561867122485163596619283719957742418751029506356875) -
               BigInt(19748368943709762722873481576568138101489320568798111296) //
               BigInt(6484554262322259071286545935997129135111813687175650625) * 6^(1 // 2)
    A[15, 10] = BigInt(97659266139124074818193264801929547781659926543786381510190954184218570746215033823993530000000) //
                BigInt(18560076654469706205963482908787056850812308205603127326855360961727608242796551101182080033599) -
                BigInt(85297084611782122474911131363078900058888025224607913745000000) //
                BigInt(69210659450201393843166746722954036326338355649915383851733911) *
                6^(1 // 2)
    A[15, 11] = BigInt(473389749049752963256114649231353822492912259509649519870869750525) //
                BigInt(35412440882360341799798842428365422941216508121322622479260846291)
    A[15, 12] = BigInt(33351439245158438248073494056784144097872912773415904536400728387690334563968394114702414108807505158106385116468732853458202899966748488718531545706559142895903144848764637) //
                BigInt(2316611025327287427714802011322252886090793904989900621592365627649097578102163572190502232425490606773312310665593424982745744299371285598588298606088543376742054644818966)
    A[15, 13] = -BigInt(38714992656958413389743252726016897599283911682945255636643554687500000) //
                BigInt(48540494926971587499294589382572212036169135429877901702347521300421767)
    A[15, 14] = BigInt(14800250200940323717124616175641261235119295795768814717803955078125) //
                BigInt(33565577125141877760287380588632421223433194078156948298488471160489)
    A[16, 1] = BigInt(2305785696086397561080858186939897173645641331085041313944389849986584101287) //
               BigInt(617508244345282265819087370078275122671246164669900462139876057008239440000) -
               BigInt(85404623305589712632165905233974183137607899140719) //
               BigInt(124822287169084833758410283469525117460541643292500) * 6^(1 // 2)
    A[16, 2] = 0
    A[16, 3] = 0
    A[16, 4] = 0
    A[16, 5] = 0
    A[16, 6] = BigInt(102903996961580448264190625267026062654799259083) //
               BigInt(5046398084890004857481629999673320438819484730) +
               BigInt(41320925487304219313300272052128374567081128125) //
               BigInt(51473260465878049546312625996667868475958744246) * 6^(1 // 2)
    A[16, 7] = BigInt(62798443349876457506718920843975661399949564598018488144466) //
               BigInt(4132553498782573324058263582553715220777051359780141380625) -
               BigInt(72308807081932961554425711089716771013571419950657300729103) //
               BigInt(12397660496347719972174790747661145662331154079340424141875) *
               6^(1 // 2)
    A[16, 8] = BigInt(1794909142126482564390848522924225553221469019751470544959297614654661293377) //
               BigInt(52596481193994264435601626109752988674679691644275456716633975785978672500) +
               BigInt(12200660472227101804595129319139169019658271305817) //
               BigInt(16931561456559959115207709344056578263397760602500) * 6^(1 // 2)
    A[16, 9] = -BigInt(2775244732780109667342845612394739319115662636371477300455747022423270475907256) //
               BigInt(228417153675584029725018045422706955827996328208181619436454383447149337555625) +
               BigInt(341618493222358850528663620935896732550431596562876) //
               BigInt(96101338378773357469245211954911505447551097205625) * 6^(1 // 2)
    A[16, 10] = -BigInt(27680554659769016623530979176727448251292244310769996015342190819068970556083063125000) //
                BigInt(3299557777429648960576561382256606844677258438797072955341581354051375036522231471437) +
                BigInt(4426552127579895373479670356100179759944766558141730312500) //
                BigInt(3077113738667320707748877199804636746494977000658967987677) *
                6^(1 // 2)
    A[16, 11] = -BigInt(292603171929706291053929402159930330736639136252680853622275) //
                BigInt(15473622826279161150227076887290262443510550964275858143964)
    A[16, 12] = -BigInt(9815717129569106988569302193220999343824932084582093647596086931754666098662594153095258988516305165794739744873539829069617203523509136682216933020431) //
                BigInt(286476991170934153076146641094402171801937250068596542931028678669501762253287693294397689327797388113854588113430063939405071979092547998950955940992)
    A[16, 13] = BigInt(2729491144709837905799148766650782532906050298971406518524169921875) //
                BigInt(2158115888622139473142775812109447802920656149243127309253686951469)
    A[16, 14] = 0
    A[16, 15] = 0
    α[1] = 8198160366203173411119943711500331 // 561057579384085860167277847128765528
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = -BigInt(455655493073428838813281446213740000000) //
           BigInt(1163808011150910561240464225837312497869)
    α[9] = BigInt(19965163648706008081135075746915614720000000) //
           BigInt(86394404190537086868394686205782432516544599)
    α[10] = BigInt(89231107919981418705566970804343750000000000000000000000) //
            BigInt(699979870988335674445594679856445060562597693583175985391)
    α[11] = 47104273954945906713184913871143492 // 209684639122339601934631113492763467
    α[12] = BigInt(20845004421404500464010584740796750650832176798370383084226351294730731196673647311062330972740734737279503119387627146381678677156136042524139311907482802844083) //
            BigInt(36670849891136373020238225328265100250605144718501926305140966586758054847604681466336103169284755987753542321202462371554120593858149755539878561976786592389608)
    α[13] = BigInt(6053037282142306509795911286909179687500000000) //
            BigInt(103899257350518063455290077573775162739725126989)
    α[14] = BigInt(917401104920993498360358406096725463867187500) //
            BigInt(6724249815911346653315790737453607382989551463)
    α[15] = 2585449557665268951371699596493957 // 84574345160764140163208606048427531
    α[16] = 0
    αEEst[1] = 552562031208180939317806684253 // 27669654257734667858523344041464
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 0
    αEEst[7] = 0
    αEEst[8] = 221223388631423597589898601690000000 // 100946136798587090054685074667127461
    αEEst[9] = BigInt(101835408791305297984657812561920000000) //
               BigInt(1149763833200743759976506650241312100139)
    αEEst[10] = BigInt(1313720309077630014453239843750000000000000000000) //
                BigInt(11518201923215510989126466531107437037395719117133)
    αEEst[11] = 4833611232701440504508086151728 // 19081321241454145230196661524503
    αEEst[12] = -BigInt(2129662374582324648106919795703373645353118273066742230724172731025813964712473647144010599206669825382719359113196238857709025512340589957) //
                BigInt(1035543739272367080885190546201097218891268728118207332592595987554851882972292670881794178380097716583123063485287435793657425889233080568)
    αEEst[13] = BigInt(1084761591753640855844358063964843750000000) //
                BigInt(3182895486031249071938549691320502488733423)
    αEEst[14] = 0
    αEEst[15] = 0
    αEEst[16] = 1839190071060649887127895100784 // 38045139523510634351420875415397

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 9, αEEst = αEEst, adaptiveorder = 8,
        stability_size = 4.476172722561913))
end

"""
Journal of Applied Mathematics & Decision Sciences, 4(2), 183-192 (2000),
 "High order explicit Runge-Kutta pairs for ephemerides of the Solar System and the Moon".
"""
function constructSharp9(T::Type = Float64)
    A = zeros(T, 16, 16)
    c = zeros(T, 16)
    α = zeros(T, 16)
    αEEst = zeros(T, 16)

    c[2] = 1 // 50
    c[3] = BigInt(3837236) // BigInt(48429375) +
           BigInt(1031368) // BigInt(145288125) * 6^(1 // 2)
    c[4] = BigInt(1918618) // BigInt(16143125) +
           BigInt(515684) // BigInt(48429375) * 6^(1 // 2)
    c[5] = 14 // 45
    c[6] = BigInt(156) // BigInt(625) + BigInt(26) // BigInt(625) * 6^(1 // 2)
    c[7] = BigInt(156) // BigInt(625) - BigInt(26) // BigInt(625) * 6^(1 // 2)
    c[8] = 52 // 125
    c[9] = 39 // 125
    c[10] = 21 // 200
    c[11] = 280 // 477
    c[12] = 3658227035053715 // 5349704719299032
    c[13] = 247 // 281
    c[14] = 229 // 250
    c[15] = 1
    c[16] = 1
    A[2, 1] = 1 // 50
    A[3, 1] = -BigInt(24000387317036) // BigInt(281448523546875) -
              BigInt(5917264532296) // BigInt(281448523546875) * 6^(1 // 2)
    A[3, 2] = BigInt(46300580261936) // BigInt(281448523546875) +
              BigInt(7915204837696) // BigInt(281448523546875) * 6^(1 // 2)
    A[4, 1] = BigInt(959309) // BigInt(32286250) +
              BigInt(128921) // BigInt(48429375) * 6^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = BigInt(2877927) // BigInt(32286250) +
              BigInt(128921) // BigInt(16143125) * 6^(1 // 2)
    A[5, 1] = BigInt(2826523628723851) // BigInt(5953434698904030) -
              BigInt(68459492317475) // BigInt(595343469890403) * 6^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -BigInt(704240024458145) // BigInt(396895646593602) +
              BigInt(91277530807085) // BigInt(198447823296801) * 6^(1 // 2)
    A[5, 4] = BigInt(958925642225180) // BigInt(595343469890403) -
              BigInt(205373100103780) // BigInt(595343469890403) * 6^(1 // 2)
    A[6, 1] = BigInt(376341108) // BigInt(9406484375) +
              BigInt(207933466) // BigInt(65845390625) * 6^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = BigInt(4343545768844529) // BigInt(27892881885795625) +
              BigInt(469265141246109) // BigInt(27892881885795625) * 6^(1 // 2)
    A[6, 5] = BigInt(1559927818449) // BigInt(28957835234375) +
              BigInt(4382126882523) // BigInt(202704846640625) * 6^(1 // 2)
    A[7, 1] = BigInt(11781705468) // BigInt(235162109375) +
              BigInt(2328587014) // BigInt(1646134765625) * 6^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = BigInt(23459106068523828440829) // BigInt(354298872323611753203125) +
              BigInt(7870375504052283205581) // BigInt(354298872323611753203125) *
              6^(1 // 2)
    A[7, 5] = BigInt(146263465360621089) // BigInt(7558718942052734375) -
              BigInt(1881455818308499953) // BigInt(52911032594369140625) * 6^(1 // 2)
    A[7, 6] = BigInt(9444124356888) // BigInt(82889304453125) -
              BigInt(2459298027368) // BigInt(82889304453125) * 6^(1 // 2)
    A[8, 1] = 52 // 1125
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 0
    A[8, 6] = BigInt(208) // BigInt(1125) - BigInt(13) // BigInt(1125) * 6^(1 // 2)
    A[8, 7] = BigInt(208) // BigInt(1125) + BigInt(13) // BigInt(1125) * 6^(1 // 2)
    A[9, 1] = 741 // 16000
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 0
    A[9, 6] = BigInt(2301) // BigInt(16000) - BigInt(897) // BigInt(32000) * 6^(1 // 2)
    A[9, 7] = BigInt(2301) // BigInt(16000) + BigInt(897) // BigInt(32000) * 6^(1 // 2)
    A[9, 8] = -351 // 16000
    A[10, 1] = 35291978967 // 748709478400
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 0
    A[10, 6] = BigInt(23154511989) // BigInt(149741895680) +
               BigInt(39398793) // BigInt(1772093440) * 6^(1 // 2)
    A[10, 7] = BigInt(23154511989) // BigInt(149741895680) -
               BigInt(39398793) // BigInt(1772093440) * 6^(1 // 2)
    A[10, 8] = -6251205429 // 149741895680
    A[10, 9] = -981041103 // 4679434240
    A[11, 1] = 1601589807329134144752443 // 16639785968494158002257920
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = 0
    A[11, 6] = -BigInt(1736562342312744743536201) // BigInt(1109319064566277200150528) -
               BigInt(360257484908262597335743) // BigInt(511993414415204861607936) *
               6^(1 // 2)
    A[11, 7] = -BigInt(1736562342312744743536201) // BigInt(1109319064566277200150528) +
               BigInt(360257484908262597335743) // BigInt(511993414415204861607936) *
               6^(1 // 2)
    A[11, 8] = 512032742176678555764127 // 369773021522092400050176
    A[11, 9] = 248233526294563631278471 // 103998662303088487514112
    A[11, 10] = -3 // 20
    A[12, 1] = -BigInt(131987017608786696357225423387594635612719389206128606880670434178321331969627889057541436355642743061150672386594396559) //
               BigInt(318753926087995555015147926201612010240447228295789462486798116221476093939683123897279961564118214685494052121351290880)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 0
    A[12, 5] = 0
    A[12, 6] = -BigInt(581038619225160876203856834629458675128926705143465192450716448466169075797359178616021045291080972121429188543592047) //
               BigInt(1011917225676176365127453733973371461080784851732664960275549575306273314094232139356444322425772110112679530543972352) +
               BigInt(71348279807898814965088233729737906656351096659541311591981953238552744805776480028282602780942939988708855996996031) //
               BigInt(51893191060316736673202755588378023645168453935008459501310234631090939184319596890074067816706262057060488745844736) *
               6^(1 // 2)
    A[12, 7] = -BigInt(581038619225160876203856834629458675128926705143465192450716448466169075797359178616021045291080972121429188543592047) //
               BigInt(1011917225676176365127453733973371461080784851732664960275549575306273314094232139356444322425772110112679530543972352) -
               BigInt(71348279807898814965088233729737906656351096659541311591981953238552744805776480028282602780942939988708855996996031) //
               BigInt(51893191060316736673202755588378023645168453935008459501310234631090939184319596890074067816706262057060488745844736) *
               6^(1 // 2)
    A[12, 8] = -BigInt(189357008262607724321683086336517345228379250897103291049044350530935228180690663776657891613652665009511679250229667441) //
               BigInt(104902085728430283184879370421906174798708029629619600881898639306750333561102065113284728091471708748347777999725133824)
    A[12, 9] = -BigInt(1618350992792815653992284152254111827399426534014847245801101845172567304269189800544372100050869595166981551925667441) //
               BigInt(19637518660778297585754649024920739916598981028937029385347383945787366501641192454385997632075140011874187139618963456)
    A[12, 10] = BigInt(6883437842714982754414155283530543027800010156600147069119889350771791431366439329656536871565378282089012991331513) //
                BigInt(1827181489551794784669860898707808352423218653885642080180242252960011545073999200946066370836641132319880849653760)
    A[12, 11] = BigInt(115590271440716912566235566233889746097162479804636463234298604185457969653794053637008425503953091180886565) //
                BigInt(315361333249071836330411464879245754163964468656963767284078726004630975015921697525870414526347596936773632)
    A[13, 1] = BigInt(5215174783558918407997583468635543407988332719241764605769949554629) //
               BigInt(20283132613214812064685094275151111714651171227532533713038580121600)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = 0
    A[13, 5] = 0
    A[13, 6] = BigInt(2843598186227456480865065344408178581293412110128603) //
               BigInt(792075053175002139265844335272820716293355019960320) -
               BigInt(18227070890226867447840942666790512323422585544257) //
               BigInt(121857700488461867579360666965049340968208464609280) * 6^(1 // 2)
    A[13, 7] = BigInt(2843598186227456480865065344408178581293412110128603) //
               BigInt(792075053175002139265844335272820716293355019960320) +
               BigInt(18227070890226867447840942666790512323422585544257) //
               BigInt(121857700488461867579360666965049340968208464609280) * 6^(1 // 2)
    A[13, 8] = BigInt(9326829464422062118248457481351539504275339476759467047326605595685633) //
               BigInt(4901901791858228863857691041029309678010355547895721285919177263022080)
    A[13, 9] = -BigInt(741604155090542466856213236072374206251235617068304762316465738169791) //
               BigInt(141551673163321136844445993892555326037025917405403892742525852712960)
    A[13, 10] = -BigInt(6058504866441219655595548618762485399974773685307046001179355536003) //
                BigInt(2252275720815396172726400694965157641073696835574259179818290290400)
    A[13, 11] = -BigInt(72917047186465183128180555150230405657138451692847535142343993) //
                BigInt(44661747288016218276854771442831738093234145203222656783563600)
    A[13, 12] = BigInt(2736153920540927643774133147635296486946660915558253285983742020488887296849241173151960647763453239551016003889152) //
                BigInt(2485672110698341015290264470463939203955869249618375406787169018009688457749866177826801192710345262847046284166825)
    A[14, 1] = BigInt(1961431625890315687063141575818232405522545898155499982338718373117379429883) //
               BigInt(480056647167077429990593568055406093586176318669944422673481728000000000000)
    A[14, 2] = 0
    A[14, 3] = 0
    A[14, 4] = 0
    A[14, 5] = 0
    A[14, 6] = -BigInt(8688525606146315530022414580346392155721271039) //
               BigInt(22386738118754433181814607602481176248320000000) -
               BigInt(10256190098435854298655077997613296122112139953) //
               BigInt(1148037852243817086246902953973393653760000000) * 6^(1 // 2)
    A[14, 7] = -BigInt(8688525606146315530022414580346392155721271039) //
               BigInt(22386738118754433181814607602481176248320000000) +
               BigInt(10256190098435854298655077997613296122112139953) //
               BigInt(1148037852243817086246902953973393653760000000) * 6^(1 // 2)
    A[14, 8] = -BigInt(108151392092290424953498836380059772609736403739434481043071361807712075869481) //
               BigInt(8600735495194563448316261478331353993230137106756553600705167485829120000000)
    A[14, 9] = BigInt(683210554257935462600257975958139742203919396113084127371502375524416129719) //
               BigInt(26895337200565243662247103690698994332502640106760065066162305761280000000)
    A[14, 10] = -BigInt(125971034051203704183074450363446847441594334546885083244594242327104115033) //
                BigInt(5066049934698363488698655054901069679084758735799331062593807151200000000)
    A[14, 11] = BigInt(4322338495495152743252505005837177994220267688026960252214552638944423) //
                BigInt(236867625787508422152958167179676757535142999000924357630500000000000)
    A[14, 12] = -BigInt(88682414394183619425441647866243388112917289239161463940944492930492112547171652363240146123589908870567811533658125375935101390832) //
                BigInt(9405104776230176067202383689444499684238639823097650089954333383377958770491999076242219980453370248022825420814384818872314453125)
    A[14, 13] = BigInt(26235475641986625187247554297838197168935151270802587) //
                BigInt(31781620957198174033817415268740604591106877500000000)
    A[15, 1] = -BigInt(2933688768685553737193922190442902414638569907165819426999847151894747) //
               BigInt(1423967854813137802350516795065201258930107696470226170813903745843200)
    A[15, 2] = 0
    A[15, 3] = 0
    A[15, 4] = 0
    A[15, 5] = 0
    A[15, 6] = -BigInt(279050827135618188106138704976571118076242172562777) //
               BigInt(26980717750745660055932121988692169249262917386240) +
               BigInt(59017804198407615229179283246229064921710388893173) //
               BigInt(17987145167163773370621414659128112832841944924160) * 6^(1 // 2)
    A[15, 7] = -BigInt(279050827135618188106138704976571118076242172562777) //
               BigInt(26980717750745660055932121988692169249262917386240) -
               BigInt(59017804198407615229179283246229064921710388893173) //
               BigInt(17987145167163773370621414659128112832841944924160) * 6^(1 // 2)
    A[15, 8] = BigInt(68240477823918559060550996013166770535743446467404475965020846786328901) //
               BigInt(69628625900822775316247857716283393138506964199995817446989013471723520)
    A[15, 9] = BigInt(48531604865335743440838806675493568975092395234916265724406574203650554879) //
               BigInt(7529075569049450724715447951105730327391304314424323559158913160761835520)
    A[15, 10] = BigInt(7315898198049114373691779027237206235234893868747090317226910860963581499) //
                BigInt(432400219379131684183655517400956805867024783959757582533074834771793600)
    A[15, 11] = -BigInt(87035912584683752124645187592152267644073875904388006117245587111831) //
                BigInt(41468532532723053663401983927439573730970639521941633396843682248800)
    A[15, 12] = BigInt(793006054328041651061360131256412474400253089909554005378332728214806108995212138291759017448087224471716436232175864384424753159293287828190208) //
                BigInt(600635305507048430170531323174007915813446499853348186301444535722932552444071789775577508083955280258070844780520371987041149053341617704219325)
    A[15, 13] = BigInt(19013238692887784267164981427867630356262081870600946422701364458) //
                BigInt(146516308633144198110735805762400905606463733191840985648075179899)
    A[15, 14] = BigInt(368176545506575596342007241113258886329861009608750000000) //
                BigInt(7515329389098801941975451526298754679007062667248055263091)
    A[16, 1] = -BigInt(36388658330162124762200023703074655379362961851837455245313588466117) //
               BigInt(299291302137383314536268908335335078201218487321808786511993451315200)
    A[16, 2] = 0
    A[16, 3] = 0
    A[16, 4] = 0
    A[16, 5] = 0
    A[16, 6] = -BigInt(5017294099975580158862668031284197043053591531405121) //
               BigInt(432874134925719951608236287739966681862399087083520) -
               BigInt(115700422823857939498444446575144266776173664871303) //
               BigInt(199788062273409208434570594341523083936491886346240) * 6^(1 // 2)
    A[16, 7] = -BigInt(5017294099975580158862668031284197043053591531405121) //
               BigInt(432874134925719951608236287739966681862399087083520) +
               BigInt(115700422823857939498444446575144266776173664871303) //
               BigInt(199788062273409208434570594341523083936491886346240) * 6^(1 // 2)
    A[16, 8] = -BigInt(91869384706617020415871523809581333688744319256441669606452442156503951243) //
               BigInt(10053990532496870785611049673131231012151168225596329651077879646797168640)
    A[16, 9] = BigInt(52862999381403119807509472978743982056878734540171187101666495333163485251) //
               BigInt(2399905281474961058988767393522446655063434938667458594324573249770618880)
    A[16, 10] = BigInt(11964965861294434337427534231330501089458731146841410298258149571218167) //
                BigInt(1974017923973741961134058915802491373563600428420464580517659076868800)
    A[16, 11] = BigInt(37167680872257703003686692191635149388479305578570534942584948859) //
                BigInt(3913390524272247558198063815920084715858046504866006220243936800)
    A[16, 12] = -BigInt(29443955867054347753341026121045589578978185460220369427665634428174791788280117223531690217195051250854008448256016995751289856) //
                BigInt(6138364406531943832091734698181735919618182422997726728173639255839442923939727454711783369110119257425239678616468261908910975)
    A[16, 13] = BigInt(407816748385172686498153181346812432791118177175769818363629863) //
                BigInt(626162397882386095196201629768303759628250075752997120162580989)
    A[16, 14] = 0
    A[16, 15] = 0
    α[1] = 30703843389361946002220520407 // 1036329015084155723633962896000
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = 1516681888913470906364013671875 // 19423768214582439936604117641536
    α[9] = BigInt(1929922737998470573359614532470703125) //
           BigInt(9295447834009061726737853188569292704)
    α[10] = BigInt(27072397368129209968072433152000000000) //
            BigInt(159540891067276798629433718421290211669)
    α[11] = BigInt(3416676287738448149119878197304164096817920457) //
            BigInt(22521752441211566270536786917243920830369456000)
    α[12] = BigInt(909034900749411645631439991260524977916886591502548355130330148829066896764151555292038222333366816993556860935646735988456500531298304) //
            BigInt(6301978749188979317659380355882211371188146506066543226107226217493443986031316306450151922600620534579104501042337690306078523205079625)
    α[13] = BigInt(9160897746149204383653282352747804858423571) //
            BigInt(54934119002888850773584011583391921191449440)
    α[14] = 3769686146953412690297035156250000 // 195792979665408643382362918863397227
    α[15] = 50782110772148063247179059 // 1538266148871578545201811280
    α[16] = 0
    αEEst[1] = 135131455470598097879473933 // 4525454214341291369580624000
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 0
    αEEst[7] = 0
    αEEst[8] = 2518169234679274570156341552734375 // 38284247150941989115046715871467456
    αEEst[9] = 13171020424136540706261627197265625 // 61559257178867958455217570785227104
    αEEst[10] = 33191111003144264098986272000000000 // 196721197370254992144801132455351679
    αEEst[11] = BigInt(98603841096694858013088556726735239713679) //
                BigInt(574051243626833692823306576536179258032000)
    αEEst[12] = BigInt(5093635768538576107415789300891145054334065554846760039946162068969260804962050605021707042060987973190001839145435136) //
                BigInt(40589106364182299226269510371318864016087966429387805447847648739509460272079710913035195393529050583646417571607929125)
    αEEst[13] = BigInt(108010721096523379193662759959856611609133) //
                BigInt(570689193181223151553200582051397411377120)
    αEEst[14] = 0
    αEEst[15] = 0
    αEEst[16] = 26859551018855966185191031 // 763900876650511794556001520

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 9, αEEst = αEEst, adaptiveorder = 8,
        stability_size = 5.191689996466664))
end

"""
Optimized explicit Runge-Kutta pairs of order 9(8), by Ch. Tsitouras,
 Applied Numerical Mathematics, 38 (2001) 123-134.
"""
function constructTsitouras9(T::Type = Float64)
    A = zeros(T, 16, 16)
    c = zeros(T, 16)
    α = zeros(T, 16)
    αEEst = zeros(T, 16)

    c[2] = 1 // 49
    c[3] = BigInt(64) // BigInt(705) - BigInt(16) // BigInt(14805) * 6^(1 // 2)
    c[4] = BigInt(32) // BigInt(235) - BigInt(8) // BigInt(4935) * 6^(1 // 2)
    c[5] = 3 // 7
    c[6] = BigInt(8) // BigInt(21) + BigInt(4) // BigInt(63) * 6^(1 // 2)
    c[7] = BigInt(8) // BigInt(21) - BigInt(4) // BigInt(63) * 6^(1 // 2)
    c[8] = 40 // 63
    c[9] = 10 // 21
    c[10] = 19 // 18
    c[11] = 7 // 9
    c[12] = 319999786 // 2170712113
    c[13] = 15 // 16
    c[14] = 39 // 40
    c[15] = 1
    c[16] = 1
    A[2, 1] = 1 // 49
    A[3, 1] = -BigInt(165952) // BigInt(1491075) +
              BigInt(38896) // BigInt(10437525) * 6^(1 // 2)
    A[3, 2] = BigInt(301312) // BigInt(1491075) -
              BigInt(7168) // BigInt(1491075) * 6^(1 // 2)
    A[4, 1] = BigInt(8) // BigInt(235) - BigInt(2) // BigInt(4935) * 6^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = BigInt(24) // BigInt(235) - BigInt(2) // BigInt(1645) * 6^(1 // 2)
    A[5, 1] = BigInt(38937) // BigInt(44800) + BigInt(171) // BigInt(5600) * 6^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -BigInt(149931) // BigInt(44800) - BigInt(81) // BigInt(700) * 6^(1 // 2)
    A[5, 4] = BigInt(65097) // BigInt(22400) + BigInt(477) // BigInt(5600) * 6^(1 // 2)
    A[6, 1] = BigInt(176) // BigInt(5103) - BigInt(29) // BigInt(5103) * 6^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = BigInt(364520) // BigInt(1674351) +
              BigInt(87715) // BigInt(5023053) * 6^(1 // 2)
    A[6, 5] = BigInt(1940224) // BigInt(15069159) +
              BigInt(779264) // BigInt(15069159) * 6^(1 // 2)
    A[7, 1] = BigInt(4336) // BigInt(127575) + BigInt(479) // BigInt(127575) * 6^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = BigInt(90731944) // BigInt(400648275) -
              BigInt(170142739) // BigInt(8413613775) * 6^(1 // 2)
    A[7, 5] = BigInt(8245504) // BigInt(62429373) -
              BigInt(22187008) // BigInt(437005611) * 6^(1 // 2)
    A[7, 6] = -BigInt(3936) // BigInt(340025) +
              BigInt(11464) // BigInt(3060225) * 6^(1 // 2)
    A[8, 1] = 40 // 567
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 0
    A[8, 6] = BigInt(160) // BigInt(567) - BigInt(10) // BigInt(567) * 6^(1 // 2)
    A[8, 7] = 10 // 567 * 6^(1 // 2) + 160 // 567
    A[9, 1] = 95 // 1344
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 0
    A[9, 6] = BigInt(295) // BigInt(1344) - BigInt(115) // BigInt(2688) * 6^(1 // 2)
    A[9, 7] = BigInt(295) // BigInt(1344) + BigInt(115) // BigInt(2688) * 6^(1 // 2)
    A[9, 8] = -15 // 448
    A[10, 1] = 52918819 // 138240000
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 0
    A[10, 6] = -BigInt(1453047743) // BigInt(103680000) -
               BigInt(4153586941) // BigInt(829440000) * 6^(1 // 2)
    A[10, 7] = -BigInt(1453047743) // BigInt(103680000) +
               BigInt(4153586941) // BigInt(829440000) * 6^(1 // 2)
    A[10, 8] = 44599023 // 5120000
    A[10, 9] = 518179039 // 25920000
    A[11, 1] = BigInt(258780283) // BigInt(8618400000) +
               BigInt(585428803) // BigInt(51710400000) * 6^(1 // 2)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = 0
    A[11, 6] = 19 // 25
    A[11, 7] = BigInt(1180508473123) // BigInt(443296800000) -
               BigInt(136404911099) // BigInt(147765600000) * 6^(1 // 2)
    A[11, 8] = -BigInt(106856621) // BigInt(190800000) +
               BigInt(585428803) // BigInt(2289600000) * 6^(1 // 2)
    A[11, 9] = -BigInt(1260561943) // BigInt(591300000) +
               BigInt(585428803) // BigInt(886950000) * 6^(1 // 2)
    A[11, 10] = BigInt(13167297224) // BigInt(792049782825) -
                BigInt(9366860848) // BigInt(2376149348475) * 6^(1 // 2)
    A[12, 1] = BigInt(307213395328582867964430765847473084972824867512957518186088963) //
               BigInt(5126364212860621132939944111710304798478633358572140981841000000) +
               BigInt(119107533326819222510639750832411974467191643469020133053) //
               BigInt(29137664905764716334007503363213354847664112851105748578125) *
               6^(1 // 2)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 0
    A[12, 5] = 0
    A[12, 6] = BigInt(10354821182100230493026667000379184955505622050895245676387169) //
               BigInt(146467548938874889512569831763151565670818095959204028052600000) -
               BigInt(3268463788087907168885902319404754223684763360609854606186699) //
               BigInt(41847871125392825575019951932329018763090884559772579443600000) *
               6^(1 // 2)
    A[12, 7] = BigInt(542371157260956891298011197530777699174343335570401188235055068891) //
               BigInt(715703215922030799396778727922657043395761853183510539934169000000) -
               BigInt(2557103919967567420571445567798380600417779349808054092568810687111) //
               BigInt(10019845022908431191554902190917198607540665944569147559078366000000) *
               6^(1 // 2)
    A[12, 8] = -BigInt(892225578009519154676238995901578841509244882985862325637048827) //
               BigInt(38813900468801845720831005417235164902766795429189067433939000000) +
               BigInt(20367388198886087049319397392342447633889771033202442752063) //
               BigInt(220613748572218566528913954035758258132313997301229239234375) *
               6^(1 // 2)
    A[12, 9] = -BigInt(550688605235770338034863642917595195825050633798073448796494914) //
               BigInt(835322740042020854251374821774223772966384453517335472487484375) +
               BigInt(72417380262706087286468968506106480476052519229164240896224) //
               BigInt(303864219731546327483221106502082129125640034018674235171875) *
               6^(1 // 2)
    A[12, 10] = BigInt(123396895115495738434549229715587040998178289648593721222723693824) //
                BigInt(16960466422214315122555088502288155618001490986047888958907461648625) -
                BigInt(8781560217119727637264447549372070053517105489683916369731584) //
                BigInt(6169685857480653009296139869875647732994358307038155314262445125) *
                6^(1 // 2)
    A[12, 11] = -BigInt(241890129426298647138485610377551406165672225318019246672) //
                BigInt(3661688723471872237814245794078789141770452398980100701315)
    A[13, 1] = BigInt(45077846760256141387004276823) //
               BigInt(110315894143992133591739924480) +
               BigInt(1493491403898138129099) // BigInt(13100021190238236835840) *
               6^(1 // 2)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = 0
    A[13, 5] = 0
    A[13, 6] = BigInt(28530732123103900185) // BigInt(9849639992660328448) -
               BigInt(15062887306567756845) // BigInt(5628365710091616256) * 6^(1 // 2)
    A[13, 7] = BigInt(530875502237315716994493) // BigInt(24064781139210466754560) -
               BigInt(8920823473649531766699837) // BigInt(1347627743795786138255360) *
               6^(1 // 2)
    A[13, 8] = -BigInt(155850251753928802974915857362119) //
               BigInt(174015086605340016477040216637440) +
               BigInt(13441422635083243161891) // BigInt(5220309196109974077440) *
               6^(1 // 2)
    A[13, 9] = -BigInt(493074073683718697930133408597) //
               BigInt(27602712116408194083051274240) +
               BigInt(1493491403898138129099) // BigInt(224694912332563742720) * 6^(1 // 2)
    A[13, 10] = BigInt(200609996314078300148532240828075) //
                BigInt(1019933691979646265167106381709312) -
                BigInt(336035565877081079047275) // BigInt(8465066424794973677551616) *
                6^(1 // 2)
    A[13, 11] = -1259978731825102407292471875 // 947642075600343143202947072
    A[13, 12] = -BigInt(193916214235317468987992391599053188049367133486207120889311375) //
                BigInt(42748092349455088111344007455417233641020280098816132254793728)
    A[14, 1] = BigInt(36716621212098036093935018687105425505961) //
               BigInt(72248275402215258274603114496000000000000) +
               1221461237263884679751555607 // 9994523002806272000000000000 * 6^(1 // 2)
    A[14, 2] = 0
    A[14, 3] = 0
    A[14, 4] = 0
    A[14, 5] = 0
    A[14, 6] = BigInt(8283471074731862302286097) // BigInt(7514678949478400000000000) -
               BigInt(15668946773152185221466849) // BigInt(4294102256844800000000000) *
               6^(1 // 2)
    A[14, 7] = BigInt(397408075485926915758262202639) //
               BigInt(18359970961922048000000000000) -
               BigInt(6493922587539771225254133441201) //
               BigInt(1028158373867634688000000000000) * 6^(1 // 2)
    A[14, 8] = BigInt(2435278493903047909370803905780425755549361787) //
               BigInt(6951941874960961660657094977773568000000000000) +
               10993151135374962117764000463 // 3982779843223552000000000000 * 6^(1 // 2)
    A[14, 9] = -BigInt(558123239069103416347126929975086912148938889) //
               BigInt(34184771406932232290260692794368000000000000) +
               1221461237263884679751555607 // 171428613534976000000000000 * 6^(1 // 2)
    A[14, 10] = BigInt(306734586161727173704146823378382330889519) //
                BigInt(1406266800626214225846914661737584000000000) -
                10993151135374962117764000463 // 258333325951995046312000000000 * 6^(1 // 2)
    A[14, 11] = -BigInt(163845778835264660255510638493965671483) //
                BigInt(114196173990354810149157741209600000000)
    A[14, 12] = -BigInt(379727098691580451304129337662817719784451899678250630685021894233169654109523405323) //
                BigInt(74980943726337976053062781324716052947047103330285761853079535691020475760640000000)
    A[14, 13] = -235412270220829707518634576 // 10004921377982463725322265625
    A[15, 1] = BigInt(246936976626965995144662615055843) //
               BigInt(458607216908592948032148676000000) +
               BigInt(1082690484492446489) // BigInt(9760276369928000000) * 6^(1 // 2)
    A[15, 2] = 0
    A[15, 3] = 0
    A[15, 4] = 0
    A[15, 5] = 0
    A[15, 6] = -BigInt(5604689614035063) // BigInt(3669276830800000) -
               BigInt(556441805223969) // BigInt(131045601100000) * 6^(1 // 2)
    A[15, 7] = BigInt(613864878280672501731) // BigInt(35859318285004000000) -
               BigInt(1204548796910252147313) // BigInt(251015227995028000000) * 6^(1 // 2)
    A[15, 8] = BigInt(2797434125723089631989036145593862912367) //
               BigInt(1454544479209361363081818944324206000000) +
               9744214360432018401 // 3889433440648000000 * 6^(1 // 2)
    A[15, 9] = -BigInt(78495184115962971170097375833699424291) //
               BigInt(6993857040773733071103139785566500000) +
               1082690484492446489 // 167410755405250000 * 6^(1 // 2)
    A[15, 10] = BigInt(38933993403647216169460437449231039304) //
                BigInt(187659712667940353060087104196403689875) -
                155907429766912294416 // 4036458217999922598625 * 6^(1 // 2)
    A[15, 11] = -4116079644901049270506242987 // 3129604261871291193409631155
    A[15, 12] = -BigInt(869293931367062297302433741250572190497190440843976042585249416347005728190566119731) //
                BigInt(185354770907948507678541561399002564941703939116570719233261621774655604595976403844)
    A[15, 13] = 37424319425692041216 // 5567956245138066768875
    A[15, 14] = -26637096887808000000 // 690666426738105277187
    A[16, 1] = BigInt(29766990313562078086727231879295211) //
               BigInt(61911974282660047984340071260000000) +
               BigInt(13591007763158148317) // BigInt(146404145548920000000) * 6^(1 // 2)
    A[16, 2] = 0
    A[16, 3] = 0
    A[16, 4] = 0
    A[16, 5] = 0
    A[16, 6] = -BigInt(583353427111301443) // BigInt(293542146464000000) -
               BigInt(40076566638916409) // BigInt(10483648088000000) * 6^(1 // 2)
    A[16, 7] = BigInt(19531933359906438277223) // BigInt(1434372731400160000000) -
               BigInt(9404496119706140958251) // BigInt(2510152279950280000000) * 6^(1 // 2)
    A[16, 8] = BigInt(117813280551292383465368402104818510183261) //
               BigInt(58181779168374454523272757772968240000000) +
               40773023289474444951 // 19447167203240000000 * 6^(1 // 2)
    A[16, 9] = -BigInt(2325450203082696961260142930045971433803) //
               BigInt(279754281630949322844125591422660000000) +
               13591007763158148317 // 2511161331078750000 * 6^(1 // 2)
    A[16, 10] = BigInt(2829276993010365332321333380796499415893) //
                BigInt(15951075576774930010107403856694313639375) -
                652368372631591119216 // 20182291089999612993125 * 6^(1 // 2)
    A[16, 11] = -64853022020814701128056854079867 // 63343190260274933754610934577200
    A[16, 12] = -BigInt(5501801061234219240967503472431237166019376078542725490351730718419538863636177307419649) //
                BigInt(1398104557705668743632427777409619346988852569336419139359459089957402274666793446137600)
    A[16, 13] = -6425519161642982676103168 // 293904570399612854395066875
    A[16, 14] = -51926882941360640000 // 2663999074561263212007
    A[16, 15] = 0
    α[1] = 173734691637390647 // 4182794002754640000
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = -72263163141715044860361 // 169939769455665013040000
    α[9] = 14586697891849999254003 // 29700462390576849520000
    α[10] = 102209317997264953344 // 225042304099487188475
    α[11] = 1883570537693211021 // 1872275755054959100
    α[12] = BigInt(17109990417889849939560223376925306674323804078983341334325755071278367152457480027) //
            BigInt(71381427125808828146076534703173195056090466069876919551824450642795947233961810560)
    α[13] = -10678264099993989152768 // 2396652442219114419375
    α[14] = 1212545712242913280000000 // 130535954653501897388343
    α[15] = -5462519910419447 // 852178998090420
    α[16] = 10 // 13
    αEEst[1] = 310663590276876647 // 7768046005115760000
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 0
    αEEst[7] = 0
    αEEst[8] = 3852995276698729629111 // 96052913170593268240000
    αEEst[9] = 444432599960137218801 // 1132275692309381360000
    αEEst[10] = -439332801641475586992 // 417935707613333350025
    αEEst[11] = -217759597939752597 // 869270886275516725
    αEEst[12] = BigInt(17659727425367474410466504999969566436559834034377723810849746362072595341) //
                BigInt(71629450771664553932203312665447432781082741824305141429680586700646203520)
    αEEst[13] = 982070268140959301632 // 143578256907135886875
    αEEst[14] = -3885243675284052992000000 // 242423915785074952292637
    αEEst[15] = 10
    αEEst[16] = 10 // 13

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 9, αEEst = αEEst, adaptiveorder = 8,
        stability_size = 3.939793900092272))
end

"""
Optimized explicit Runge-Kutta pairs of order 9(8), by Ch. Tsitouras,
 Applied Numerical Mathematics, 38 (2001) 123-134.
"""
function constructTsitouras92(T::Type = Float64)
    A = zeros(T, 16, 16)
    c = zeros(T, 16)
    α = zeros(T, 16)
    αEEst = zeros(T, 16)

    c[2] = 1 // 46
    c[3] = BigInt(96755252944) // BigInt(718444993695) -
           BigInt(11256225944) // BigInt(718444993695) * 6^(1 // 2)
    c[4] = BigInt(48377626472) // BigInt(239481664565) -
           BigInt(5628112972) // BigInt(239481664565) * 6^(1 // 2)
    c[5] = 71 // 136
    c[6] = BigInt(276) // BigInt(715) - BigInt(46) // BigInt(715) * 6^(1 // 2)
    c[7] = BigInt(276) // BigInt(715) + BigInt(46) // BigInt(715) * 6^(1 // 2)
    c[8] = 92 // 143
    c[9] = 69 // 143
    c[10] = 3 // 44
    c[11] = 103 // 411
    c[12] = 30258248819701 // 45339732981913
    c[13] = 59 // 69
    c[14] = 44 // 49
    c[15] = 1
    c[16] = 1
    A[2, 1] = 1 // 46
    A[3, 1] = -BigInt(163287951175938724532816) // BigInt(516163208965408589753025) +
              BigInt(42011574289334042817176) // BigInt(516163208965408589753025) *
              6^(1 // 2)
    A[3, 2] = BigInt(232801278267248934720896) // BigInt(516163208965408589753025) -
              BigInt(50098553466700618240256) // BigInt(516163208965408589753025) *
              6^(1 // 2)
    A[4, 1] = BigInt(12094406618) // BigInt(239481664565) -
              BigInt(1407028243) // BigInt(239481664565) * 6^(1 // 2)
    A[4, 2] = 0
    A[4, 3] = BigInt(36283219854) // BigInt(239481664565) -
              BigInt(4221084729) // BigInt(239481664565) * 6^(1 // 2)
    A[5, 1] = BigInt(450479172821804238979159483) // BigInt(489985471732935255816699904) +
              BigInt(65404175703680378526395577) // BigInt(244992735866467627908349952) *
              6^(1 // 2)
    A[5, 2] = 0
    A[5, 3] = -BigInt(1663285823745576633021875313) // BigInt(489985471732935255816699904) -
              BigInt(258991054585998425691922779) // BigInt(244992735866467627908349952) *
              6^(1 // 2)
    A[5, 4] = BigInt(734303944921586208649981787) // BigInt(244992735866467627908349952) +
              BigInt(96793439441159023582763601) // BigInt(122496367933233813954174976) *
              6^(1 // 2)
    A[6, 1] = BigInt(188634486760257) // BigInt(2753187875656075) -
              BigInt(40451003556679) // BigInt(5506375751312150) * 6^(1 // 2)
    A[6, 2] = 0
    A[6, 3] = 0
    A[6, 4] = BigInt(890541395040155939974909749) // BigInt(3404930508779360011084250045) -
              BigInt(235414842445143790083329443) // BigInt(6809861017558720022168500090) *
              6^(1 // 2)
    A[6, 5] = BigInt(127509164130554343284736) // BigInt(2278805333809176804299525) -
              BigInt(51090254569210884816896) // BigInt(2278805333809176804299525) *
              6^(1 // 2)
    A[7, 1] = BigInt(523150756520001) // BigInt(5294592068569375) +
              BigInt(372205675002861) // BigInt(137659393782803750) * 6^(1 // 2)
    A[7, 2] = 0
    A[7, 3] = 0
    A[7, 4] = BigInt(121832502441158811994748302664452173) //
              BigInt(6319431229672072722127362725145820625) -
              BigInt(12054008141355156662680357922224203047) //
              BigInt(164305211971473890775311430853791336250) * 6^(1 // 2)
    A[7, 5] = -BigInt(7345188891123909155979140554752) //
              BigInt(52428978281511938535235507146875) +
              BigInt(71382195182457889488943971467264) //
              BigInt(681576717659655200958061592909375) * 6^(1 // 2)
    A[7, 6] = BigInt(84211752143498940768) // BigInt(206389046233053165625) +
              BigInt(567839841668979868) // BigInt(18762640566641196875) * 6^(1 // 2)
    A[8, 1] = 92 // 1287
    A[8, 2] = 0
    A[8, 3] = 0
    A[8, 4] = 0
    A[8, 5] = 0
    A[8, 6] = BigInt(368) // BigInt(1287) + BigInt(23) // BigInt(1287) * 6^(1 // 2)
    A[8, 7] = BigInt(368) // BigInt(1287) - BigInt(23) // BigInt(1287) * 6^(1 // 2)
    A[9, 1] = 1311 // 18304
    A[9, 2] = 0
    A[9, 3] = 0
    A[9, 4] = 0
    A[9, 5] = 0
    A[9, 6] = BigInt(4071) // BigInt(18304) + BigInt(1587) // BigInt(36608) * 6^(1 // 2)
    A[9, 7] = BigInt(4071) // BigInt(18304) - BigInt(1587) // BigInt(36608) * 6^(1 // 2)
    A[9, 8] = -621 // 18304
    A[10, 1] = 2451872601 // 50434064384
    A[10, 2] = 0
    A[10, 3] = 0
    A[10, 4] = 0
    A[10, 5] = 0
    A[10, 6] = BigInt(84329349) // BigInt(1146228736) -
               BigInt(1383050643) // BigInt(100868128768) * 6^(1 // 2)
    A[10, 7] = BigInt(84329349) // BigInt(1146228736) +
               BigInt(1383050643) // BigInt(100868128768) * 6^(1 // 2)
    A[10, 8] = -1098320769 // 50434064384
    A[10, 9] = -333490521 // 3152129024
    A[11, 1] = -BigInt(11290810941252792923651) // BigInt(1669469461414577748900000) -
               BigInt(76218489460616423924209) // BigInt(10016816768487466493400000) *
               6^(1 // 2)
    A[11, 2] = 0
    A[11, 3] = 0
    A[11, 4] = 0
    A[11, 5] = 0
    A[11, 6] = 1 // 30
    A[11, 7] = -BigInt(44608220078798131601386867) // BigInt(1778327431680661626219300000) -
               BigInt(302663621648107819403033939) // BigInt(5334982295041984878657900000) *
               6^(1 // 2)
    A[11, 8] = BigInt(4768550623191902657077) // BigInt(335320789258483564950000) +
               BigInt(76218489460616423924209) // BigInt(9388982099237539818600000) *
               6^(1 // 2)
    A[11, 9] = BigInt(76371166597983496297729) // BigInt(1268154687036073482337500) +
               BigInt(76218489460616423924209) // BigInt(1902232030554110223506250) *
               6^(1 // 2)
    A[11, 10] = BigInt(12837092726068800321242176) // BigInt(73489499260117750229428125) +
                BigInt(224387232972054752032871296) //
                BigInt(13889515360162254793361915625) * 6^(1 // 2)
    A[12, 1] = -BigInt(355843792738780589211336013211266011894384892859307606673616682840321037292934264483298005378753134077) //
               BigInt(6114980501338466999761271314291724602877525435623381536957589416556451501095664913528943447961301150000) +
               BigInt(1634253006867467370609841039730482111683170725738988407797492220074165441589604495769296514008733014761) //
               BigInt(62508689569237662664226328990537629273859148897483455711122025147021504233422352449406977468048856200000) *
               6^(1 // 2)
    A[12, 2] = 0
    A[12, 3] = 0
    A[12, 4] = 0
    A[12, 5] = 0
    A[12, 6] = -BigInt(17594918026040863488775323023688367263938903884097169692619657825340943285779242173824532021369934883) //
               BigInt(34620197299885474706340736056297763905521990158606221624621429312196525421587764433517710597688597280) -
               BigInt(1718800244468971293911190106580576895962584281439946985382458420029691073118839588510586218092229) //
               BigInt(134447368154895047403264994393389374390376660810121249027656036163870001637234036635020235330829504) *
               6^(1 // 2)
    A[12, 7] = -BigInt(138207307502597872414466248711907117101238289950946956575872497287994333538943652935288845283125425999567221) //
               BigInt(449445892634106961164222874370868126852218761612180707797390011848931329871568905581758203549930436500950000) +
               BigInt(373422609934822782117113512174693455251869921393106586844571126171756173532479385325874622983683540733325899) //
               BigInt(1797783570536427844656891497483472507408875046448722831189560047395725319486275622327032814199721746003800000) *
               6^(1 // 2)
    A[12, 8] = BigInt(12632445836279619113294210780168171565726892406599945732333835965116776912265868042131556830549327654393) //
               BigInt(112996477298237313277639902405971868302745384545450862247028276227308103806571175581620305423011393900000) -
               BigInt(6303547312203088429495101153246145287920801370707526715790327134571780988988474483681572268319398771221) //
               BigInt(225992954596474626555279804811943736605490769090901724494056552454616207613142351163240610846022787800000) *
               6^(1 // 2)
    A[12, 9] = BigInt(1115379785277627713874987636432084522404410826010464628546537384975239192450316138834225701046102632948257) //
               BigInt(3347520639960280405850082108776916598468832017158981794068212683234002575269671076605501548156712544287500) -
               BigInt(1634253006867467370609841039730482111683170725738988407797492220074165441589604495769296514008733014761) //
               BigInt(11870640567235036900177596130414597866910751833897098560525576890900718352020110200728728894172739518750) *
               6^(1 // 2)
    A[12, 10] = BigInt(11195395793619792818419572687595618271094922253700867106612850891912759357989068567623193453186950145251776) //
                BigInt(37146762745465379216384266737567863919875050348378042013698628944987922182903823702648363618722076517534375) -
                BigInt(687320121745403419867910288709505619542179230939368838936545299414049008577113662220686991034529999350912) //
                BigInt(12382254248488459738794755579189287973291683449459347337899542981662640727634607900882787872907358839178125) *
                6^(1 // 2)
    A[12, 11] = BigInt(73825104187474768875967421005730375156586805528272061538441201406502469378399602291332635087412) //
                BigInt(92871272539029599143289144740668914525146202470165909630170359355065425993729905138177287271317)
    A[13, 1] = -BigInt(64719070744144335733962214412431202035561004419937320357188840883) //
               BigInt(101209800434111984325751891003355158288493635359529771992774000000) +
               BigInt(3559372342256314491923633965576963561357143961289) //
               BigInt(36357243706101884581258626118802982597639184500000) * 6^(1 // 2)
    A[13, 2] = 0
    A[13, 3] = 0
    A[13, 4] = 0
    A[13, 5] = 0
    A[13, 6] = -BigInt(66005588850553492194962577164587068367312330247) //
               BigInt(19603612723583140697294028442079849957408380800) -
               BigInt(13176568020345600854085799476910283462657777) //
               BigInt(76130534848866565814734091037203300805469440) * 6^(1 // 2)
    A[13, 7] = -BigInt(9318570701703691930191534222789494006237103047579881711) //
               BigInt(3562968864769188061196125478115789029665196590638000000) +
               BigInt(6434757944786940764537553522803508307647498263055753667) //
               BigInt(7125937729538376122392250956231578059330393181276000000) * 6^(1 // 2)
    A[13, 8] = -BigInt(3761276943318114593157887595063197569838926771255702620302334773297) //
               BigInt(1464216542935037622889021617525045499738187320641028236799386000000) -
               BigInt(3559372342256314491923633965576963561357143961289) //
               BigInt(34078442106289230846820763513066043517445275500000) * 6^(1 // 2)
    A[13, 9] = BigInt(438696881717332693164574683605145782902979477404432064959580227217) //
               BigInt(190320793597985643524049406938231176153947641989529318077977593750) -
               BigInt(3559372342256314491923633965576963561357143961289) //
               BigInt(6904380415341463658460171786984220252917056671875) * 6^(1 // 2)
    A[13, 10] = BigInt(39393077397330623465044422844130393549390670768068628724675776122752) //
                BigInt(24104666116749847055078895726017945791587677982957405550119962671875) -
                BigInt(911199319617616509932450295187702671707428854089984) //
                BigInt(4383797390367410875150919593877143523760184037765625) * 6^(1 // 2)
    A[13, 11] = BigInt(1455014345065512890761705230599774567929341275785651) //
                BigInt(437274533681980646175086038705021590800143506216585)
    A[13, 12] = BigInt(43491235157516875762138510242422928845963515012067399050743561004230490804423974878600834753885909) //
                BigInt(15650073716933473532700213239415706508578020603344410725067680963543589421556571669667160279586515)
    A[14, 1] = BigInt(1113523038973067377822299106095464437335258246585030516055024186482995053) //
               BigInt(1172081821952666268571554825436238158838606320428900711229598548800375000) -
               BigInt(9656355719858106181793772673915704646989403181) //
               BigInt(42650123991481333817759653132492492758879187500) * 6^(1 // 2)
    A[14, 2] = 0
    A[14, 3] = 0
    A[14, 4] = 0
    A[14, 5] = 0
    A[14, 6] = BigInt(65156558806314041332627176859309009839731119) //
               BigInt(11498348446055040912040697694327646399463400) +
               BigInt(2893908603034923903437694020377944093049) //
               BigInt(44653780373029285095303680366320956891120) * 6^(1 // 2)
    A[14, 7] = BigInt(8209640140354849710164186531509948506760631106345347) //
               BigInt(2089832016538282470163966831380108687881472614625000) -
               BigInt(7326425889801057568096943993446608268878122947239059) //
               BigInt(4179664033076564940327933662760217375762945229250000) * 6^(1 // 2)
    A[14, 8] = BigInt(358387301461393123075253933751833574944722054371964074460462507967392441417) //
               BigInt(600381220734980367702253537950121782364502552531871330639897526611625875000) +
               BigInt(9656355719858106181793772673915704646989403181) //
               BigInt(39976896846716122011518307352108347913593312500) * 6^(1 // 2)
    A[14, 9] = -BigInt(642851404125311605101723817200220613720981758568362376903035228606467992291) //
               BigInt(274647144865445671171423994203344292862210985060954852952360191727428234375) +
               BigInt(77250845758864849454350181391325637175915225448) //
               BigInt(64795380679365872530827165335902056306758765625) * 6^(1 // 2)
    A[14, 10] = -BigInt(369969589301763518639761285997824520507173970793320519753635225871563412913152) //
                BigInt(180850870268987312974476662632186804296491491034759756207691410652763681078125) +
                BigInt(454852979828196233587213868032125351691788847437824) //
                BigInt(946232027174365750372191179393662736769436409578125) * 6^(1 // 2)
    A[14, 11] = -BigInt(7715303458281199041757411869319392732268755939408305903151657039) //
                BigInt(1483062857858019956599354292891821472831958379628421680653323310)
    A[14, 12] = -BigInt(945432908618140420875128019003679821034774098061354410087876969999388668430388640158852256997513450369741874915618358283) //
                BigInt(1183336256073681379271114843949311985569932050547818490887959987766109714104308279282406274069868573023638087923030436570)
    A[14, 13] = BigInt(1056911827593717127690972016166243945915857152155) //
                BigInt(7359790814129537930306952068888958882079118176267)
    A[15, 1] = BigInt(71647807109611556622880429107737597894249449144677487291739567) //
               BigInt(4608441750702552040544087416037317498817362012884252849000000) -
               751354288295256883750046307992744439 //
               260856700304180835921321781903000000 * 6^(1 // 2)
    A[15, 2] = 0
    A[15, 3] = 0
    A[15, 4] = 0
    A[15, 5] = 0
    A[15, 6] = BigInt(538527061627015039205734500467757677) //
               BigInt(6645826087749591758241674935867200) +
               BigInt(55735563892040165710948254310427) //
               BigInt(25809033350483851488317184216960) * 6^(1 // 2)
    A[15, 7] = BigInt(10166737271161511304691200043681679377006693) //
               BigInt(172554720727113306507896188662956931000000) -
               BigInt(8156608278080671316516960851836659322130571) //
               BigInt(345109441454226613015792377325913862000000) * 6^(1 // 2)
    A[15, 8] = BigInt(13730153329629135442521913725273883916343982389395960187482153016633) //
               BigInt(175471617700441378138284361033681728617635976027122756488121000000) +
               BigInt(20286565783971935861251250315804099853) //
               BigInt(6601681107698115001393451249699000000) * 6^(1 // 2)
    A[15, 9] = -BigInt(653480004921889631988496090206918830045333930618354572374567471359) //
               BigInt(12778455326756234342932466086484629156552273764794479235501062500) +
               751354288295256883750046307992744439 // 49537690682765110667751011467156250 *
               6^(1 // 2)
    A[15, 10] = -BigInt(1892005381291509192972446465852606492709607810610330309174017655027008) //
                BigInt(53004882228844808583820366645125601151772763811007845509116404390625) +
                BigInt(2211987024741236265760136330730639628416) //
                BigInt(361709036854355745656111465938134640625) * 6^(1 // 2)
    A[15, 11] = -BigInt(754694656733111129471952033526995126769302954114714957) //
                BigInt(10457763407796883263416278675753537185966936455426905)
    A[15, 12] = -BigInt(3809310391558597134511116925648017476134448177941143033030907811018093065572169647629059619278510306257177668207871699977487) //
                BigInt(49657222601249321414413002320858790760827862034335478005117886015451993139010886447467038877038198587307390592603767686355)
    A[15, 13] = BigInt(1329021384022690762809893094055902242790) //
                BigInt(314695772164147033634887576401647768207)
    A[15, 14] = -6893352883273198849981789770 // 5486160832491820655265629723
    A[16, 1] = -BigInt(2570623236993305537882133348905686861198442251997239696046711165213) //
               BigInt(1788444074612646395894349444415762174941041849960120845639920000000) +
               BigInt(17049571162138593122630888759407239103637) //
               BigInt(41413609740291749510869046094920280000000) * 6^(1 // 2)
    A[16, 2] = 0
    A[16, 3] = 0
    A[16, 4] = 0
    A[16, 5] = 0
    A[16, 6] = -BigInt(104754859887173191244557803812288138363) //
               BigInt(11164987827419314153846013892256896000) +
               3784288812661530935801759817347107 // 43359176028812870500372869484492800 *
               6^(1 // 2)
    A[16, 7] = -BigInt(12625417179938563922326639910886861363642731119) //
               BigInt(2029243515750852484532859178676373508560000000) +
               BigInt(4034430827439578769139663757402876281578669781) //
               BigInt(1352829010500568323021906119117582339040000000) * 6^(1 // 2)
    A[16, 8] = -BigInt(121185791557797038984865049921261438811435078875133421073861312913308541) //
               BigInt(761448556714003333948396687643999840432514259811137574414811792240000000) -
               BigInt(17049571162138593122630888759407239103637) //
               BigInt(38817884913264916208193493348230120000000) * 6^(1 // 2)
    A[16, 9] = BigInt(7831706161943944732903576802308126564075661689777008313099672046941217) //
               BigInt(1912115178728933570934305539114995143003765478134474309438025195000000) -
               BigInt(17049571162138593122630888759407239103637) //
               BigInt(7864603772795788969612150580525726250000) * 6^(1 // 2)
    A[16, 10] = BigInt(1231862339518930958878836935481237970791679023430521469669885456183927336) //
                BigInt(350627295943808408781971725357505851618976832609816898042805015043984375) -
                BigInt(6274242187667002269128167063461863990138416) //
                BigInt(7178115836374689772545532041542281943203125) * 6^(1 // 2)
    A[16, 11] = BigInt(23975119191644794266115435124846718653445865470965352804548709) //
                BigInt(2618162142480251258794523717542364435163988588525025359875200)
    A[16, 12] = BigInt(54918621392363748304153623111076473105288800876284517914641326502791860135537567319289496873218100975783896487633335215553248687) //
                BigInt(42045763520929825428011777325117555313008167341712535936493416447003511630663297772799291157965783507844913762569462175390505600)
    A[16, 13] = -BigInt(35648602464057004250287591332213307928385) //
                BigInt(140983705929537871068429634227938200156736)
    A[16, 14] = 180996536632507119654644161465573 // 475057638807131734180761408974016
    A[16, 15] = 0
    α[1] = 385924436255198461459913 // 25885297292164750617319296
    α[2] = 0
    α[3] = 0
    α[4] = 0
    α[5] = 0
    α[6] = 0
    α[7] = 0
    α[8] = -81508791888782942071778080019859673 // 399395400777999552787672219983509760
    α[9] = BigInt(935936315524449978576662571821361001) //
           BigInt(4086801409502444670525951816974611200)
    α[10] = 2172547024243858864854526674870272 // 16972283408414027681310867548214525
    α[11] = BigInt(1193746724713997342094811077921918673563219) //
            BigInt(5333839543124606397779553737025343484595200)
    α[12] = BigInt(78998885900843720607956191287550920802893766832745691186069543513026890766039702331356142187984337094261434894685526783) //
            BigInt(199728353759405693383344357357050366806065734043901244628659238274456345277244226082448918708806038102409307770530790400)
    α[13] = 11409994679937666036993318622713183 // 210646829819316523523317257829664000
    α[14] = 10048608923923592706010638721995991 // 79176273134521955696793568162336000
    α[15] = -491671864801784912209409 // 935818304406566046193728000
    α[16] = 1 // 31
    αEEst[1] = 252989724430399818032879507 // 2563363467960203776409535840
    αEEst[2] = 0
    αEEst[3] = 0
    αEEst[4] = 0
    αEEst[5] = 0
    αEEst[6] = 0
    αEEst[7] = 0
    αEEst[8] = BigInt(29177952541973984236460298638317525091) //
               BigInt(1395926082130949417341226925726678720)
    αEEst[9] = -481359954433374971820713358519126547 // 262512559006778653160811139234495296
    αEEst[10] = -325585565780649471399541319800832 // 3935361974048520890008015792480311
    αEEst[11] = BigInt(828906176711503347732519444777342468824091) //
                BigInt(1234749218911637779745721368343853865102720)
    αEEst[12] = -BigInt(19819005742497004038604841557068415948625832456498098356731269808855703554130425034149167680385863035602494067) //
                BigInt(944246035063763492934067294951661726695697507391691339725080641234241699449384342323015093676255343392536960)
    αEEst[13] = 96762135771288342741313072067956263 // 15762193817514374346399946534150720
    αEEst[14] = -250908066606160220133619476267717281 // 58399326977497401115669462861805760
    αEEst[15] = 33 // 92
    αEEst[16] = 1 // 31

    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 9, αEEst = αEEst, adaptiveorder = 8,
        stability_size = 4.501341330915361))
end

"""
Tom Baker, University of Teeside.
Part of RK-Aid
http://www.scm.tees.ac.uk/users/u0000251/research/researcht.htm
http://www.scm.tees.ac.uk/users/u0000251/j.r.dormand/t.baker/rk10921m/rk10921m
"""
function constructBaker10(T::Type = BigFloat)
    A = zeros(T, 21, 21)
    c = zeros(T, 21)
    α = zeros(T, 21)
    αEEst = zeros(T, 21)

    c[2] = convert(T,
        big".2232129192123735665527132860096509572180113541116006004801774831897914255866976731834")
    c[3] = convert(T,
        big".3348193788185603498290699290144764358270170311674009007202662247846871383800465097751")
    c[4] = convert(T,
        big".5022290682278405247436048935217146537405255467511013510803993371770307075700697646627")
    c[5] = convert(T,
        big".1176948756548443596524604744040430749724814058786280779261159420482531581524139111558")
    c[6] = convert(T,
        big".6425923677604462086219057105219705193480713109262431864517703346402040737649802005250")
    c[7] = convert(T, big".181826565311210497")
    c[8] = convert(T,
        big".4341610334077954337462856744299909865481150413847250582144329859884319977831972426112")
    c[9] = convert(T,
        big".7122335424182418496009012272374473135415214567547251841147975874613222331084545787886")
    c[10] = convert(T,
        big".1894400309592476109308204841234726559856268372193596062465558450034952083309882121192")
    c[11] = convert(T,
        big".4943921538691045853333333333333333333333333333333333333333333333333333333333333333333")
    c[12] = convert(T,
        big".6403413702953303301268868454436798781085931758963391614454137298592697657577711636312")
    c[13] = convert(T, big".741588230803656878")
    c[14] = convert(T, big".382911329862520855")
    c[15] = convert(T, big".107157755822487322")
    c[16] = convert(T, big".875691376241245711")
    c[17] = convert(T, big".964069299370187816")
    c[18] = convert(T, big".281729610523717922")
    c[19] = convert(T, big".631145492177176748")
    c[20] = convert(T, big".973803039377034607")
    c[21] = convert(T, big"1.")
    A[2, 1] = convert(T,
        big".2232129192123735665527132860096509572180113541116006004801774831897914255866976731834")
    A[3, 1] = convert(T,
        big".8370484470464008745726748225361910895675425779185022518006655619617178459501162744378e-1")
    A[3, 2] = convert(T,
        big".2511145341139202623718024467608573268702627733755506755401996685885153537850348823314")
    A[4, 1] = convert(T,
        big".1255572670569601311859012233804286634351313866877753377700998342942576768925174411657")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T,
        big".3766718011708803935577036701412859903053941600633260133102995028827730306775523234970")
    A[5, 1] = convert(T,
        big".8645012631860283043095827588443080342589988168209518203030156172654040053726457363040e-1")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T,
        big".5236243791730284147504959465381950604441622286481307774905930763200738501679225239292e-1")
    A[5, 4] = convert(T,
        big"-.2111768858106131225354739613420723449783469866828018185324492731029462740164291486751e-1")
    A[6, 1] = convert(T,
        big"-.2639257721366917655020719039008537814312675653905044467083967854693494521584367325145e-1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T,
        big".3321586951309942525210016931443520748342752809336835485070926725007327816351124496012")
    A[6, 5] = convert(T,
        big".3368262498431211326511112077677038226569227865316100826155173406864062373457114241753")
    A[7, 1] = convert(T,
        big".4214481430128891924590207161639307925193144098983055869572587559196013842334308634337e-1")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T, big"0.")
    A[7, 5] = convert(T,
        big".1395091008993089985837339976891755843541033762656629420011519949721775047502000865276")
    A[7, 6] = convert(T,
        big".1726501106125791703639306944313363939651827445064993031221294358623568264568271289984e-3")
    A[8, 1] = convert(T,
        big".1913223517063317850561878098464623128757461736096094048517648270556161690028960671678")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big"-.6549476653170516368338607340883454058106339891146285998978549868910881698171713079544")
    A[8, 6] = convert(T,
        big".1755876323742502725238471568246235150972597859217101866902108545773998108136024750683e-1")
    A[8, 7] = convert(T,
        big".8802275837810902582715738829894117279732768782975732345915020603661640175161122358910")
    A[9, 1] = convert(T,
        big".5819519506251173970176101222224200992345451475572200041442763494089781743712740682029e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big".1647908558976784212367930935320516004070340774715814512653291320117464889726338996259")
    A[9, 7] = convert(T,
        big".2562758618556534830885936612274458691344630164445656645146446091037555911278059725788")
    A[9, 8] = convert(T,
        big".2329716296023982055737534602557078340765698480828560679203962114049223355708872997637")
    A[10, 1] = convert(T,
        big".6609808969751805186432114363691927852621669632797821827493854746683225501745446560252e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big".4993798303607936315852000088060192246997357096577012552548527217200430403080630524414e-1")
    A[10, 7] = convert(T,
        big".1536216094454890707317822140262379594167730891287251273659790090180235843883503581847")
    A[10, 8] = convert(T,
        big"-.5414130381994316899529893137177555939937633423733000345305400547297139497578508331128e-1")
    A[10, 9] = convert(T,
        big"-.2607634739989570582850394304851094502796018496578386146679297818039354012983783360089e-1")
    A[11, 1] = convert(T,
        big".5977281413067455622264838600535046960155059547948827360548395666199015131925709118733e-1")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T, big"0.")
    A[11, 7] = convert(T, big"0.")
    A[11, 8] = convert(T,
        big".1656649460691327819663527933964608353435263909407310250007028968559622442637053740488")
    A[11, 9] = convert(T,
        big"-.1270467077219322361480177163110966635379558631317440513037999854071349138024611223099e-2")
    A[11, 10] = convert(T,
        big".2702248607465165695058123310946329950236359055444314752401844796694522868883954793203")
    A[12, 1] = convert(T,
        big".5714646346927726543617292683289117330020257083064077657640682066128476784510949443172e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T, big"0.")
    A[12, 7] = convert(T, big"0.")
    A[12, 8] = convert(T, big"0.")
    A[12, 9] = convert(T,
        big".1706647333279140031893907789070497592665517669717541507696626200796144815409471923569e-1")
    A[12, 10] = convert(T,
        big".2853767151165243146804641383356191775396044057521502467001987179873305495674085297549")
    A[12, 11] = convert(T,
        big".2807517183767373496913107023844645513421310226163727230918419292026930001911584202089")
    A[13, 1] = convert(T,
        big".5713690770667411245888754904665701448619979494564888202366599688508571201235483609660e-1")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T, big"0.")
    A[13, 5] = convert(T, big"0.")
    A[13, 6] = convert(T, big"0.")
    A[13, 7] = convert(T, big"0.")
    A[13, 8] = convert(T, big"0.")
    A[13, 9] = convert(T,
        big".8662146075905090149738288126241361791710263955784173047391548833627524888649185757901e-1")
    A[13, 10] = convert(T,
        big".2854329582389984323634810972884390808440109342965556282630983683305184760043391707767")
    A[13, 11] = convert(T,
        big".2797633923959254653770985744082321055561046971559230633717557229339774187704436927395")
    A[13, 12] = convert(T,
        big".3263351170300796630314989799425818119658193404403069586756442351414314432637044280815e-1")
    A[14, 1] = convert(T,
        big".5792154217958510846632846508298421800727675597066414244319766319648769731268367294802e-1")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T, big"0.")
    A[14, 5] = convert(T, big"0.")
    A[14, 6] = convert(T, big"0.")
    A[14, 7] = convert(T, big"0.")
    A[14, 8] = convert(T, big"0.")
    A[14, 9] = convert(T,
        big"-.1168598764394854203108670771355968386139754171212251879659829443523204385801129100146")
    A[14, 10] = convert(T,
        big".2791270024477207066381851603610170508129292534731657465038373367821947287200953472810")
    A[14, 11] = convert(T,
        big".5523482312158897030618178268050450731680950103624425704577340542168627357243106364454e-1")
    A[14, 12] = convert(T,
        big".3313821420230448163586175200236513161592367149984143944502756986720155621978125806284e-1")
    A[14, 13] = convert(T,
        big".7434962435080700826430991700872593086103623514130960252814696908475018275512156807813e-1")
    A[15, 1] = convert(T,
        big".5424349131266048906033564077410202025889503282676783408313456162421294081865905389248e-1")
    A[15, 2] = convert(T, big"0.")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T, big"0.")
    A[15, 7] = convert(T, big"0.")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T,
        big".2555744036484019017045493997261865838594175905367142496400185786785191604031397243218e-1")
    A[15, 10] = convert(T,
        big".1065318948685044129706483445621162791089602467505369099004177225283624953902959013228")
    A[15, 11] = convert(T,
        big".1494739176822617600298548318084537294306548147599349812604239821749187453609545343658")
    A[15, 12] = convert(T,
        big"-.7541418150821686840498891896239704084250667222662001839556328840119093341810397465136e-1")
    A[15, 13] = convert(T,
        big".1183662716974706173695161845106353658054818835708868187585164205844835807880512638150e-2")
    A[15, 14] = convert(T, big"-.154418469614537368")
    A[16, 1] = convert(T,
        big"-.3130073355625165399871063999971367682939799979170688711523373232868473829204122305445e-1")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T, big"0.")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T,
        big"-.4393041016682237597697452523147415408933247130414835042917723489831082747083637622576")
    A[16, 10] = convert(T,
        big"-.2971894957150042139937620479069507765251785957115322181603409909075380714667735160813")
    A[16, 11] = convert(T,
        big".1262980195755398451585609791771777378208362550125140868265029793320289593624402155837")
    A[16, 12] = convert(T,
        big"-.1352235663389501393036813956309918309028161801076417057198393439731017100387729485379e-1")
    A[16, 13] = convert(T,
        big".6695605432945227295340251006073274395173466715429726933128280272846122961086155806635")
    A[16, 14] = convert(T, big".398926850454282556")
    A[16, 15] = convert(T, big".462222650490275222")
    A[17, 1] = convert(T,
        big".4254608456033344945622257854593736284131499735662276010271273463216872625702497770866e-1")
    A[17, 2] = convert(T, big"0.")
    A[17, 3] = convert(T, big"0.")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T, big"0.")
    A[17, 7] = convert(T, big"0.")
    A[17, 8] = convert(T, big"0.")
    A[17, 9] = convert(T,
        big"-.1821182296512653493976221650402057582117417514147759470622260585285199838506585092559")
    A[17, 10] = convert(T,
        big".1693488568788993457692776359468516487310835997985566768493328170663073221084671579068")
    A[17, 11] = convert(T,
        big"-.2812932311380423097505532359005039285971305910469541561543126052387690343228512975388e-2")
    A[17, 12] = convert(T,
        big".4367359500435692592480947784835035112081291953371300785838045346974343617081011595250")
    A[17, 13] = convert(T,
        big".9649623253501190961728232014648167901724237884864786722867265364689659268793239634042e-2")
    A[17, 14] = convert(T,
        big".1827541189026459104590716992261562936435133038172843022420870171219768684076842803782")
    A[17, 15] = convert(T,
        big".8163992629788446963341645082372801271393454416361108721409416282456258720143539344240e-1")
    A[17, 16] = convert(T,
        big".2263259013959999629673163223583858004580131789671757969088706528737681492423808136362")
    A[18, 1] = convert(T,
        big"-.6162786199706763955052568808743247751165729542444646544961722501628952288996567177280e-1")
    A[18, 2] = convert(T, big"0.")
    A[18, 3] = convert(T, big"0.")
    A[18, 4] = convert(T, big"0.")
    A[18, 5] = convert(T, big"0.")
    A[18, 6] = convert(T, big"0.")
    A[18, 7] = convert(T, big"0.")
    A[18, 8] = convert(T, big"0.")
    A[18, 9] = convert(T,
        big".2630819816424975496840706912577088603030984130761613342651746722632362416592729338065e-1")
    A[18, 10] = convert(T,
        big"-.6563819353303824389148917307603424726525093364227772615481314402265589817962663323807")
    A[18, 11] = convert(T,
        big"-.4996969738639946601661934398570589720521936732243676171454742217015940903119194293034")
    A[18, 12] = convert(T,
        big"-.7123015625596798250364056005999154454586809433543859495632863935368185749281956104361e-1")
    A[18, 13] = convert(T,
        big".1945887064918253295880838290230513152388881501421435771759995471289411971821731631008")
    A[18, 14] = convert(T,
        big".7502438829684451516052952921664213906375154664124515160998662893556068266252405773932")
    A[18, 15] = convert(T,
        big".6778525665536439867432200838066450041417881571145524785032579156045245492715999778952")
    A[18, 16] = convert(T,
        big"-.1046057647455310786652080097009039359249473747848876107492887123406055065619567484248")
    A[18, 17] = convert(T,
        big".2627894853849749889545315434384080663867415921515384464319901932333376180798673115539e-1")
    A[19, 1] = convert(T,
        big".3796521632874261201186203436841209430559695655768644122346210733272382403437031805294e-1")
    A[19, 2] = convert(T, big"0.")
    A[19, 3] = convert(T, big"0.")
    A[19, 4] = convert(T, big"0.")
    A[19, 5] = convert(T, big"0.")
    A[19, 6] = convert(T, big"0.")
    A[19, 7] = convert(T, big"0.")
    A[19, 8] = convert(T, big"0.")
    A[19, 9] = convert(T,
        big"-.2594056416058702079475722192315214454017994671360402614116801042085357278273610676615")
    A[19, 10] = convert(T,
        big".6023088303591701879982724849459819049819908747310071139535368513067871259409896536125")
    A[19, 11] = convert(T,
        big".6119387226913542078978281038711879095845432514371708320113689988639552031445275756169")
    A[19, 12] = convert(T,
        big"-.8557593281472216275173535332798856826091404689475696899046474643882545069592793853088")
    A[19, 13] = convert(T,
        big"1.016014044479150879819181655935044451622873646716581194764492422092270586935018600847")
    A[19, 14] = convert(T,
        big".4174635465318728144170599907012271233953357864579855698710555478932012900105054820193")
    A[19, 15] = convert(T,
        big"-.3142068748294011550388298801522627398568106507178091710884531834553635025761917127745e-2")
    A[19, 16] = convert(T,
        big"-.3113859536805386540643212859648902182769708893494791481799285293081168071417215791223")
    A[19, 17] = convert(T,
        big".7475323182002611486490630088842356498085791607202230555314272159767158496743797941577e-1")
    A[19, 18] = convert(T,
        big"-.6996051078512155679294752334324570751847195000321862661699180193471489380787256603442")
    A[20, 1] = convert(T,
        big".6354247858004162057833209973904648790423789887640110847597201944718467473335513915457e-1")
    A[20, 2] = convert(T, big"0.")
    A[20, 3] = convert(T, big"0.")
    A[20, 4] = convert(T, big"0.")
    A[20, 5] = convert(T, big"0.")
    A[20, 6] = convert(T, big"0.")
    A[20, 7] = convert(T, big"0.")
    A[20, 8] = convert(T, big"0.")
    A[20, 9] = convert(T,
        big".1335709427214855386564319367816597150656631989678512452281616080142611136536750546091")
    A[20, 10] = convert(T,
        big".3491939710685402856556370077059343902437022368395394555525545826876029022896488736825")
    A[20, 11] = convert(T,
        big"-.1284299219085288104043584513269729325317465049017238832115899338924554224999281557779")
    A[20, 12] = convert(T,
        big".5820143454609020234416665664957204705761540086346612646427774305741185500271652747165")
    A[20, 13] = convert(T,
        big"-.7032295749061967719609757517752619736970110102217160910226632819061710474763494562046e-1")
    A[20, 14] = convert(T,
        big".4108572358925710179476579003323188290642373757595878674150744283920916547885055706399")
    A[20, 15] = convert(T,
        big"-.2510654779239962466269627907503763505655390551209676232633905831376503530736867220512e-1")
    A[20, 16] = convert(T,
        big".1223166254706007785807460575472930729088496060945489634909944334258869465746143128690")
    A[20, 17] = convert(T,
        big".4712374572639742166850752217533416565767093458666637954607896698933218387654817641078e-1")
    A[20, 18] = convert(T,
        big"-.2035773172368828280351383367281992122282011746329188150975183458730713499142844139626")
    A[20, 19] = convert(T,
        big"-.3073795611150731392306884484695711542343125736903452146138998032605691134742962145162")
    A[21, 1] = convert(T,
        big".1189924137880676460852663367282125023346761041034714988771622602337896848257576311721")
    A[21, 2] = convert(T, big"0.")
    A[21, 3] = convert(T, big"0.")
    A[21, 4] = convert(T, big"0.")
    A[21, 5] = convert(T, big"0.")
    A[21, 6] = convert(T, big"0.")
    A[21, 7] = convert(T, big"0.")
    A[21, 8] = convert(T, big"0.")
    A[21, 9] = convert(T,
        big"1.567065389218647018255284948985537561396036403483411695550932448674279128530361565244")
    A[21, 10] = convert(T,
        big"1.062742213102872753657741254837796588064400661866785937264672961506012193460062190291")
    A[21, 11] = convert(T,
        big"-.7051696022721412819771965938099389611228930171636494130149946862206926743626267591091")
    A[21, 12] = convert(T,
        big".9907825344767892999287937933800084396119601938859814415204705439463102044422773242413")
    A[21, 13] = convert(T,
        big".2507980962931107922266871779870423357231706420300879746540980799200371363012966217081")
    A[21, 14] = convert(T,
        big"1.972591757153629542226038011120263810650524345125694323573613060238829589123920677616")
    A[21, 15] = convert(T,
        big"-.3262628387941532893880299499584174546880007304318604226454131079685364715701264155684")
    A[21, 16] = convert(T,
        big"-.7195961905900189101450624277119790470362703727002334732811111315403329989928142118605")
    A[21, 17] = convert(T,
        big".5802031535291661528133435810283051023192180079596550160030993013203826047471987749923")
    A[21, 18] = convert(T,
        big"-1.503104563178041688806298917914932279967077722444117887454135595533548208509663052219")
    A[21, 19] = convert(T,
        big"-2.082917406986733532241895178393857460921252779279252368270349439324592549474879051363")
    A[21, 20] = convert(T,
        big"-.2061249557411945026346720362780411363644917364359743227780446952519376385207652951444")
    α[1] = convert(T,
        big".3074440935793207675883517938984149786411995438646211366049687901929864736726360396701e-1")
    α[2] = convert(T, big"0.")
    α[3] = convert(T, big"0.")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big"0.")
    α[6] = convert(T, big"0.")
    α[7] = convert(T, big"0.")
    α[8] = convert(T, big"0.")
    α[9] = convert(T, big"0.")
    α[10] = convert(T, big"0.")
    α[11] = convert(T, big"0.")
    α[12] = convert(T,
        big".4447886124335158463373029422601699728444313324875120673928896655863659883998986985001")
    α[13] = convert(T,
        big"-.7714472898672551027236427412075689933241937304249667830486089002567926592070578328099e-1")
    α[14] = convert(T,
        big".1895969236666960081782143661904924178003392273539755384045555627359634132139704039829")
    α[15] = convert(T,
        big".1696794742660741325357704245298007854728500376148224262050684467711119689665579207267")
    α[16] = convert(T,
        big".3107725327623489988155072749206969138501610275590999188855016092183529050525521871700")
    α[17] = convert(T,
        big"-1.052173416444910172223599400939351222118848154980691993896969290833000795415995203779")
    α[18] = convert(T,
        big".1098353169596310609348506603862869818226618225743191040640071521542528953059954349550")
    α[19] = convert(T,
        big"-.1334572232187388442842841059278351057902328711498735111118805829160094435254500793407")
    α[20] = convert(T,
        big"1.154026328574134208065260150796827370454610532524502151460002781580501491609214462821")
    α[21] = convert(T,
        big"-.1466682293699578048454932174861727128676735353276311367588113332911578050533016457222")
    αEEst[1] = convert(T,
        big".3069935889151056339784198488152917237594701959064431113839566366053496143740989806580e-1")
    αEEst[2] = convert(T, big"0.")
    αEEst[3] = convert(T, big"0.")
    αEEst[4] = convert(T, big"0.")
    αEEst[5] = convert(T, big"0.")
    αEEst[6] = convert(T, big"0.")
    αEEst[7] = convert(T, big"0.")
    αEEst[8] = convert(T, big"0.")
    αEEst[9] = convert(T, big"0.")
    αEEst[10] = convert(T, big"0.")
    αEEst[11] = convert(T, big"0.")
    αEEst[12] = convert(T,
        big".4501037370765134863000400416812391552231563060566733706571606677666390164501322407860")
    αEEst[13] = convert(T,
        big"-.6514757933632400121177566572810133388907129562810217527103583351082868427677317818171e-1")
    αEEst[14] = convert(T,
        big".1916476309177840038043731680068457542540807928558914659725060136615214781805707317846")
    αEEst[15] = convert(T,
        big".1699093037015861786030207592719576231809404332930270025745281206261397945871067257969")
    αEEst[16] = convert(T,
        big".2831822563446878149008449047225439184357774091336334653238641676396106848038134341285")
    αEEst[17] = convert(T,
        big"-.7938273328067640993878535581971311480056976834520424653124835645605083935506842817243")
    αEEst[18] = convert(T,
        big".1085805168491408111596918257428235537727685528809234551377846603272747280143198404809")
    αEEst[19] = convert(T,
        big"-.1442902422009048623094173280045594665841429401482516755746485620108546739010869727583")
    αEEst[20] = convert(T,
        big".8718101111217405681350791198631736702436128801469513492598526664004710882551915616216")
    αEEst[21] = convert(T,
        big"-.102667760558970463391845252240320899007371474729348103905924")

    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 10, αEEst = αEEst, adaptiveorder = 9,
        stability_size = 5.827746350596612))
end

"""
Ono10
"""
function constructOno10(T::Type = BigFloat)
    A = zeros(T, 17, 17)
    c = zeros(T, 17)
    α = zeros(T, 17)

    c[2] = convert(T,
        big".3357505083417036184129939488963472892525539577157696170900805997207182929518116563365")
    c[3] = convert(T,
        big".5263563553500217821118595221339767434067222948292617539626688113669112965383197614021")
    c[4] = convert(T,
        big".7895345330250326731677892832009651151100834422438926309440032170503669448074796421031")
    c[5] = convert(T,
        big".1852155685265047076970656199437875208987680209449849941608024849050604421517806064228")
    c[6] = convert(T,
        big".2895345330250326731677892832009651151100834422438926309440032170503669448074796421031")
    c[7] = convert(T,
        big".7659879027055932401898717206953675811356898851157703878511505293452324618973555436361")
    c[8] = convert(T,
        big".1080739009578824490100240661758266719014599665988034811870817465871751958605805741080")
    c[9] = convert(T,
        big".3573842417596774518429245029795604640404982636367873040901247917361510345429002009092")
    c[10] = convert(T,
        big".8825276619647323464255014869796690751828678442680521196637911779185276585194132570617")
    c[11] = convert(T,
        big".6426157582403225481570754970204395359595017363632126959098752082638489654570997990908")
    c[12] = convert(T,
        big".1174723380352676535744985130203309248171321557319478803362088220814723414805867429383")
    c[13] = convert(T,
        big".7659879027055932401898717206953675811356898851157703878511505293452324618973555436361")
    c[14] = convert(T,
        big".2895345330250326731677892832009651151100834422438926309440032170503669448074796421031")
    c[15] = convert(T,
        big".5263563553500217821118595221339767434067222948292617539626688113669112965383197614021")
    c[16] = convert(T,
        big".3357505083417036184129939488963472892525539577157696170900805997207182929518116563365")
    c[17] = convert(T, big"1.")
    A[2, 1] = convert(T,
        big".3357505083417036184129939488963472892525539577157696170900805997207182929518116563365")
    A[3, 1] = convert(T,
        big".1137717040478783056449396184425018992536986324028992632761012508174661576282091865815")
    A[3, 2] = convert(T,
        big".4125846513021434764669199036914748441530236624263624906865675605494451389101105748206")
    A[4, 1] = convert(T,
        big".1973836332562581682919473208002412787775208605609731577360008042625917362018699105258")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T,
        big".5921508997687745048758419624007238363325625816829194732080024127877752086056097315774")
    A[5, 1] = convert(T,
        big".1360001717992528326665261557210531883399922632350163155393510462488093550801937826726")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T,
        big".8247208052028112223790799775745048180379960003630910417779053602715946236310693542837e-1")
    A[5, 4] = convert(T,
        big"-.3325668379302924720736853353471614924502384232634042555633909737090837529152011167815e-1")
    A[6, 1] = convert(T,
        big".6546785199481110579756630335868870724165088583024509548938434492848601278421108871818e-1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T,
        big".6858715620423033993966826640550281985093888681236135859756994188871385539662876478859e-3")
    A[6, 5] = convert(T,
        big".2233808094681792639708262971782213796699231675455239218686431727029937934693022657371")
    A[7, 1] = convert(T,
        big".2379439999135994223360182301765245366877534054337310606780367192080205084533052094979")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".1285793626493546102687479560301629068933431333159680850577212924587597635005521433930")
    A[7, 5] = convert(T,
        big"-.7303763776380165799305971651206183012312272113455920027358381669467410091843571268521")
    A[7, 6] = convert(T,
        big"1.129840917780655787515702699609298438785820557711663244851230684625193199127855317597")
    A[8, 1] = convert(T,
        big".6062780896441720418288517632894493561013839292138719356621888281245740094338544165775e-1")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big".7888222248692022386398029561573557081956734138939980123787302990185683517283562747179e-1")
    A[8, 6] = convert(T,
        big"-.3213213504125939005101626782074483911778772517654658420946839341553511609992151857815e-1")
    A[8, 7] = convert(T,
        big".6960045478044110141748620518910045895419574645630705924582272883960758442810235566394e-3")
    A[9, 1] = convert(T,
        big".3104415865438721889084097538819995930358479868352200893591644921124828667005686536920e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big".1571656120644442171909172634188670081467914525944196307657031602142991040495820182979")
    A[9, 7] = convert(T,
        big".1117638541384084302989687338271505642661663918417037071041752801553969715970722155609e-3")
    A[9, 8] = convert(T,
        big".1690627071867076073308672954386663460258558459670039606814010070304482468516642450265")
    A[10, 1] = convert(T,
        big".1430614201422677732943411456202705981582929884075538519224332811759391256045628709601e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big"-.3914725335338579386439799990026912235591541327277623987271673700062264924050571982953")
    A[10, 7] = convert(T,
        big".2848400967322984698670378405298733357388045123736554424201702352859211869408079104816")
    A[10, 8] = convert(T,
        big".2559426777170804941167590501634194225203620020234561863154326261989141389268124787652")
    A[10, 9] = convert(T,
        big".7189112790349845437562504807270404806670261637579475044631123583223249124963937790142")
    A[11, 1] = convert(T,
        big".1006668457426631788172270822620227260906121023931371723741841060852730989522549641648e-1")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T,
        big"-.3683271809355001194608455001594963737510628435069995678571139196388186857259494651043")
    A[11, 7] = convert(T,
        big".8575502624022299903546939609254513705068159712672738891650679949767497156363396523266e-1")
    A[11, 8] = convert(T,
        big".2633792404483482304677774849110341109669789276529274236259927084605103685123897043694")
    A[11, 9] = convert(T,
        big".6783120744401823352568978305425238981627472606453781308422134412090277276134222894570")
    A[11, 10] = convert(T,
        big"-.2657008652719721502394642259236950907890441579413439685514223187307272640162219128045e-1")
    A[12, 1] = convert(T,
        big".4212573015701458022875431871453412827491044618882839077863171149554249469159410384126e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T,
        big"-.1140730257394986628323587415013320453327540775054758894945046140858856585220720793160e-1")
    A[12, 7] = convert(T,
        big"-.3577323583881041642639854910390843838477094752248779125994771974507179107770816673622e-2")
    A[12, 8] = convert(T,
        big".8478085069047181595883795669333443274198965196004767054783136203447118275436721594692e-1")
    A[12, 9] = convert(T,
        big"-.8236897403845190534516723638464155343698796681468739910260266595384301707513187151823e-4")
    A[12, 10] = convert(T,
        big".7920174936649160445251535005941104644271462692404428580786438040911208546142377592512e-3")
    A[12, 11] = convert(T,
        big".4840734825985701173601980408776943260994401783442431626214940796417131157064341867563e-2")
    A[13, 1] = convert(T,
        big".1866168584194642618748985099397434698267694490527249444362505121514680479505442836574")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T,
        big".1285793626493546102687479560301629068933431333159680850577212924587597635005521433930")
    A[13, 5] = convert(T,
        big"-.7303763776380165799305971651206183012312272113455920027358381669467410091843571268521")
    A[13, 6] = convert(T,
        big".3365620746651354218129539461265307620933082614096471567813938850533189571406833598679")
    A[13, 7] = convert(T,
        big".3439152195696137654283681672732838006806635738118480999635325742827372086353048385257")
    A[13, 8] = convert(T,
        big".6116809614283795957124860804122425718138169583157463948158033894717678649195939132136")
    A[13, 9] = convert(T,
        big".8555642651045021202243366091904917631763787495189342130253621820535903092162656518481")
    A[13, 10] = convert(T,
        big"-.8913928524866960048456401557668818573165120479499861285407465844438950550960191679472e-1")
    A[13, 11] = convert(T,
        big"-.4263317531098201582091498458499961586420389075833366648914840378325891404712881232497")
    A[13, 12] = convert(T,
        big"-.4510834231343501965076085217297850477436729165851712257475164429026900343003414799732")
    A[14, 1] = convert(T,
        big".8132722174581177948072809791683141502457873136125978954741528825822682868334020231374e-1")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T,
        big".6858715620423033993966826640550281985093888681236135859756994188871385539662876478859e-3")
    A[14, 5] = convert(T,
        big".2233808094681792639708262971782213796699231675455239218686431727029937934693022657371")
    A[14, 6] = convert(T,
        big"-.3543152669595614715438074770044728248424634682417120916726014150358797064660842419218")
    A[14, 7] = convert(T,
        big"-.1614228337471357340842753214971350545154909424968569413717974406833347832344548170007")
    A[14, 8] = convert(T,
        big"-1.267989518319081979614362221100086229540210555088918164010191288219362963604567702803")
    A[14, 9] = convert(T,
        big".2482505660965057784745955362600844051763121832196360142449926255422280563592831868088")
    A[14, 10] = convert(T,
        big".1894015072204329862304305088267495664543558350848310970103331708369628629196601961145e-2")
    A[14, 11] = convert(T,
        big"-.2367108504250404697299275875477514623653536202773466252372015412173547424356865161882e-1")
    A[14, 12] = convert(T,
        big"1.376373544047006195976245547649688122109906103559164922810806152088122419417100442009")
    A[14, 13] = convert(T,
        big".1650212091015662542191305948002865244010106371945579174943772453918520072439660689694")
    A[15, 1] = convert(T,
        big".1137717040478783056449396184425018992536986324028992632761012508174661576282091865815")
    A[15, 2] = convert(T,
        big".4125846513021434764669199036914748441530236624263624906865675605494451389101105748206")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T,
        big"-.5266141894222268880524004361760397955773267936583064281992495328485812564705031059898")
    A[15, 7] = convert(T,
        big".4270299995350464166415042303891280056028480962227556089080302280248036279775533169729")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T, big"0.")
    A[15, 10] = convert(T, big"0.")
    A[15, 11] = convert(T, big"0.")
    A[15, 12] = convert(T, big"0.")
    A[15, 13] = convert(T,
        big"-.4270299995350464166415042303891280056028480962227556089080302280248036279775533169729")
    A[15, 14] = convert(T,
        big".5266141894222268880524004361760397955773267936583064281992495328485812564705031059898")
    A[16, 1] = convert(T,
        big".3357505083417036184129939488963472892525539577157696170900805997207182929518116563365")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T,
        big"-.4368367083891407610190158064846822677554138099678168504498524386288907219928146135373")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T, big"0.")
    A[16, 10] = convert(T, big"0.")
    A[16, 11] = convert(T, big"0.")
    A[16, 12] = convert(T, big"0.")
    A[16, 13] = convert(T, big"0.")
    A[16, 14] = convert(T, big"0.")
    A[16, 15] = convert(T,
        big".4368367083891407610190158064846822677554138099678168504498524386288907219928146135373")
    A[17, 1] = convert(T,
        big".3528400914539065757695374250020526983792042094778706633574048575138598395445375459964e-1")
    A[17, 2] = convert(T,
        big"-.4368588562339554617519938616027395568153523718226810275461527253101515442644071737416")
    A[17, 3] = convert(T,
        big"-.5185253025751911700815354370658368214653404020923461413204927007907407344767758710502")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T,
        big".8353882146318347708823979907871444293693024897533008980286851134950209556868285481998e-1")
    A[17, 7] = convert(T,
        big".3357324883823615793514438806416832609180686549969009829861766772233990660034579152038")
    A[17, 8] = convert(T,
        big"-.1180346853290197678148440843338719378456952414935754110833043790640005049977678056530")
    A[17, 9] = convert(T,
        big"-.2028121524999718341813208304021045923675061150408674489110338814798518734278169865166")
    A[17, 10] = convert(T,
        big".4018734442526045674154192467995252408036818247069631361197899796078897150680466096292")
    A[17, 11] = convert(T,
        big".6982808681238487516481594542254805552510541794271997632640921307356093400075120505480")
    A[17, 12] = convert(T,
        big".2745953340127460413630142920026546525788337913278252157685398122819348957290706395825")
    A[17, 13] = convert(T,
        big"-.8071051158813288531951186779278104056732084407804403857002131168492985419909905100904")
    A[17, 14] = convert(T,
        big".2986469883301853807480531774155235135599206769328769914173437804434298240853514778770")
    A[17, 15] = convert(T,
        big".5185253025751911700815354370658368214653404020923461413204927007907407344767758710502")
    A[17, 16] = convert(T,
        big".4368588562339554617519938616027395568153523718226810275461527253101515442644071737416")
    α[1] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    α[2] = convert(T,
        big"-.2192242833052276559865092748735244519392917369308600337268128161888701517706576728499e-1")
    α[3] = convert(T,
        big"-.5671077504725897920604914933837429111531190926275992438563327032136105860113421550095e-1")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big"0.")
    α[6] = convert(T,
        big"-.5604719764011799410029498525073746312684365781710914454277286135693215339233038348083e-1")
    α[7] = convert(T,
        big".1789297658862876254180602006688963210702341137123745819397993311036789297658862876254")
    α[8] = convert(T, big"0.")
    α[9] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[10] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[11] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[12] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[13] = convert(T,
        big"-.1789297658862876254180602006688963210702341137123745819397993311036789297658862876254")
    α[14] = convert(T,
        big".5604719764011799410029498525073746312684365781710914454277286135693215339233038348083e-1")
    α[15] = convert(T,
        big".5671077504725897920604914933837429111531190926275992438563327032136105860113421550095e-1")
    α[16] = convert(T,
        big".2192242833052276559865092748735244519392917369308600337268128161888701517706576728499e-1")
    α[17] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 10, stability_size = 3.381557905977185))
end

"""
Feagin10 in Tableau form
"""
function constructFeagin10(T::Type = BigFloat)
    A = zeros(T, 17, 17)
    c = zeros(T, 17)
    α = zeros(T, 17)

    c[2] = convert(T, big".1")
    c[3] = convert(T,
        big".5393578408029817875324851978813024368572734497010090155054997959606637421764517472534")
    c[4] = convert(T,
        big".8090367612044726812987277968219536552859101745515135232582496939409956132646776208802")
    c[5] = convert(T,
        big".3090367612044726812987277968219536552859101745515135232582496939409956132646776208802")
    c[6] = convert(T,
        big".9810741902197952682548795483105620804890567461187248820278053092163063076419899423642")
    c[7] = convert(T,
        big".8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[8] = convert(T,
        big".3540173658568023763292641859487967421158240538073739683241838300209849970156463319588")
    c[9] = convert(T,
        big".8825276619647323464255014869796690751828678442680521196637911779185276585194132570617")
    c[10] = convert(T,
        big".6426157582403225481570754970204395359595017363632126959098752082638489654570997990908")
    c[11] = convert(T,
        big".3573842417596774518429245029795604640404982636367873040901247917361510345429002009092")
    c[12] = convert(T,
        big".1174723380352676535744985130203309248171321557319478803362088220814723414805867429383")
    c[13] = convert(T,
        big".8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[14] = convert(T,
        big".3090367612044726812987277968219536552859101745515135232582496939409956132646776208802")
    c[15] = convert(T,
        big".5393578408029817875324851978813024368572734497010090155054997959606637421764517472534")
    c[16] = convert(T, big".1")
    c[17] = convert(T, big"1.")
    A[2, 1] = convert(T, big".1")
    A[3, 1] = convert(T,
        big"-.9151765613752914405200150192753421543189513876643697205646603979963048584489560617612")
    A[3, 2] = convert(T,
        big"1.454534402178273228052500217156644591176224837365378736070160193956968600625407809015")
    A[4, 1] = convert(T,
        big".2022591903011181703246819492054884138214775436378783808145624234852489033161694052200")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T,
        big".6067775709033545109740458476164652414644326309136351424436872704557467099485082156601")
    A[5, 1] = convert(T,
        big".1840247147086435751491006934711206642167740479795914178446349661389134467299750502903")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T,
        big".1979668312271923690681417705103887933706372874633604015557458756877941908534751780490")
    A[5, 4] = convert(T,
        big"-.7295478473136326291851466715955580230150116089143829614213114788571202431877260745922e-1")
    A[6, 1] = convert(T,
        big".8790073402066813373197770941321254759188868249445485340413779789589432902346930000504e-1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T,
        big".4104597025202606453181748959204534260880353259028486952104060382836232634567776433876")
    A[6, 5] = convert(T,
        big".4827137536788664892047269429768961068091327377214213334132614730367887151617429989716")
    A[7, 1] = convert(T,
        big".8597005049024603021884802259458084014111326156366002225938797696846728899962000680512e-1")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".3308859630407221839488840576587531736482401548384020334486324456561774132656637353804")
    A[7, 5] = convert(T,
        big".4896629573094501928445070111358982011780154784337900972107904505844628620805256247950")
    A[7, 6] = convert(T,
        big"-.7318563750708507367890575805589888163403556150251881958547753987577423101247603364721e-1")
    A[8, 1] = convert(T,
        big".1209304491253337206603788549276689539589389969997036788126207242479388194626784582422")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big".2601246757582956228090076178383351743681087564846933618878385480736948525315211889536")
    A[8, 6] = convert(T,
        big".3254026215490913301588993343912312593327166759927000007761014763410936811599410398613e-1")
    A[8, 7] = convert(T,
        big"-.5957802118173610015601222025633051214449536727629307245388558993475804309454741922317e-1")
    A[9, 1] = convert(T,
        big".1108543795803914835089361710102184419094257801686565598070379712238628334864034048948")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big"-.6057614882550055876209249536555168755263444153543392346194660745470782044984741351417e-1")
    A[9, 7] = convert(T,
        big".3217637056017783901008987990498789040814043686030771292511102949020610017274377990772")
    A[9, 8] = convert(T,
        big".5104857256080630315777590122851234167446721370317523540675895192473116437554194666039")
    A[10, 1] = convert(T,
        big".1120544147528790048297150027618023630037176111581722293293926514934695608734670818246")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big"-.1449427759028659156723498283409807771816684997485068388761852225169215686873306825460")
    A[10, 7] = convert(T,
        big"-.3332697190962567065897052114157468717094674239921154979687242221418137519998082528658")
    A[10, 8] = convert(T,
        big".4992692295568800613533168439699785678602768165926732012403315488354044828766525571429")
    A[10, 9] = convert(T,
        big".5095046089296861042360986900453862539866432323529896021850604525937102423941190955351")
    A[11, 1] = convert(T,
        big".1139767839641859861380041867369011638907247525414868316403412210309290664271336257765")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T,
        big"-.7688133642033569385862142891208952708213490233909229874063835373263841258915727093324e-1")
    A[11, 7] = convert(T,
        big".2395273603243906491077114552718823730197413112010041193395628785781224443177793593046")
    A[11, 8] = convert(T,
        big".3977746623680946390478304624889521045647164163434546399026133008178160445613368390966")
    A[11, 9] = convert(T,
        big".1075589568736074555506091474414774502571367828232808385470239734312131460138378663157e-1")
    A[11, 10] = convert(T,
        big"-.3277691241640188741470610873502333953782629923923940719064566523011994227755761389669")
    A[12, 1] = convert(T,
        big".7983145282801960463514268644864003227587376304234139453562837151408388507380447344593e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T,
        big"-.5203296868006030765149498876129590687213114438816835269372978568212814716048923869522e-1")
    A[12, 7] = convert(T,
        big"-.5769541461685488817327843552834335090661592871529687230218639793413055668210844841955e-1")
    A[12, 8] = convert(T,
        big".1947819157121041649763062621473828711561429213544093647380902228866390247033946064390")
    A[12, 9] = convert(T,
        big".1453849231883250697275248259770711948592034675682365238665823980908913337790067458849")
    A[12, 10] = convert(T,
        big"-.7829427103516707775539867297256924472520770472391605513350159142016458207969788779419e-1")
    A[12, 11] = convert(T,
        big"-.1145032993610989121843031642905546709701332184056581226746743953737186161533235079226")
    A[13, 1] = convert(T,
        big".9851156101648572801200415003065172784136466773141955595205285811315908348605429600665")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T,
        big".3308859630407221839488840576587531736482401548384020334486324456561774132656637353804")
    A[13, 5] = convert(T,
        big".4896629573094501928445070111358982011780154784337900972107904505844628620805256247950")
    A[13, 6] = convert(T,
        big"-1.378964865748435675821127209307519023539043271485594715263967279333987427873524280873")
    A[13, 7] = convert(T,
        big"-.8611641950276356666739169996655345733510260609874270933144115020405879169494864941030")
    A[13, 8] = convert(T,
        big"5.784288136375372200229997854865784360068727896894991726018561768761366739727928637196")
    A[13, 9] = convert(T,
        big"3.288077619851035668904606159373148054772682529033423565819249490150659566267799708492")
    A[13, 10] = convert(T,
        big"-2.386339050931363840134223252155278661484014659759541045858065415400776902303795495458")
    A[13, 11] = convert(T,
        big"-3.254793424836439186545893675877887267477115046747806802699112121425479969838061889986")
    A[13, 12] = convert(T,
        big"-2.163435416864229823539542113000548208896780364201099991548873084750091865904259172177")
    A[14, 1] = convert(T,
        big".8950802957716328910496131323365851381481562792415613459917095170757672839845525640729")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T,
        big".1979668312271923690681417705103887933706372874633604015557458756877941908534751780490")
    A[14, 4] = convert(T,
        big"-.7295478473136326291851466715955580230150116089143829614213114788571202431877260745922e-1")
    A[14, 5] = convert(T, big"0.")
    A[14, 6] = convert(T,
        big"-.8512362396620076197390493714459667932893597228757022271661048887740645775128087411721")
    A[14, 7] = convert(T,
        big".3983201123185333017197186141743736433364809181037739042318564985414294618462922613810")
    A[14, 8] = convert(T,
        big"3.639372631810356060294129200470900441320273878939778041762287245257885919339182367288")
    A[14, 9] = convert(T,
        big"1.548228770398303223653016630751745649199817363489734963130650682036414979412951081310")
    A[14, 10] = convert(T,
        big"-2.122217147040537160260624274604272610253184611462601244015607790107056265271775849169")
    A[14, 11] = convert(T,
        big"-1.583503985453261727133843496257532127572691889344342379752907276826529676767705433797")
    A[14, 12] = convert(T,
        big"-1.715616082859362649220318197513490989126158808275519929730335787990138709296569362588")
    A[14, 13] = convert(T,
        big"-.2440364057501274521354154444122168754655935983709105660691323307479496900414383703505e-1")
    A[15, 1] = convert(T,
        big"-.9151765613752914405200150192753421543189513876643697205646603979963048584489560617612")
    A[15, 2] = convert(T,
        big"1.454534402178273228052500217156644591176224837365378736070160193956968600625407809015")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T,
        big"-.7773336436449682335389312285753021378033510536295472863344690943163799350772376765098")
    A[15, 6] = convert(T, big"0.")
    A[15, 7] = convert(T,
        big"-.9108956621551760695932035558074842001118890917701017996479849775333987275797724426330e-1")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T, big"0.")
    A[15, 10] = convert(T, big"0.")
    A[15, 11] = convert(T, big"0.")
    A[15, 12] = convert(T, big"0.")
    A[15, 13] = convert(T,
        big".9108956621551760695932035558074842001118890917701017996479849775333987275797724426330e-1")
    A[15, 14] = convert(T,
        big".7773336436449682335389312285753021378033510536295472863344690943163799350772376765098")
    A[16, 1] = convert(T, big".1")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T,
        big"-.1571786657997711633670589982731289218671837541267094194096536643481972116830948646637")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T, big"0.")
    A[16, 10] = convert(T, big"0.")
    A[16, 11] = convert(T, big"0.")
    A[16, 12] = convert(T, big"0.")
    A[16, 13] = convert(T, big"0.")
    A[16, 14] = convert(T, big"0.")
    A[16, 15] = convert(T,
        big".1571786657997711633670589982731289218671837541267094194096536643481972116830948646637")
    A[17, 1] = convert(T,
        big".1817813007000952838884720625822623796504438314631995216649446510802956214784300321820")
    A[17, 2] = convert(T, big".675")
    A[17, 3] = convert(T,
        big".3427581598471898399422205534138508717423387347039589199372599557781883490612271042488")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T,
        big".2591112145483227445129770761917673792677836845431824287781563647721266450257458921699")
    A[17, 6] = convert(T,
        big"-.3582789667179520890489612767219793977397506346732688024842708024524021503320635232761")
    A[17, 7] = convert(T,
        big"-1.045948959408833060950500687564099051315881231723784892860799571162295692535474761226")
    A[17, 8] = convert(T,
        big".9303278454156269832923005644324287771376016511829657946803974409790117863060139928499")
    A[17, 9] = convert(T,
        big"1.779509594317081024461421067948244539262757432433277905359997500386920301161540307119")
    A[17, 10] = convert(T, big".1")
    A[17, 11] = convert(T,
        big"-.2825475695390440816124777852222872764084893759762111899528770188407477720704209073106")
    A[17, 12] = convert(T,
        big"-.1593273501199725491692619843734858592780315421275519314618208481319295132597866290481")
    A[17, 13] = convert(T,
        big"-.1455158946470015108609919610810841113086501305786264049455713518588525807482385112893")
    A[17, 14] = convert(T,
        big"-.2591112145483227445129770761917673792677836845431824287781563647721266450257458921699")
    A[17, 15] = convert(T,
        big"-.3427581598471898399422205534138508717423387347039589199372599557781883490612271042488")
    A[17, 16] = convert(T, big"-.675")
    α[1] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    α[2] = convert(T, big".25e-1")
    α[3] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big".5e-1")
    α[6] = convert(T, big"0.")
    α[7] = convert(T, big".4e-1")
    α[8] = convert(T, big"0.")
    α[9] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[10] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[11] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[12] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[13] = convert(T, big"-.4e-1")
    α[14] = convert(T, big"-.5e-1")
    α[15] = convert(T,
        big"-.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    α[16] = convert(T, big"-.25e-1")
    α[17] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 10, stability_size = 2.527944696214181))
end

"""
A Runge-Kutta Method of Order 10, E. Hairer, J. Inst. Maths Applics (1978) 21, 47-59.
"""
function constructHairer10(T::Type = BigFloat)
    A = zeros(T, 17, 17)
    c = zeros(T, 17)
    α = zeros(T, 17)

    c[2] = convert(T,
        big".5233584004620047139632937023215170497515953383496781610502942213792083195573343614542")
    c[3] = convert(T,
        big".5265091001416125727329516734775408259523008085442588621240322448552569517961160776999")
    c[4] = convert(T,
        big".7897636502124188590994275102163112389284512128163882931860483672828854276941741165498")
    c[5] = convert(T,
        big".3939235701256720143227738119996521367488350094732736567444747389961242969755195414769")
    c[6] = convert(T,
        big".7666539862535505911932668693560686601417881683220066628197494388337786524595845606945")
    c[7] = convert(T,
        big".2897636502124188590994275102163112389284512128163882931860483672828854276941741165498")
    c[8] = convert(T,
        big".1084776892195672933536461100396272126536495082872530661076180009074872219675135232227")
    c[9] = convert(T,
        big".3573842417596774518429245029795604640404982636367873040901247917361510345429002009092")
    c[10] = convert(T,
        big".8825276619647323464255014869796690751828678442680521196637911779185276585194132570617")
    c[11] = convert(T,
        big".6426157582403225481570754970204395359595017363632126959098752082638489654570997990908")
    c[12] = convert(T,
        big".1174723380352676535744985130203309248171321557319478803362088220814723414805867429383")
    c[13] = convert(T,
        big".7666539862535505911932668693560686601417881683220066628197494388337786524595845606945")
    c[14] = convert(T,
        big".2897636502124188590994275102163112389284512128163882931860483672828854276941741165498")
    c[15] = convert(T,
        big".5265091001416125727329516734775408259523008085442588621240322448552569517961160776999")
    c[16] = convert(T,
        big".5233584004620047139632937023215170497515953383496781610502942213792083195573343614542")
    c[17] = convert(T, big"1.")
    A[2, 1] = convert(T,
        big".5233584004620047139632937023215170497515953383496781610502942213792083195573343614542")
    A[3, 1] = convert(T,
        big".2616697163778127283312402097548997641973627614039327706102102491953717704187699219879")
    A[3, 2] = convert(T,
        big".2648393837637998444017114637226410617549380471403260915138219956598851813773461557120")
    A[4, 1] = convert(T,
        big".1974409125531047147748568775540778097321128032040970732965120918207213569235435291375")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T,
        big".5923227376593141443245706326622334291963384096122912198895362754621640707706305874124")
    A[5, 1] = convert(T,
        big".1973205486287023067036649485978952117578825584337199991656132317825743833021058860612")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T,
        big".2950833340926721918228255598274359230145509784596048092593866299668586683096812321792")
    A[5, 4] = convert(T,
        big"-.9848031259570248420371669642567899802359852742005115168052512275330875463626757676351e-1")
    A[6, 1] = convert(T,
        big".1313134173444616536130177999345909470542768366990469474228776563495603985387834598888")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T,
        big".1101544395386396206773696967716892932905881833462590372574439152868807925960672597269")
    A[6, 5] = convert(T,
        big".5251861293704493169028793726497884197969231482767006781394278671973374613247338410788")
    A[7, 1] = convert(T,
        big".1342003418463226002727476951680931444018781918996634475042938727232827990267510411664")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".6960887032881160802299824047678314828416281106463280160258778625630871045892118300249")
    A[7, 5] = convert(T,
        big".2504977215703398097125518092509800218006823479445679226772611578145189247821334145350")
    A[7, 6] = convert(T,
        big"-.7910231164923596311158543989705934101157374376741710930213845258180034007039221691766")
    A[8, 1] = convert(T,
        big".7221827418966261942005081845517049770474117717860435538167259617193885317787239789517e-1")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big"-.5833632293645610716380606138930651874161115931850840194371530896953184153028603826158e-1")
    A[8, 6] = convert(T,
        big".3047557668574525220174950070294036092519082552287015783049982467857131629284936232933e-2")
    A[8, 7] = convert(T,
        big".9154818029778625587722640290346919759800040787487009688661073123722307869064222735619e-1")
    A[9, 1] = convert(T,
        big".3125500813516617947050120528263476612150928811800844295622424453129979711857118942140e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big".1091238215424128929294834955207716494619546474783251185719596191189550956073995218558e-3")
    A[9, 7] = convert(T,
        big".1567257586309938356246107465648794708917441337700138495832023584660674808897051732313")
    A[9, 8] = convert(T,
        big".1692943511719750238548830676365254553777828871012866864321262291196648014390164387346")
    A[10, 1] = convert(T,
        big".1190660441466861924216884258080925099509546061156163781761478570338134316043617322173e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big".2834370820246027860255992266982188665428702252125274101591962479489267620092644600260")
    A[10, 7] = convert(T,
        big"-.4163121675706282353724276181356300212164192944619403444080980777942870933794472685216")
    A[10, 8] = convert(T,
        big".2646463339497663668210902091085361027053564926326924812994390366777834273576276339310")
    A[10, 9] = convert(T,
        big".7388498091463228097090708267277348761559649602732109347956391853827232193715322584046")
    A[11, 1] = convert(T,
        big".2340657369133197891470838377984007842503946857756845416362339865999662377057916960290e-1")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T,
        big".9449313018949365401300253095605614324982516623777346096448800625130954018295939047740e-1")
    A[11, 7] = convert(T,
        big"-.2728720559019952606363092580665963250433705067252372208829562550632597611313198545757")
    A[11, 8] = convert(T,
        big".2240220461156057997944315522518131846261124699330640294272458923970695162204120486767")
    A[11, 9] = convert(T,
        big".6043814410751657569719347222576085340011863610739072982875214128849988434755301302413")
    A[11, 10] = convert(T,
        big"-.3081537692927938090069243415828207929929122273386332605004724686626579706106108533173e-1")
    A[12, 1] = convert(T,
        big".4544377531017616315765389908153096498645890941991961177836633137902353666760633687088e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T,
        big"-.1187996671864028586765254219285356343376285990176386474891150929245660355742130860183e-2")
    A[12, 7] = convert(T,
        big".1203565499092261097966188217234362058515446695047694116231895919296441480090125575596e-1")
    A[12, 8] = convert(T,
        big".7512690298764966821627521371565572140275315500656240513504528112598150181393445048666e-1")
    A[12, 9] = convert(T,
        big"-.1822092409888012403141186105974838892761579859192182074696828369088959767419077567178e-1")
    A[12, 10] = convert(T,
        big"-.2571528540841043468806376221771396205460354181513400383106090430345956162981031278207e-3")
    A[12, 11] = convert(T,
        big".4532078371347468185965270952011502734303744355238469520648294046672741844375709484540e-2")
    A[13, 1] = convert(T,
        big".1767137782592772030958798765711993346076326211800572275450227165783753236705910865492")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T,
        big".1101544395386396206773696967716892932905881833462590372574439152868807925960672597269")
    A[13, 5] = convert(T,
        big".5251861293704493169028793726497884197969231482767006781394278671973374613247338410788")
    A[13, 6] = convert(T,
        big"-.4716207672801957948798217912152359376250630852495511063738116933651587031904328351457")
    A[13, 7] = convert(T,
        big".8990310498491875266368990071875152922763468480002185650326986125011485318362907529907")
    A[13, 8] = convert(T,
        big"-.7467230306916289638599602008088168117750310724922743198498253813592425510843163068237")
    A[13, 9] = convert(T,
        big"-1.017101516756146040853186972006065972987027196800421553809421717321497529906933631477")
    A[13, 10] = convert(T,
        big".1263508715195988962951307827687648346421985369266969430473204298972536422365713122404")
    A[13, 11] = convert(T,
        big".5660138272355064270682732249907470012763799581315503842554078250210353407723389384909")
    A[13, 12] = convert(T,
        big".5986492052088624001098038724464832066388402270027708075754868643976463442046741430643")
    A[14, 1] = convert(T,
        big".1277534947480869822694777006880571541639616513225826576695303067404023367054772185702")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T,
        big".6960887032881160802299824047678314828416281106463280160258778625630871045892118300249")
    A[14, 5] = convert(T,
        big".2504977215703398097125518092509800218006823479445679226772611578145189247821334145350")
    A[14, 6] = convert(T,
        big"-.7368246436028416867609246757454535374296880219263938462439002090823944915566264811824")
    A[14, 7] = convert(T,
        big"-.2778578777108241826773273374900723250222301109862216853553157201018147214465526588169")
    A[14, 8] = convert(T,
        big"-.5997526313598403501296884799197753021563938240370770948150479630779446286262003432092")
    A[14, 9] = convert(T,
        big".2024692338910704693500237585621903123505161701229471467587157451308903694383321235511")
    A[14, 10] = convert(T,
        big".5432036982363849780600684652634443601468189969678666775046718813224989883416104871445e-2")
    A[14, 11] = convert(T,
        big"-.1074472474155047920101206919894381337125444664272205024314798936418733258920769563370e-1")
    A[14, 12] = convert(T,
        big".6951688484570234004700591858164146072357628221597117426434839740273190245052113679250")
    A[14, 13] = convert(T,
        big"-.6246651130952503394431547116755180508600167575701318270645551618021614799102076408562e-1")
    A[15, 1] = convert(T,
        big".2616697163778127283312402097548997641973627614039327706102102491953717704187699219879")
    A[15, 2] = convert(T,
        big".2648393837637998444017114637226410617549380471403260915138219956598851813773461557120")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T,
        big"-.1998011270205324791079663580830885049848273745422651189682301346802905866051733476638")
    A[15, 7] = convert(T,
        big"-.6510499873052827124921914489683813643155863882516440645794556633240216912803403931627")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T, big"0.")
    A[15, 10] = convert(T, big"0.")
    A[15, 11] = convert(T, big"0.")
    A[15, 12] = convert(T, big"0.")
    A[15, 13] = convert(T,
        big".1998011270205324791079663580830885049848273745422651189682301346802905866051733476638")
    A[15, 14] = convert(T,
        big".6510499873052827124921914489683813643155863882516440645794556633240216912803403931627")
    A[16, 1] = convert(T,
        big".5233584004620047139632937023215170497515953383496781610502942213792083195573343614542")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T,
        big"-.5558812136754302060726143105309293455559184141943321053532734480099926250948077261183")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T, big"0.")
    A[16, 10] = convert(T, big"0.")
    A[16, 11] = convert(T, big"0.")
    A[16, 12] = convert(T, big"0.")
    A[16, 13] = convert(T, big"0.")
    A[16, 14] = convert(T, big"0.")
    A[16, 15] = convert(T,
        big".5558812136754302060726143105309293455559184141943321053532734480099926250948077261183")
    A[17, 1] = convert(T,
        big".5732079543206559103114261705103983656495216504867462310285994428078568043160654439795e-1")
    A[17, 2] = convert(T,
        big"-.5499710763899945608115841896290187887481592249811405834035066676393750158953834290913")
    A[17, 3] = convert(T,
        big"-.6499374174008749135116607420010890619711618624173024222960650740195874521599402439688")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T,
        big"-1.061667370401756207240019539023157074172524666307437022389776456477183230723296269940")
    A[17, 7] = convert(T,
        big"-.4040156689806358294269682234212183308262562023912486365220642577870402491555711062480e-1")
    A[17, 8] = convert(T,
        big"-.1828302366407607254710272774065261039379052622607190097473388370699414811305446343873")
    A[17, 9] = convert(T,
        big"-.3336592706492786845666575661828162687906558601961826440714525336287466822150370633233")
    A[17, 10] = convert(T,
        big".3956485423760567568801345107166015519577734440834727480004748180136901286634710478955")
    A[17, 11] = convert(T,
        big".6950570494599735891002099282005158129027126868215679095299345058137097320818106877162")
    A[17, 12] = convert(T,
        big".2714873764573748588377263058539220945263829691804714618529052530298982146739754552950")
    A[17, 13] = convert(T,
        big".6071810560414041202873774349794680164722661545496003750296400378855628528787164400954")
    A[17, 14] = convert(T,
        big".5918636248229842840838104081530739675596239893196764223449596939309288102548549028752")
    A[17, 15] = convert(T,
        big".6499374174008749135116607420010890619711618624173024222960650740195874521599402439688")
    A[17, 16] = convert(T,
        big".5499710763899945608115841896290187887481592249811405834035066676393750158953834290913")
    α[1] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    α[2] = convert(T,
        big"-.3846153846153846153846153846153846153846153846153846153846153846153846153846153846154e-1")
    α[3] = convert(T,
        big"-.9090909090909090909090909090909090909090909090909090909090909090909090909090909090909e-1")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big"0.")
    α[6] = convert(T,
        big"-.1348314606741573033707865168539325842696629213483146067415730337078651685393258426966")
    α[7] = convert(T,
        big"-.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
    α[8] = convert(T, big"0.")
    α[9] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[10] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[11] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[12] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[13] = convert(T,
        big".1348314606741573033707865168539325842696629213483146067415730337078651685393258426966")
    α[14] = convert(T,
        big".1111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
    α[15] = convert(T,
        big".9090909090909090909090909090909090909090909090909090909090909090909090909090909090909e-1")
    α[16] = convert(T,
        big".3846153846153846153846153846153846153846153846153846153846153846153846153846153846154e-1")
    α[17] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 10, stability_size = 2.7046790686948667))
end

"""
High-order Explicit Runge-Kutta Formulae, Their uses, and Limitations, A.R.Curtis, J. Inst. Maths Applics (1975) 16, 35-55.
"""
function constructCurtis10(T::Type = BigFloat)
    A = zeros(T, 18, 18)
    c = zeros(T, 18)
    α = zeros(T, 18)
    c[2] = convert(T,
        big".1452518960316150517617548528770033320314511251329947060838468741983976455607179673401")
    c[3] = convert(T,
        big".1452518960316150517617548528770033320314511251329947060838468741983976455607179673401")
    c[4] = convert(T,
        big".2178778440474225776426322793155049980471766876994920591257703112975964683410769510101")
    c[5] = convert(T,
        big".5446946101185564441065806982887624951179417192487301478144257782439911708526923775252")
    c[6] = convert(T,
        big".6536335321422677329278968379465149941415300630984761773773109338927894050232308530303")
    c[7] = convert(T,
        big".2746594919905254008808021630247618520892150865127407293922085868737635475402543533498")
    c[8] = convert(T,
        big".7735775201106609448405825008093973718589542913426807556412662673054607938029043386501")
    c[9] = convert(T,
        big".5801831400829957086304368756070480288942157185070105667309497004790955953521782539876")
    c[10] = convert(T,
        big".1174723380352676535744985130203309248171321557319478803362088220814723414805867429383")
    c[11] = convert(T,
        big".3573842417596774518429245029795604640404982636367873040901247917361510345429002009092")
    c[12] = convert(T,
        big".6426157582403225481570754970204395359595017363632126959098752082638489654570997990908")
    c[13] = convert(T,
        big".1174723380352676535744985130203309248171321557319478803362088220814723414805867429383")
    c[14] = convert(T,
        big".8825276619647323464255014869796690751828678442680521196637911779185276585194132570617")
    c[15] = convert(T,
        big".3573842417596774518429245029795604640404982636367873040901247917361510345429002009092")
    c[16] = convert(T,
        big".6426157582403225481570754970204395359595017363632126959098752082638489654570997990908")
    c[17] = convert(T,
        big".8825276619647323464255014869796690751828678442680521196637911779185276585194132570617")
    c[18] = convert(T, big"1.")
    A[2, 1] = convert(T,
        big".1452518960316150517617548528770033320314511251329947060838468741983976455607179673401")
    A[3, 1] = convert(T,
        big".7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1")
    A[3, 2] = convert(T,
        big".7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1")
    A[4, 1] = convert(T,
        big".5446946101185564441065806982887624951179417192487301478144257782439911708526923775252e-1")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T,
        big".1634083830355669332319742094866287485353825157746190443443277334731973512558077132576")
    A[5, 1] = convert(T,
        big".5446946101185564441065806982887624951179417192487301478144257782439911708526923775252")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T,
        big"-2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720")
    A[5, 4] = convert(T,
        big"2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720")
    A[6, 1] = convert(T,
        big".6536335321422677329278968379465149941415300630984761773773109338927894050232308530303e-1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T,
        big".3268167660711338664639484189732574970707650315492380886886554669463947025116154265151")
    A[6, 5] = convert(T,
        big".2614534128569070931711587351786059976566120252393904709509243735571157620092923412121")
    A[7, 1] = convert(T,
        big".8233707757482716585173454344310125296066814318521742241762319051772963627695955263034e-1")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".2119171963202803561687843468555305553175658807629274312902985594840086570224567152664")
    A[7, 5] = convert(T,
        big"-.3997343508054218311577932550061320162379840049816347807630118786107674477850206579628e-1")
    A[7, 6] = convert(T,
        big".2037865317596006197606259822674324543477946306275935376058802473310199901934015124941e-1")
    A[8, 1] = convert(T,
        big".8595305779007343831562027786771081909543936570474230618236291858949564375587825985001e-1")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T, big"0.")
    A[8, 6] = convert(T,
        big".2911769478058850960337179621761553399856026049598393013981874594942289837064329700000")
    A[8, 7] = convert(T,
        big".3964475145147024104912442607655312127779123206780991480607158892217361663405931088001")
    A[9, 1] = convert(T,
        big".8612093485606967549983047372292119178898514571588438099912534616486575243508895957628e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big".1397464826824442089036313891001189801074425314582326737716288563521183595455090268480")
    A[9, 7] = convert(T,
        big".3951098495815674599900526056001284215294125840404176924334653987770478924197803010468")
    A[9, 8] = convert(T,
        big"-.4079412703708563576307759281612056453162454270752418047326990081493640904820003348350e-1")
    A[10, 1] = convert(T,
        big".7233144422337948077616348229119326315582930871089020733092900891206129381937795204778e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big".2200276284689998102140972735735070061373242800181187459951219347361114857342828430157")
    A[10, 7] = convert(T,
        big".8789533425436734013369780264792573637952226487753296416823846876217040795688489371334e-1")
    A[10, 8] = convert(T,
        big"-.4445383996260350863990674880611108986832860648196030000580004690002268108984238641730e-1")
    A[10, 9] = convert(T,
        big"-.2183282289488754689095532966861839909872150913926337371522805434288481649401165594213")
    A[11, 1] = convert(T,
        big".8947100936731114228785441966773836169071038390882857211057269158522704971585365845223e-1")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T,
        big".3946008170285561860741397654755022300929434262701385530048127140223687993778661654316")
    A[11, 7] = convert(T,
        big".3443011367963333487713764986067104675654371857504670290688086760696354596195596354011")
    A[11, 8] = convert(T,
        big"-.7946682664292661290694938113119430997053815140863772328764150866582492425892231395780e-1")
    A[11, 9] = convert(T,
        big"-.3915218947895966123834967996391962853380545808840091268064277812752553499114569444180")
    A[11, 10] = convert(T, big"0.")
    A[12, 1] = convert(T,
        big".3210006877963209212945282736072241886741425314298532400216927262619488479186214523312e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T, big"0.")
    A[12, 7] = convert(T, big"0.")
    A[12, 8] = convert(T,
        big"-.1846375997512050141835163881753227910996323204749769226655464078048769505209525299752e-3")
    A[12, 9] = convert(T,
        big".1560894025313219860759149162557283383430181475726228517203663063649626288079337909898")
    A[12, 10] = convert(T,
        big".1934496857654560252749984220385188727138526287670744309970093278715606577140084022992")
    A[12, 11] = convert(T,
        big".2611612387636636496908928477536452288263163392010050661129958478089356710938164130987")
    A[13, 1] = convert(T,
        big".4423749328524996327035388417792688154433173133294892285295756457561276315648477233732e-1")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T, big"0.")
    A[13, 5] = convert(T, big"0.")
    A[13, 6] = convert(T, big"0.")
    A[13, 7] = convert(T, big"0.")
    A[13, 8] = convert(T,
        big".4640774434539039636406222168781981616534115643208114455689698789119941732444857047798e-2")
    A[13, 9] = convert(T,
        big".4704660282615136532130927218172390570903230981414159347904277946537920001824903276586e-1")
    A[13, 10] = convert(T,
        big".8620749948011488160369445167416002799205317397013619044391270706339561700281526529703e-1")
    A[13, 11] = convert(T,
        big"-.2607983024682138093233254079066687623148682426317395111719299641390118652802949600035e-1")
    A[13, 12] = convert(T,
        big"-.3858020174396621532493277639159499581333235076531298977820093139813399390137768850940e-1")
    A[14, 1] = convert(T,
        big".2318046717429411567006043539613275607940758021709332569729352990777336390158311630529e-1")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T, big"0.")
    A[14, 5] = convert(T, big"0.")
    A[14, 6] = convert(T, big"0.")
    A[14, 7] = convert(T, big"0.")
    A[14, 8] = convert(T,
        big".3197856784116367067302124322582100058864027838197120089129330601737324659881765852593")
    A[14, 9] = convert(T,
        big".5933233331841898686063939886797828376866051205773280426848164018120869674204443797948")
    A[14, 10] = convert(T,
        big"-.1162415082107004073646376038875773625236111451677368940868844877968414717849631693075")
    A[14, 11] = convert(T,
        big".1803950557030502357344063195737827904476240180662764468232042537858892203518134072359")
    A[14, 12] = convert(T,
        big"-.4554014298857220726863505256926549022316460712353658688873150702827663762861750674926")
    A[14, 13] = convert(T,
        big".3374860655879838997354164406519929498380855579907450585197434903186534889285340052666")
    A[15, 1] = convert(T,
        big".2624364325798105891527733985858552391723553030719144065844544880498188553839263944447e-1")
    A[15, 2] = convert(T, big"0.")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T, big"0.")
    A[15, 7] = convert(T, big"0.")
    A[15, 8] = convert(T,
        big".4863139423867266106526843913609225996253073727381961544415263239431571586043622332760e-1")
    A[15, 9] = convert(T,
        big".4274382538346478867636942429421724367591866585774144180215122660980822123988151132213e-1")
    A[15, 10] = convert(T,
        big"-.1565711812511167545605069814291800290911012856968190479166777907500203971025172981530")
    A[15, 11] = convert(T,
        big".1326047194917652331781527125743684254490968718259563958293167893998110899691451568372")
    A[15, 12] = convert(T,
        big"-.9402962152946515651634831658142934852383791641671387741034606371378082209616938685225e-1")
    A[15, 13] = convert(T,
        big".3697316622986642308496397344700288190173622449975653151703435941016557663496526141090")
    A[15, 14] = convert(T,
        big"-.1197020013028860976492784934312243036670658451195397948726104511062042521592125912599e-1")
    A[16, 1] = convert(T,
        big".5568066641536216461090823068917803436066365804361903532125349474551476120813558125830e-1")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T, big"0.")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T,
        big"-.4324853319508358432896036654421685136736530810118924113940744870078036705505610668088")
    A[16, 9] = convert(T,
        big"-.9979726994172038714656907882931844552238093285811791155499130927685987422432191170216")
    A[16, 10] = convert(T,
        big".1129287214043638339947432805172320041719492778714292355381358939222504897262200136889")
    A[16, 11] = convert(T,
        big"-1.024823023512132929313567156576969954855232272749038347671818195935585095295127839150")
    A[16, 12] = convert(T,
        big"1.334565206642246959252239602313589265188981560552694580059808406200559397799055652161")
    A[16, 13] = convert(T,
        big".7216035483871342125753076728585010714020628779364042084073624629288704870682937447565e-2")
    A[16, 14] = convert(T,
        big".8992773696348355846430438306111181223414632598285854300924423251352733205187087732678e-1")
    A[16, 15] = convert(T,
        big"1.497578446211167333777988534023066333042434967475357134513165331964695787890042760189")
    A[17, 1] = convert(T,
        big"-.8434891199686377639125188391985671318383858641413517143104162188088468627447515172982e-3")
    A[17, 2] = convert(T, big"0.")
    A[17, 3] = convert(T, big"0.")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T, big"0.")
    A[17, 7] = convert(T, big"0.")
    A[17, 8] = convert(T,
        big".7602144218856081893754106886111596435015500427480120290148318740899211421773423234728")
    A[17, 9] = convert(T,
        big"1.769083927820959377467464871522349066447068428702073590698445112684989184432409492025")
    A[17, 10] = convert(T,
        big"-.3420660144261856568640856557101550842129621862743150134804423843180623377713466374754e-1")
    A[17, 11] = convert(T,
        big"1.490558190212043468817221563278239942209691100326719140478588601720867838040211450448")
    A[17, 12] = convert(T,
        big"-2.552203480132132516997563217309689292804518121743365818482497611667126218719069737195")
    A[17, 13] = convert(T,
        big".3301343553488974584509065658432587594272577682203521183331319470982274757841566234948")
    A[17, 14] = convert(T,
        big"-.9161854401769482236671414092387917470686251714192236693920061138984202381209109248553e-1")
    A[17, 15] = convert(T,
        big"-1.525735678746850818217653470352135651821164556169070505816135230784807058389577753184")
    A[17, 16] = convert(T,
        big".7371445601564892133467497107205798584829803038168267854389817508169123996459113657504")
    A[18, 1] = convert(T,
        big".1017366974111576638766809656369828971944080018220332809259398740674738807023371082700")
    A[18, 2] = convert(T, big"0.")
    A[18, 3] = convert(T, big"0.")
    A[18, 4] = convert(T, big"0.")
    A[18, 5] = convert(T, big"0.")
    A[18, 6] = convert(T, big"0.")
    A[18, 7] = convert(T, big"0.")
    A[18, 8] = convert(T,
        big"-1.696217553209432810711666838709742166182992092906177246174096517233561845662947862824")
    A[18, 9] = convert(T,
        big"-3.825235846211624254528740857512255693551264719132875740261231165548583482101116676418")
    A[18, 10] = convert(T,
        big"-.3420660144261856568640856557101550842129621862743150134804423843180623377713466374754e-1")
    A[18, 11] = convert(T,
        big"-2.520767789227152291196336314591227486393143379933686189126240710041836742414125694941")
    A[18, 12] = convert(T,
        big"5.472417145227780046950992000565734793413395536531652419585004300790370984185945495978")
    A[18, 13] = convert(T,
        big".7877467869749093540456207796605763462528869036641654587754177400832643963961769634278e-2")
    A[18, 14] = convert(T,
        big".3189152692455334369024560213486753019540464785641163242047782111839399471147176681561")
    A[18, 15] = convert(T,
        big"3.447227036527756718156475010324322155277035924051392880570525223655410460762027138915")
    A[18, 16] = convert(T,
        big"-.6051983612219277832241707671295607127814820499715293613761402732652780120810041653591")
    A[18, 17] = convert(T,
        big".3334525350307787459202631378414806560287636505658634784117511174230383993073398823363")
    α[1] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    α[2] = convert(T, big"0.")
    α[3] = convert(T, big"0.")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big"0.")
    α[6] = convert(T, big"0.")
    α[7] = convert(T, big"0.")
    α[8] = convert(T, big"0.")
    α[9] = convert(T, big"0.")
    α[10] = convert(T, big"0.")
    α[11] = convert(T, big"0.")
    α[12] = convert(T,
        big".1387145942588715882541801312803271702142521598590204181697361204933422401935856968980")
    α[13] = convert(T,
        big".1892374781489234901583064041060123262381623469486258303271944256799821862794952728707")
    α[14] = convert(T,
        big".9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1")
    α[15] = convert(T,
        big".2774291885177431765083602625606543404285043197180408363394722409866844803871713937960")
    α[16] = convert(T,
        big".1387145942588715882541801312803271702142521598590204181697361204933422401935856968980")
    α[17] = convert(T,
        big".9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1")
    α[18] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 10, stability_size = 3.826924804386427))
end

"""
On the 25 stage 12th order explicit Runge-Kutta method, by Hiroshi Ono.
Transactions of the Japan Society for Industrial and applied Mathematics, Vol. 6, No. 3, (2006) pages 177 to 186
"""
function constructOno12(T::Type = BigFloat)
    A = zeros(T, 25, 25)
    c = zeros(T, 25)
    α = zeros(T, 25)

    c[2] = convert(T, big".25")
    c[3] = convert(T,
        big".4444444444444444444444444444444444444444444444444444444444444444444444444444444444444")
    c[4] = convert(T,
        big".6666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    c[5] = convert(T,
        big".1083333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[6] = convert(T,
        big".1666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    c[7] = convert(T,
        big".5851851851851851851851851851851851851851851851851851851851851851851851851851851851852")
    c[8] = convert(T,
        big".6752380952380952380952380952380952380952380952380952380952380952380952380952380952381e-1")
    c[9] = convert(T, big".2")
    c[10] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[11] = convert(T,
        big".9446116054065563496847375720237340591364440847545279548637407620203329280139199114088")
    c[12] = convert(T,
        big".5179584680428461035835615460185207449756996828425551866437683588600807071775316523299e-1")
    c[13] = convert(T,
        big".8488805186071653506398389301626743020641481756400195420459339398355773991365476236893e-1")
    c[14] = convert(T,
        big".2655756032646428930981140590456168352972012641640776214486652703185222349414361456016")
    c[15] = convert(T, big".5")
    c[16] = convert(T,
        big".7344243967353571069018859409543831647027987358359223785513347296814777650585638543984")
    c[17] = convert(T,
        big".9151119481392834649360161069837325697935851824359980457954066060164422600863452376311")
    c[18] = convert(T,
        big".9446116054065563496847375720237340591364440847545279548637407620203329280139199114088")
    c[19] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[20] = convert(T, big".2")
    c[21] = convert(T,
        big".5851851851851851851851851851851851851851851851851851851851851851851851851851851851852")
    c[22] = convert(T,
        big".1666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    c[23] = convert(T,
        big".4444444444444444444444444444444444444444444444444444444444444444444444444444444444444")
    c[24] = convert(T, big".25")
    c[25] = convert(T, big"1.")
    A[2, 1] = convert(T, big".25")
    A[3, 1] = convert(T,
        big".4938271604938271604938271604938271604938271604938271604938271604938271604938271604938e-1")
    A[3, 2] = convert(T,
        big".3950617283950617283950617283950617283950617283950617283950617283950617283950617283951")
    A[4, 1] = convert(T,
        big".1666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T, big".5")
    A[5, 1] = convert(T,
        big".8775846354166666666666666666666666666666666666666666666666666666666666666666666666667e-1")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T, big".35318359375e-1")
    A[5, 4] = convert(T,
        big"-.1474348958333333333333333333333333333333333333333333333333333333333333333333333333333e-1")
    A[6, 1] = convert(T,
        big".3899572649572649572649572649572649572649572649572649572649572649572649572649572649573e-1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T,
        big".1036484245439469320066334991708126036484245439469320066334991708126036484245439469320e-3")
    A[6, 5] = convert(T,
        big".1275672917463962240081643066717693583365225156269932389335374410001275672917463962240")
    A[7, 1] = convert(T,
        big".5105992363578097460127638454113076883995951211314722974530930635183035731732576725718")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".7543245465503858040489014562406177717822511685842494078406423402976986226718361654108e-1")
    A[7, 5] = convert(T,
        big"-2.075547971680477873791142836860749355886837204528583027037257940355620678287181163745")
    A[7, 6] = convert(T,
        big"2.074701465852814732558674031010565075494202151723870973985285827992732428031925059817")
    A[8, 1] = convert(T,
        big".3676892528603925219048664448524211597614557351797636700024418222004670689026890820447e-1")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big".4990421592006925486260552960107757259429499555984070931321874496846090397321073859571e-1")
    A[8, 6] = convert(T,
        big"-.1930463579861287136353203957778747781438400759996112964824379316128697839120856079129e-1")
    A[8, 7] = convert(T,
        big".1553041163138881199636750152773130534672480459535771443046754965888913372527235149185e-3")
    A[9, 1] = convert(T,
        big".2067486743675349484922604488404063487529235328774704968666869007873453428790773240971e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big".7957154880008160975236540766621611282548264517609854377598122975695595623677029405014e-1")
    A[9, 7] = convert(T,
        big".4507214087766488746701444843200598485120376014363895694204088380480167395865728052023e-5")
    A[9, 8] = convert(T,
        big".9974907654907712890966184600490005170073988116014004264165587607592902930792610781210e-1")
    A[10, 1] = convert(T,
        big"-.3322925324853517529956086668161082115483929413336137015414300233353633578589390840556e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big"-.6436081333941628246797714549110840967717864507030812404304536614347636764815617688536")
    A[10, 7] = convert(T,
        big".1537717626809222300334752193656044639866567030541528516754813628383191013261883003503e-2")
    A[10, 8] = convert(T,
        big".3004530256573834493677742186357887570531778112266912246408744303615696525403823773529")
    A[10, 9] = convert(T,
        big".7081799766918386616445566840965834495669146999125431907603007531116805020471447502361")
    A[11, 1] = convert(T,
        big"5.940512035627220006057067930851135219656791728697945883200194895040343032390095577622")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T, big"0.")
    A[11, 7] = convert(T,
        big"1.962055480391860623022348248133501763323649793356240465634479033599078618610316304517")
    A[11, 8] = convert(T,
        big"-13.64213386305695955284608890764824855130586422219060958260554419417201101166023190424")
    A[11, 9] = convert(T,
        big"15.06683044587338301779467734916899621382077320594120734305207687034839214126066753426")
    A[11, 10] = convert(T,
        big"-8.382652493428947744343267048481650586358906421050256154417465842795469852586927600745")
    A[12, 1] = convert(T,
        big".2644349979861653392555978654560316035400311240360155334959299612007139374157863919974e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T, big"0.")
    A[12, 7] = convert(T, big"0.")
    A[12, 8] = convert(T,
        big".2907477949164468491495277367249908074184860785029007469445512693351379740267954217878e-1")
    A[12, 9] = convert(T,
        big"-.4670146655526017247312261556616256909459346454842633293588922375484653326459852610286e-2")
    A[12, 10] = convert(T,
        big".9537840586679400528486700806869537239400914426861191427349413827712021873350961548845e-3")
    A[12, 11] = convert(T,
        big"-.6069889118531287692814140320863412762496957479595228817306174863669287380259690136170e-5")
    A[13, 1] = convert(T,
        big".1731635805893703684448829089734942070445562382116745420114094326722705203174775140332e-1")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T, big"0.")
    A[13, 5] = convert(T, big"0.")
    A[13, 6] = convert(T, big"0.")
    A[13, 7] = convert(T, big"0.")
    A[13, 8] = convert(T, big"0.")
    A[13, 9] = convert(T,
        big".8690027159880925811795093256330430220184719404902985010762237927304116502630885425916e-3")
    A[13, 10] = convert(T,
        big"-.9198044746158460357290146202828770811642343120551843217420771155968903404182588708338e-4")
    A[13, 11] = convert(T,
        big".1833594777349928706130754092439075395240645268074916590065144720493485758190652981058e-6")
    A[13, 12] = convert(T,
        big".6679448817377525524901838117990401028051762116902291244289142812068791591710992924480e-1")
    A[14, 1] = convert(T,
        big".1497702285787817250603134142375551527881374815295532946701659971640303704469333256722e-1")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T, big"0.")
    A[14, 5] = convert(T, big"0.")
    A[14, 6] = convert(T, big"0.")
    A[14, 7] = convert(T, big"0.")
    A[14, 8] = convert(T, big"0.")
    A[14, 9] = convert(T,
        big".1299053198125687130806056561360945267563487713700848629960379786258472215118648828622")
    A[14, 10] = convert(T,
        big".4444252090193421415880477613428889941728587016143899156135432416265141679149665438354e-2")
    A[14, 11] = convert(T,
        big"-.2185380330434085597487818850115824768660155325792369439283645364044817815712239904765e-5")
    A[14, 12] = convert(T,
        big".6235770138025566019351747531432668514647398310437277172705875380249580436468314885283e-1")
    A[14, 13] = convert(T,
        big".5389349250407735998767659637686133399860483467584655047185578940287507515886082812093e-1")
    A[15, 1] = convert(T,
        big"1.188479339801401774163727997421987897149748610296498687871610457754515001479951581977")
    A[15, 2] = convert(T, big"0.")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T, big"0.")
    A[15, 7] = convert(T, big"0.")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T,
        big"-2.153577860743923948574230282195458010685777028877789457391059928453625236658492400496")
    A[15, 10] = convert(T,
        big".4084894136568439103913845855092578246970391902744397258692599626546701227557648623350")
    A[15, 11] = convert(T,
        big".7865354687677762592501346934798315168840512712394797311164536002596101630809266673394e-3")
    A[15, 12] = convert(T,
        big"-5.831493996616730192196764781137026841656419605282914766077967205018020861054852689274")
    A[15, 13] = convert(T,
        big"6.133497357760298833907642604705486519224977667578896200018420236452356862070881047498")
    A[15, 14] = convert(T,
        big".7538192106733418460489897410022727797535471147396301299786200230098445012436666712925")
    A[16, 1] = convert(T,
        big"-1.155141455453753782409702338155440396898375536349635780553293168727187831583039150444")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T, big"0.")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T,
        big"5.230035453923383710410708524533553450562429691771028687028452868714084651714450832431")
    A[16, 10] = convert(T,
        big"-.6054719632283663895583349987822330203147735874309119730944517939675864870705141356280")
    A[16, 11] = convert(T,
        big".9098600291918418987613542477924427728076737025978819058456352939211248672787562673474e-2")
    A[16, 12] = convert(T,
        big"7.389890569427277711001889579558669404987538661035016661449808071907466968267996694225")
    A[16, 13] = convert(T,
        big"-8.586861200440706242153541145492405736845705766762479792176742622565505838929692855382")
    A[16, 14] = convert(T,
        big"-2.289142858952360074158926453709367741105307029478199282759359773927787664228731981389")
    A[16, 15] = convert(T,
        big".7420172511679637547821792305236827765889155660251250395984647953087827182153068879117")
    A[17, 1] = convert(T,
        big"-.3878704284298115050812248647831359622498421729413088096991246117344300538467378525660")
    A[17, 2] = convert(T, big"0.")
    A[17, 3] = convert(T, big"0.")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T, big"0.")
    A[17, 7] = convert(T, big"0.")
    A[17, 8] = convert(T, big"0.")
    A[17, 9] = convert(T,
        big"-52.79686296396567660865119351124975527608768115963057325690048513861113171732667960572")
    A[17, 10] = convert(T,
        big"3.185273134611601001666377466850930616911166170460075757653190749060280689442040942426")
    A[17, 11] = convert(T,
        big"-.1777687036844069342486311680075793729237679916897779073551385167909574116891757259578")
    A[17, 12] = convert(T,
        big"-24.01291564936409810460023986270697223960990857367094000364464827510022831834721791049")
    A[17, 13] = convert(T,
        big"44.59765117143356991664568418613870421297485805339864304716728296273494210914975546336")
    A[17, 14] = convert(T,
        big"34.64254823700930167755101801088405791740501138964826361776800135713319689123911661552")
    A[17, 15] = convert(T,
        big"-5.881644810570531178309390232014921075961483123099551413638107203551765688978516422159")
    A[17, 16] = convert(T,
        big"1.746701961099335199963616081872403749335232589961167014444435282876535760443759733215")
    A[18, 1] = convert(T,
        big"7.233966671341530645490952497192110761349200861822849620574804240945984259787128547890")
    A[18, 2] = convert(T, big"0.")
    A[18, 3] = convert(T, big"0.")
    A[18, 4] = convert(T, big"0.")
    A[18, 5] = convert(T, big"0.")
    A[18, 6] = convert(T, big"0.")
    A[18, 7] = convert(T,
        big"1.962055480391860623022348248133501763323649793356240465634479033599078618610316304517")
    A[18, 8] = convert(T,
        big"-13.64213386305695955284608890764824855130586422219060958260554419417201101166023190424")
    A[18, 9] = convert(T,
        big"43.38367605377242869214585350579838046755997844867070816778334123860217800450037322178")
    A[18, 10] = convert(T,
        big"-9.834713911815061145761755868590694958877038434912171502071857284257585145054082454689")
    A[18, 11] = convert(T,
        big".1029076524695381286377868412846724978000944487949691840993018167527636501337920551920")
    A[18, 12] = convert(T,
        big"8.468423106151264851550086751522745687779805425203659086774340333971055043511878909430")
    A[18, 13] = convert(T,
        big"-20.03922776536488408316747608806415623443728880893196401395082287814408336477871062892")
    A[18, 14] = convert(T,
        big"-19.26671999607368316162814844363603221030906811653440141113230651624038103492726683046")
    A[18, 15] = convert(T,
        big"3.391598559113951777436829245234173685429846516750495577069567102933678571370358241684")
    A[18, 16] = convert(T,
        big"-.8461837296739539683436767973109554382501241857205608248520628152174372445350689400749")
    A[18, 17] = convert(T,
        big".3096334815052354314802658810823658907325235844531318754050068324709258105543338930304e-1")
    A[19, 1] = convert(T,
        big"-.8082042305712324316381428055181002921423821181387879506081394201551828774972588869862")
    A[19, 2] = convert(T, big"0.")
    A[19, 3] = convert(T, big"0.")
    A[19, 4] = convert(T, big"0.")
    A[19, 5] = convert(T, big"0.")
    A[19, 6] = convert(T,
        big"-.6436081333941628246797714549110840967717864507030812404304536614347636764815617688536")
    A[19, 7] = convert(T,
        big".1537717626809222300334752193656044639866567030541528516754813628383191013261883003503e-2")
    A[19, 8] = convert(T,
        big".3004530256573834493677742186357887570531778112266912246408744303615696525403823773529")
    A[19, 9] = convert(T,
        big"1.629527090754742423152222936373756585439149301882472540854080995266699866827769539264")
    A[19, 10] = convert(T,
        big"-.2429977955477512310051787601585012098624988470344665355076878808184221276200961070016")
    A[19, 11] = convert(T,
        big".1937321543128929118688572792462119078051789173168830514027044834102257771758364823616e-1")
    A[19, 12] = convert(T,
        big"3.652148675291515690405068402145314798647603571421741897364878082465880262359243365600")
    A[19, 13] = convert(T,
        big"-3.575712494072609116114226433534392208450672236215592934178802498800227318056660241871")
    A[19, 14] = convert(T, big"0.")
    A[19, 15] = convert(T,
        big".1424995469019583371886649267486631800285185459043846216158470864873378043488691875332e-1")
    A[19, 16] = convert(T,
        big".7783761260627937338705544952370637303437265442198231246348540202794422441986844495010e-2")
    A[19, 17] = convert(T,
        big".1680830476266175444138000910160343423421707378563384687234091879582812751321589095322e-4")
    A[19, 18] = convert(T,
        big"-.2123426209823757245364666745406479474016549497429582971324756529195024847371745454980e-1")
    A[20, 1] = convert(T,
        big".3030500118630462312535985041299218980792680301551577592229328311072204102313787331391")
    A[20, 2] = convert(T, big"0.")
    A[20, 3] = convert(T, big"0.")
    A[20, 4] = convert(T, big"0.")
    A[20, 5] = convert(T, big"0.")
    A[20, 6] = convert(T,
        big".7957154880008160975236540766621611282548264517609854377598122975695595623677029405014e-1")
    A[20, 7] = convert(T,
        big".4507214087766488746701444843200598485120376014363895694204088380480167395865728052023e-5")
    A[20, 8] = convert(T,
        big".9974907654907712890966184600490005170073988116014004264165587607592902930792610781210e-1")
    A[20, 9] = convert(T,
        big"-.3171617769987166876020401487839312465336733321856991443597739463291092021830833497646")
    A[20, 10] = convert(T,
        big".4097746751846734839576474906854331980839883999186796914837456369263814852323254415473")
    A[20, 11] = convert(T,
        big"-.1037579600906669233865265654007540316189053766128282845577718841615974679757246082771e-1")
    A[20, 12] = convert(T,
        big"-1.317367560379299247049913838592471137800861178079091445393764303876611018306956677528")
    A[20, 13] = convert(T,
        big"1.279912157519227083039575094951746274090584680488104164669880623922417279542229865347")
    A[20, 14] = convert(T, big"0.")
    A[20, 15] = convert(T, big"0.")
    A[20, 16] = convert(T,
        big"-.3909856506984242164691157929655038995858933888334921275260056794557284704385445582377e-2")
    A[20, 17] = convert(T,
        big"-.9117410119332842109125953032188840693633119317454401715027062207348542639139803659651e-5")
    A[20, 18] = convert(T,
        big".1121029070910995545199917388600394550626832254247437200910915070169440797714588293937e-1")
    A[20, 19] = convert(T,
        big"-.3344481605351170568561872909698996655518394648829431438127090301003344481605351170569")
    A[21, 1] = convert(T,
        big".5105992363578097460127638454113076883995951211314722974530930635183035731732576725718")
    A[21, 2] = convert(T, big"0.")
    A[21, 3] = convert(T, big"0.")
    A[21, 4] = convert(T,
        big".7543245465503858040489014562406177717822511685842494078406423402976986226718361654108e-1")
    A[21, 5] = convert(T,
        big"-2.075547971680477873791142836860749355886837204528583027037257940355620678287181163745")
    A[21, 6] = convert(T,
        big"2.074701465852814732558674031010565075494202151723870973985285827992732428031925059817")
    A[21, 7] = convert(T, big"0.")
    A[21, 8] = convert(T, big"0.")
    A[21, 9] = convert(T,
        big"-.6056348202998013364637890554451869243272530251038468484739028354704713743904641502619")
    A[21, 10] = convert(T,
        big"-.6962774898931662522054431030410802850752281852155429905113849488059348994873647212459")
    A[21, 11] = convert(T,
        big"-.2530464169447410633467649077644850846834331528772883593493989655560821683241998040578e-1")
    A[21, 12] = convert(T, big"0.")
    A[21, 13] = convert(T, big"0.")
    A[21, 14] = convert(T, big"0.")
    A[21, 15] = convert(T, big"0.")
    A[21, 16] = convert(T, big"0.")
    A[21, 17] = convert(T, big"0.")
    A[21, 18] = convert(T,
        big".2530464169447410633467649077644850846834331528772883593493989655560821683241998040578e-1")
    A[21, 19] = convert(T,
        big".6962774898931662522054431030410802850752281852155429905113849488059348994873647212459")
    A[21, 20] = convert(T,
        big".6056348202998013364637890554451869243272530251038468484739028354704713743904641502619")
    A[22, 1] = convert(T,
        big".3899572649572649572649572649572649572649572649572649572649572649572649572649572649573e-1")
    A[22, 2] = convert(T, big"0.")
    A[22, 3] = convert(T, big"0.")
    A[22, 4] = convert(T,
        big".1036484245439469320066334991708126036484245439469320066334991708126036484245439469320e-3")
    A[22, 5] = convert(T,
        big".1275672917463962240081643066717693583365225156269932389335374410001275672917463962240")
    A[22, 6] = convert(T, big"0.")
    A[22, 7] = convert(T,
        big".1146176189444805429210887300555916540321998411667027651433109522778138762544220633889")
    A[22, 8] = convert(T, big"0.")
    A[22, 9] = convert(T,
        big"-.1316367023754075454122030740568234746157428970656730321378667908709827666511411271542")
    A[22, 10] = convert(T,
        big"-.3845365626455519329296693060083837913367489520260829063809967396367023754075454122031")
    A[22, 11] = convert(T, big"0.")
    A[22, 12] = convert(T, big"0.")
    A[22, 13] = convert(T, big"0.")
    A[22, 14] = convert(T, big"0.")
    A[22, 15] = convert(T, big"0.")
    A[22, 16] = convert(T, big"0.")
    A[22, 17] = convert(T, big"0.")
    A[22, 18] = convert(T, big"0.")
    A[22, 19] = convert(T,
        big".3845365626455519329296693060083837913367489520260829063809967396367023754075454122031")
    A[22, 20] = convert(T,
        big".1316367023754075454122030740568234746157428970656730321378667908709827666511411271542")
    A[22, 21] = convert(T,
        big"-.1146176189444805429210887300555916540321998411667027651433109522778138762544220633889")
    A[23, 1] = convert(T,
        big".4938271604938271604938271604938271604938271604938271604938271604938271604938271604938e-1")
    A[23, 2] = convert(T,
        big".3950617283950617283950617283950617283950617283950617283950617283950617283950617283951")
    A[23, 3] = convert(T, big"0.")
    A[23, 4] = convert(T, big"0.")
    A[23, 5] = convert(T, big"0.")
    A[23, 6] = convert(T,
        big"-.6985294117647058823529411764705882352941176470588235294117647058823529411764705882353")
    A[23, 7] = convert(T,
        big"-.3834558992979908012587751149842653110626966836117162914548535463568143306705398208666")
    A[23, 8] = convert(T, big"0.")
    A[23, 9] = convert(T, big"0.")
    A[23, 10] = convert(T, big"0.")
    A[23, 11] = convert(T, big"0.")
    A[23, 12] = convert(T, big"0.")
    A[23, 13] = convert(T, big"0.")
    A[23, 14] = convert(T, big"0.")
    A[23, 15] = convert(T, big"0.")
    A[23, 16] = convert(T, big"0.")
    A[23, 17] = convert(T, big"0.")
    A[23, 18] = convert(T, big"0.")
    A[23, 19] = convert(T, big"0.")
    A[23, 20] = convert(T, big"0.")
    A[23, 21] = convert(T,
        big".3834558992979908012587751149842653110626966836117162914548535463568143306705398208666")
    A[23, 22] = convert(T,
        big".6985294117647058823529411764705882352941176470588235294117647058823529411764705882353")
    A[24, 1] = convert(T, big".25")
    A[24, 2] = convert(T, big"0.")
    A[24, 3] = convert(T,
        big"-.3179947624392068836513280957725402169846614291058735503179947624392068836513280957725")
    A[24, 4] = convert(T, big"0.")
    A[24, 5] = convert(T, big"0.")
    A[24, 6] = convert(T, big"0.")
    A[24, 7] = convert(T, big"0.")
    A[24, 8] = convert(T, big"0.")
    A[24, 9] = convert(T, big"0.")
    A[24, 10] = convert(T, big"0.")
    A[24, 11] = convert(T, big"0.")
    A[24, 12] = convert(T, big"0.")
    A[24, 13] = convert(T, big"0.")
    A[24, 14] = convert(T, big"0.")
    A[24, 15] = convert(T, big"0.")
    A[24, 16] = convert(T, big"0.")
    A[24, 17] = convert(T, big"0.")
    A[24, 18] = convert(T, big"0.")
    A[24, 19] = convert(T, big"0.")
    A[24, 20] = convert(T, big"0.")
    A[24, 21] = convert(T, big"0.")
    A[24, 22] = convert(T, big"0.")
    A[24, 23] = convert(T,
        big".3179947624392068836513280957725402169846614291058735503179947624392068836513280957725")
    A[25, 1] = convert(T,
        big"-8.372198884106881534693978918377517566370093845755858261918976097289568611513794096080")
    A[25, 2] = convert(T, big"-3.465")
    A[25, 3] = convert(T,
        big"-2.497530864197530864197530864197530864197530864197530864197530864197530864197530864198")
    A[25, 4] = convert(T, big"0.")
    A[25, 5] = convert(T, big"0.")
    A[25, 6] = convert(T, big"-1.6625")
    A[25, 7] = convert(T,
        big"-1.835440144855967078189300411522633744855967078189300411522633744855967078189300411523")
    A[25, 8] = convert(T, big"0.")
    A[25, 9] = convert(T,
        big"-73.24196602168055417777290626459588169783506828087515911825626050817981149401335285663")
    A[25, 10] = convert(T,
        big".8828637638234425252214582684159081671841568649549683123402374725184174104751169772484")
    A[25, 11] = convert(T,
        big"-.8814190544577027565498395466693637610438359285341071841488888129094808522593697679385")
    A[25, 12] = convert(T,
        big"-.4803279491030982366096634651792419931607860681361863738107180364716221880854580858817")
    A[25, 13] = convert(T,
        big"32.87910868410832718408338785621195992650771605408789380060168750394871935374300100483")
    A[25, 14] = convert(T,
        big"52.97402230205695709773274200747992549975863857695768897511204041642006328750715675109")
    A[25, 15] = convert(T,
        big"-8.886835693767662867753935454859468248687886743734693869236791808373822253213426903999")
    A[25, 16] = convert(T,
        big"2.509716434086569985363077956674238286724267256774032757480361264001771307725305753519")
    A[25, 17] = convert(T,
        big".1162475886937088360535666476065014131326443901611275096836079807311057386735832665899")
    A[25, 18] = convert(T,
        big".5839488303468939449260909132929399737902477241002934521537006256042283009612379572628")
    A[25, 19] = convert(T, big"1.581")
    A[25, 20] = convert(T, big"1.33584")
    A[25, 21] = convert(T,
        big"1.835440144855967078189300411522633744855967078189300411522633744855967078189300411523")
    A[25, 22] = convert(T, big"1.6625")
    A[25, 23] = convert(T,
        big"2.497530864197530864197530864197530864197530864197530864197530864197530864197530864198")
    A[25, 24] = convert(T, big"3.465")
    α[1] = convert(T,
        big".2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1")
    α[2] = convert(T, big"-.11")
    α[3] = convert(T, big"-.17")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big"0.")
    α[6] = convert(T, big"-.19")
    α[7] = convert(T, big"-.21")
    α[8] = convert(T, big"0.")
    α[9] = convert(T, big"-.23")
    α[10] = convert(T, big"-.27")
    α[11] = convert(T, big"-.29")
    α[12] = convert(T, big"0.")
    α[13] = convert(T,
        big".1384130236807829740053502031450331467488136400899412345912671194817223119377730668077")
    α[14] = convert(T,
        big".2158726906049313117089355111406811389654720741957730511230185948039919737765126474781")
    α[15] = convert(T,
        big".2438095238095238095238095238095238095238095238095238095238095238095238095238095238095")
    α[16] = convert(T,
        big".2158726906049313117089355111406811389654720741957730511230185948039919737765126474781")
    α[17] = convert(T,
        big".1384130236807829740053502031450331467488136400899412345912671194817223119377730668077")
    α[18] = convert(T, big".29")
    α[19] = convert(T, big".27")
    α[20] = convert(T, big".23")
    α[21] = convert(T, big".21")
    α[22] = convert(T, big".19")
    α[23] = convert(T, big".17")
    α[24] = convert(T, big".11")
    α[25] = convert(T,
        big".2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 12, stability_size = 3.028073429132015))
end

"""
Tableau form of Feagin12
"""
function constructFeagin12(T::Type = BigFloat)
    A = zeros(T, 25, 25)
    c = zeros(T, 25)
    α = zeros(T, 25)
    αEEst = zeros(T, 25)
    c[2] = convert(T, big"2")
    c[3] = convert(T,
        big".5555555555555555555555555555555555555555555555555555555555555555555555555555555555556")
    c[4] = convert(T,
        big".8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[5] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[6] = convert(T, big"1.")
    c[7] = convert(T,
        big".6718357091705138127122456610027975704389534205686825507102221981691883968854641411340")
    c[8] = convert(T,
        big".2887249411106202019354584889670249769081185983418069764696736835244986590217861541568")
    c[9] = convert(T, big".5625")
    c[10] = convert(T,
        big".8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[11] = convert(T,
        big".9476954311791992875623801621018367216495893258927406464583220401251379590957558619519")
    c[12] = convert(T,
        big".5481128768638026438877536748107544758421536129311287850283685092941115229776952973311e-1")
    c[13] = convert(T,
        big".8488805186071653506398389301626743020641481756400195420459339398355773991365476236893e-1")
    c[14] = convert(T,
        big".2655756032646428930981140590456168352972012641640776214486652703185222349414361456016")
    c[15] = convert(T, big".5")
    c[16] = convert(T,
        big".7344243967353571069018859409543831647027987358359223785513347296814777650585638543984")
    c[17] = convert(T,
        big".9151119481392834649360161069837325697935851824359980457954066060164422600863452376311")
    c[18] = convert(T,
        big".9476954311791992875623801621018367216495893258927406464583220401251379590957558619519")
    c[19] = convert(T,
        big".8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[20] = convert(T,
        big".2887249411106202019354584889670249769081185983418069764696736835244986590217861541568")
    c[21] = convert(T,
        big".6718357091705138127122456610027975704389534205686825507102221981691883968854641411340")
    c[22] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[23] = convert(T,
        big".5555555555555555555555555555555555555555555555555555555555555555555555555555555555556")
    c[24] = convert(T, big".2")
    c[25] = convert(T, big"1.")
    A[2, 1] = convert(T, big".2")
    A[3, 1] = convert(T,
        big"-.2160493827160493827160493827160493827160493827160493827160493827160493827160493827160")
    A[3, 2] = convert(T,
        big".7716049382716049382716049382716049382716049382716049382716049382716049382716049382716")
    A[4, 1] = convert(T,
        big".2083333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T, big".625")
    A[5, 1] = convert(T,
        big".1933333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T, big".22")
    A[5, 4] = convert(T, big"-.8e-1")
    A[6, 1] = convert(T, big".1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T, big".4")
    A[6, 5] = convert(T, big".5")
    A[7, 1] = convert(T,
        big".1033644716500104775703954356904817915433427083303498792441969232007020836692076572831")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".1240530945289467610615818892371153282110747849551802980440742772315189064222869368349")
    A[7, 5] = convert(T,
        big".4831711675610328992888364804519625087241092575172891773023800341818228904190296838151")
    A[7, 6] = convert(T,
        big"-.3875302456947632520856814437676205803957333023413680388042903644485548362506013679914e-1")
    A[8, 1] = convert(T,
        big".1240382614318333240819045859801751681400246706986336122924798440448279924947470898167")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big".2170506321979584863178462569531599428759163537577341676846566744497954872579456743665")
    A[8, 6] = convert(T,
        big".1374557920759667598129078018350481905944439909394085308429177535450188051849310154998e-1")
    A[8, 7] = convert(T,
        big"-.6610953172676828444558313414981495316726682520850165659175461032462670124939971157633e-1")
    A[9, 1] = convert(T,
        big".9147748948568829831449918469804321970888320999766601000904864070801948953424558529918e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big"-.5443485237174696899657549441448386113461568738470091780683178340728768152294577966066e-2")
    A[9, 7] = convert(T,
        big".6807168016884535185785151208951038631127517307587943722039518677930013251388913149245e-1")
    A[9, 8] = convert(T,
        big".4083943155826410467273068526538947800933031856649246445512393508534091461041598611744")
    A[10, 1] = convert(T,
        big".8900136525025510189545093554238417801432326974034341186926989993223983196257058017745e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big".4995282266455323601977934084206928004058911494068140919558100299611985086558381414011e-2")
    A[10, 7] = convert(T,
        big".3979182388198289973417396030013471560834350609314249708263044040250700525391549655752")
    A[10, 8] = convert(T,
        big".4279302107525766110681926083008979815582407305803964063123594678423525630853815189010")
    A[10, 9] = convert(T,
        big"-.8651176375578270057402774759550291032672463941289959659415853876594109934033211273435e-1")
    A[11, 1] = convert(T,
        big".6950876241349075431126939064098098227060210616855446152557575792953755238456142607555e-1")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T,
        big".1291469419001764619707595794827465511228717515014826340454871396367885195565200865413")
    A[11, 7] = convert(T,
        big"1.530736381023112950763425661432149390311775041124338743130106194678206265388420534343")
    A[11, 8] = convert(T,
        big".5778747611291400525467513494545767153348921004185718827180365180125884787253599132116")
    A[11, 9] = convert(T,
        big"-.9512947723210889805323408373888594539309244987992286480509495701319828569591060982197")
    A[11, 10] = convert(T,
        big"-.408276642965631951497484981519757463459627174520978426909934")
    A[12, 1] = convert(T,
        big".4448614032951358662694535070924635816201655010186841529333130115732839067844631360662e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T,
        big"-.3804768670569617319842326865745472030163315636268560657179642939163812931823178090201e-2")
    A[12, 7] = convert(T,
        big".1069550640296242007212626028090591544692060776449573995939720279979491092361490435771e-1")
    A[12, 8] = convert(T,
        big".2096162444999043332966742059289199208067346506600398980746522493272461502026769432938e-1")
    A[12, 9] = convert(T,
        big"-.2331460232593217866485614315519780776653378187560536038988470502127295139273620447039e-1")
    A[12, 10] = convert(T,
        big".263265981064536974369934736325334761174975280887405725010964e-2")
    A[12, 11] = convert(T,
        big".315472768977025060103545855572111407955208306374459723959783e-2")
    A[13, 1] = convert(T,
        big".1945888151197554755888010965253177612420737620162731862312147414806601584184678965302e-1")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T, big"0.")
    A[13, 5] = convert(T, big"0.")
    A[13, 6] = convert(T, big"0.")
    A[13, 7] = convert(T, big"0.")
    A[13, 8] = convert(T, big"0.")
    A[13, 9] = convert(T,
        big".6785129491718125093061216534523674761943647812591653323215344531729041102383826068881e-4")
    A[13, 10] = convert(T,
        big"-.4297958590492736232710053302301623435688633877248836036755502109954489170536483550959e-4")
    A[13, 11] = convert(T,
        big".1763589822602851554074859289533021399375534428299757341489807141089039042639321498044e-4")
    A[13, 12] = convert(T,
        big".6538666274150270510095952313851810335495113587873820983519242333986308816206310607574e-1")
    A[14, 1] = convert(T,
        big".2068368356642771059168281747982723610789091960434464115982313836045688158444651960319")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T, big"0.")
    A[14, 5] = convert(T, big"0.")
    A[14, 6] = convert(T, big"0.")
    A[14, 7] = convert(T, big"0.")
    A[14, 8] = convert(T, big"0.")
    A[14, 9] = convert(T,
        big".1667960671041564728280458666646964503063265050947925052155140101473800468812591484420e-1")
    A[14, 10] = convert(T,
        big"-.8795015632007102144570241782499865911302349902199592087049794487571848499522935439738e-2")
    A[14, 11] = convert(T,
        big".3466754553624639108244623152463792094275136540985964036372314865670539326325668110382e-2")
    A[14, 12] = convert(T,
        big"-.8612644601057176781614325622583512420302704989668912017992248333482802709963824294591")
    A[14, 13] = convert(T,
        big".9086518820740502810962394784692621450349571299392567891787847986693969945784247315140")
    A[15, 1] = convert(T,
        big".2039260846544840100915113146769256860385044495624130045623817918441487821545443195970e-1")
    A[15, 2] = convert(T, big"0.")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T, big"0.")
    A[15, 7] = convert(T, big"0.")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T,
        big".8694693920166859486754005555839475058339544609309409595773466506310797056354130304882e-1")
    A[15, 10] = convert(T,
        big"-.1916496304101498422864366117914050532871700766023376735876807254839218918835669149770e-1")
    A[15, 11] = convert(T,
        big".6556291594936632873648715732442445160348287552537460240988384514573072919062060289680e-2")
    A[15, 12] = convert(T,
        big".9874761281274347809037985286740338997389249680066322014454609882039267630836494384605e-1")
    A[15, 13] = convert(T,
        big".5353646955249960550832601736155674087171102472740210561183290254962131450211315663807e-2")
    A[15, 14] = convert(T,
        big".3011678640109679168370913038170516769200592297849574799980774547109414597317226366896")
    A[16, 1] = convert(T,
        big".2284104339177780995471154128930043987791369945969485457222835639366182967824509778177")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T, big"0.")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T,
        big"-.498707400793025250635016567442511512138603770959682292383042")
    A[16, 10] = convert(T,
        big".134841168335724478552596703792570104791700727205981058201689")
    A[16, 11] = convert(T,
        big"-.3874582440558341584399042269240292309351610591428068056743603625175611891589978761354e-1")
    A[16, 12] = convert(T,
        big"-1.274732574734748442403884308249089523809792927132503501996417388213689856389637269712")
    A[16, 13] = convert(T,
        big"1.439163644628771652011844524370380818752993035779118396305246406270252915801577043206")
    A[16, 14] = convert(T,
        big"-.2140074679679902542195035408273495696390280923448127954990268736031526728285351786283")
    A[16, 15] = convert(T,
        big".9582024177544302398927241391097813710599088746051536487680380575432052006086080693288")
    A[17, 1] = convert(T,
        big"2.002224776559742036142496460125067471214403062257117212098013934570904909077811687956")
    A[17, 2] = convert(T, big"0.")
    A[17, 3] = convert(T, big"0.")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T, big"0.")
    A[17, 7] = convert(T, big"0.")
    A[17, 8] = convert(T, big"0.")
    A[17, 9] = convert(T,
        big"2.067018099615249120919546564381385958254118596733416006795539073419461189789965414880")
    A[17, 10] = convert(T,
        big".6239781360861395419574712798314944661552923161670210806631400641515876319986835645504")
    A[17, 11] = convert(T,
        big"-.4622836855003114302832035541290620693919471018801127231857717578220394222019200167241e-1")
    A[17, 12] = convert(T,
        big"-8.849732883626496148600752467271189492866048354570927010946452349296049980461005208984")
    A[17, 13] = convert(T,
        big"7.742577078508559762274372257918355895601885907850371974336285358245514323924977674267")
    A[17, 14] = convert(T,
        big"-.5883585192508692109933533141277117456441258821309412028964590026164853179425282383760")
    A[17, 15] = convert(T,
        big"-1.106837333623806493957047080169530561761957696170148994429009662680169614076285221420")
    A[17, 16] = convert(T,
        big"-.9295290375792039997783972382912332142207880575118997475070736339961169400050824335694")
    A[18, 1] = convert(T,
        big"3.137895334120734429344516089898887968081612593303221002683164397131264854080655887290")
    A[18, 2] = convert(T, big"0.")
    A[18, 3] = convert(T, big"0.")
    A[18, 4] = convert(T, big"0.")
    A[18, 5] = convert(T, big"0.")
    A[18, 6] = convert(T,
        big".1291469419001764619707595794827465511228717515014826340454871396367885195565200865413")
    A[18, 7] = convert(T,
        big"1.530736381023112950763425661432149390311775041124338743130106194678206265388420534343")
    A[18, 8] = convert(T,
        big".5778747611291400525467513494545767153348921004185718827180365180125884787253599132116")
    A[18, 9] = convert(T,
        big"5.420882630551266830500568408918574219413005588518621564033611906468713101542334579174")
    A[18, 10] = convert(T,
        big".2315469260348293048726638008776436609048801808359456938369354937028526451467932046110")
    A[18, 11] = convert(T,
        big".7592929955789135601623013117852518735618013423331948952920597912045137042432676518814e-1")
    A[18, 12] = convert(T,
        big"-12.37299733801865132874145534025958065913498226175359059762562192965479091087793342733")
    A[18, 13] = convert(T,
        big"9.854558834647695439359572093173692020803677657217771019069796164343105374428968973275")
    A[18, 14] = convert(T,
        big".8591114313704365295793577090523677728899804951223296011591248114960893346142179981950e-1")
    A[18, 15] = convert(T,
        big"-5.652427528626439211171820900817627611803926026441892186739654808639279953585825515205")
    A[18, 16] = convert(T,
        big"-1.943009352428196108838337767823642877287248991241669204778734219860256733548749791369")
    A[18, 17] = convert(T,
        big"-.1283526018494045420184287143193446207421464913356123535599232759641139856465371475959")
    A[19, 1] = convert(T,
        big"1.383600544321960148785381182981677168251632684899225199955650618232490683056428985354")
    A[19, 2] = convert(T, big"0.")
    A[19, 3] = convert(T, big"0.")
    A[19, 4] = convert(T, big"0.")
    A[19, 5] = convert(T, big"0.")
    A[19, 6] = convert(T,
        big".4995282266455323601977934084206928004058911494068140919558100299611985086558381414011e-2")
    A[19, 7] = convert(T,
        big".3979182388198289973417396030013471560834350609314249708263044040250700525391549655752")
    A[19, 8] = convert(T,
        big".4279302107525766110681926083008979815582407305803964063123594678423525630853815189010")
    A[19, 9] = convert(T,
        big"-1.302991074244757709165514391230475733420714759983996459821461482811069921735213561975")
    A[19, 10] = convert(T,
        big".6612922786693770290971125281075130727345734122940080715006994611067693544004895081155")
    A[19, 11] = convert(T,
        big"-.1445597743069543497659693936887034639005858224415456555301446467509926053698427625331")
    A[19, 12] = convert(T,
        big"-6.965760347317982034678538674610839193567922481059192554608250635478358752496696344633")
    A[19, 13] = convert(T,
        big"6.658085432359917483534082955422104506321931975769351207164429345917517820065930603070")
    A[19, 14] = convert(T,
        big"-1.669973751088414864046958057255108450498079691992362275757972756183370386412619176110")
    A[19, 15] = convert(T,
        big"2.064137023180352638322890403018326471306046512239864521700896369142928855915137997534")
    A[19, 16] = convert(T,
        big"-.6747439626443064718629581295708377231920798759984050586488916899925211549437997031744")
    A[19, 17] = convert(T,
        big"-.1156188347949395004907036084359076100596057549353055820457292017095159857577078205324e-2")
    A[19, 18] = convert(T,
        big"-.544057908677007389319819914241631024660726585015012485938593e-2")
    A[20, 1] = convert(T,
        big".9512362970482876694746379758949735521669033789834754257582261243292060560082303478949")
    A[20, 2] = convert(T, big"0.")
    A[20, 3] = convert(T, big"0.")
    A[20, 4] = convert(T, big"0.")
    A[20, 5] = convert(T,
        big".2170506321979584863178462569531599428759163537577341676846566744497954872579456743665")
    A[20, 6] = convert(T,
        big".1374557920759667598129078018350481905944439909394085308429177535450188051849310154998e-1")
    A[20, 7] = convert(T,
        big"-.6610953172676828444558313414981495316726682520850165659175461032462670124939971157633e-1")
    A[20, 8] = convert(T, big"0.")
    A[20, 9] = convert(T,
        big".1522816967364144471366046970407471319214864326994221120996167853723749050279525002079")
    A[20, 10] = convert(T,
        big"-.3377410183575998408023007931339980043546434244575396676700801526788878440490902531919")
    A[20, 11] = convert(T,
        big"-.1928259816339957815349491992868244004693531106307879821211327084555949417498597974185e-1")
    A[20, 12] = convert(T,
        big"-3.682592696968668099324090155354996035763121207468888802018816928929650528231675747053")
    A[20, 13] = convert(T,
        big"3.161978704069820635415335284196838540183520803428870023313121375688872348961613906410")
    A[20, 14] = convert(T,
        big"-.3704625221068852907169918560220511254779434822840805691773856863151016773010882921618")
    A[20, 15] = convert(T,
        big"-.5149742003654404349964344566981279849411686164743168710203142306658077550892325779417e-1")
    A[20, 16] = convert(T,
        big"-.8296255321201529467870435417928484166593826752027206775365542784220425091994969746694e-3")
    A[20, 17] = convert(T,
        big".2798010414192785989865865890700275839613554026408795032135027653873641738192519718714e-5")
    A[20, 18] = convert(T,
        big".4186039164123602879698410207767884617941194406893561789422525187837946844315554254800e-1")
    A[20, 19] = convert(T,
        big".2790842550908773559156608745553796499662821675601262692902221878625437021870196271539")
    A[21, 1] = convert(T,
        big".1033644716500104775703954356904817915433427083303498792441969232007020836692076572831")
    A[21, 2] = convert(T, big"0.")
    A[21, 3] = convert(T, big"0.")
    A[21, 4] = convert(T,
        big".1240530945289467610615818892371153282110747849551802980440742772315189064222869368349")
    A[21, 5] = convert(T,
        big".4831711675610328992888364804519625087241092575172891773023800341818228904190296838151")
    A[21, 6] = convert(T,
        big"-.3875302456947632520856814437676205803957333023413680388042903644485548362506013679914e-1")
    A[21, 7] = convert(T, big"0.")
    A[21, 8] = convert(T,
        big"-.4383138203611224203910597889409601764206828366526006985800905251155989894984839457596")
    A[21, 9] = convert(T, big"0.")
    A[21, 10] = convert(T,
        big"-.2186366337216766476851114850171511993625093736982883305934862323914026063358863684281")
    A[21, 11] = convert(T,
        big"-.3123347643947192299816349952064403497661747596265781223230155751192843872703400896300e-1")
    A[21, 12] = convert(T, big"0.")
    A[21, 13] = convert(T, big"0.")
    A[21, 14] = convert(T, big"0.")
    A[21, 15] = convert(T, big"0.")
    A[21, 16] = convert(T, big"0.")
    A[21, 17] = convert(T, big"0.")
    A[21, 18] = convert(T,
        big".3123347643947192299816349952064403497661747596265781223230155751192843872703400896300e-1")
    A[21, 19] = convert(T,
        big".2186366337216766476851114850171511993625093736982883305934862323914026063358863684281")
    A[21, 20] = convert(T,
        big".4383138203611224203910597889409601764206828366526006985800905251155989894984839457596")
    A[22, 1] = convert(T,
        big".1933333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[22, 2] = convert(T, big"0.")
    A[22, 3] = convert(T, big".22")
    A[22, 4] = convert(T, big"-.8e-1")
    A[22, 5] = convert(T, big"0.")
    A[22, 6] = convert(T, big"0.")
    A[22, 7] = convert(T,
        big".9842561304993159281529002868560482433482025214912885759521434145952725510591247703880e-1")
    A[22, 8] = convert(T,
        big"-.1964108892230546534465265043901004176775390953401355324188487070169517581089448241223")
    A[22, 9] = convert(T, big"0.")
    A[22, 10] = convert(T,
        big".4364579304930687293918261225879491376096706767125250347633166283025775876654576372905")
    A[22, 11] = convert(T,
        big".6526137216757210985603709398055556983505438107084147167302701337621822207129608838308e-1")
    A[22, 12] = convert(T, big"0.")
    A[22, 13] = convert(T, big"0.")
    A[22, 14] = convert(T, big"0.")
    A[22, 15] = convert(T, big"0.")
    A[22, 16] = convert(T, big"0.")
    A[22, 17] = convert(T, big"0.")
    A[22, 18] = convert(T,
        big"-.6526137216757210985603709398055556983505438107084147167302701337621822207129608838308e-1")
    A[22, 19] = convert(T,
        big"-.4364579304930687293918261225879491376096706767125250347633166283025775876654576372905")
    A[22, 20] = convert(T,
        big".1964108892230546534465265043901004176775390953401355324188487070169517581089448241223")
    A[22, 21] = convert(T,
        big"-.9842561304993159281529002868560482433482025214912885759521434145952725510591247703880e-1")
    A[23, 1] = convert(T,
        big"-.2160493827160493827160493827160493827160493827160493827160493827160493827160493827160")
    A[23, 2] = convert(T,
        big".7716049382716049382716049382716049382716049382716049382716049382716049382716049382716")
    A[23, 3] = convert(T, big"0.")
    A[23, 4] = convert(T, big"0.")
    A[23, 5] = convert(T,
        big"-.6666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    A[23, 6] = convert(T, big"0.")
    A[23, 7] = convert(T,
        big"-.3906964692959784514469998022584959812490996652943959455591631110212407281436339846013")
    A[23, 8] = convert(T, big"0.")
    A[23, 9] = convert(T, big"0.")
    A[23, 10] = convert(T, big"0.")
    A[23, 11] = convert(T, big"0.")
    A[23, 12] = convert(T, big"0.")
    A[23, 13] = convert(T, big"0.")
    A[23, 14] = convert(T, big"0.")
    A[23, 15] = convert(T, big"0.")
    A[23, 16] = convert(T, big"0.")
    A[23, 17] = convert(T, big"0.")
    A[23, 18] = convert(T, big"0.")
    A[23, 19] = convert(T, big"0.")
    A[23, 20] = convert(T, big"0.")
    A[23, 21] = convert(T,
        big".3906964692959784514469998022584959812490996652943959455591631110212407281436339846013")
    A[23, 22] = convert(T,
        big".6666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    A[24, 1] = convert(T, big".2")
    A[24, 2] = convert(T, big"0.")
    A[24, 3] = convert(T,
        big"-.1646090534979423868312757201646090534979423868312757201646090534979423868312757201646")
    A[24, 4] = convert(T, big"0.")
    A[24, 5] = convert(T, big"0.")
    A[24, 6] = convert(T, big"0.")
    A[24, 7] = convert(T, big"0.")
    A[24, 8] = convert(T, big"0.")
    A[24, 9] = convert(T, big"0.")
    A[24, 10] = convert(T, big"0.")
    A[24, 11] = convert(T, big"0.")
    A[24, 12] = convert(T, big"0.")
    A[24, 13] = convert(T, big"0.")
    A[24, 14] = convert(T, big"0.")
    A[24, 15] = convert(T, big"0.")
    A[24, 16] = convert(T, big"0.")
    A[24, 17] = convert(T, big"0.")
    A[24, 18] = convert(T, big"0.")
    A[24, 19] = convert(T, big"0.")
    A[24, 20] = convert(T, big"0.")
    A[24, 21] = convert(T, big"0.")
    A[24, 22] = convert(T, big"0.")
    A[24, 23] = convert(T,
        big".1646090534979423868312757201646090534979423868312757201646090534979423868312757201646")
    A[25, 1] = convert(T,
        big"1.471787248811104084529495509890236112935353155185716919393995504060418834507026416303")
    A[25, 2] = convert(T, big".7875")
    A[25, 3] = convert(T,
        big".4212962962962962962962962962962962962962962962962962962962962962962962962962962962963")
    A[25, 4] = convert(T, big"0.")
    A[25, 5] = convert(T,
        big".2916666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    A[25, 6] = convert(T, big"0.")
    A[25, 7] = convert(T,
        big".3486007176283295632068544216296575692746899473678474657537568981834986072873525336473")
    A[25, 8] = convert(T,
        big".2294995447689948495828902337105554470738235696665067006625099245761293614167793479615")
    A[25, 9] = convert(T,
        big"5.790464857904819791598319781770034710982795060367224113331903192121494621770290062110")
    A[25, 10] = convert(T,
        big".4185875118565068688740737594265962072264614476042481510800164169962235305164482745364")
    A[25, 11] = convert(T,
        big".3070398802224740026496538174901066903892514823132139993866510031445598321274580577686")
    A[25, 12] = convert(T,
        big"-4.687009053506033322142563446838532480655744157947420404703040830648169300485343958299")
    A[25, 13] = convert(T,
        big"3.135716655938022621520381523998738565543954361999629154290906331430073621513072293611")
    A[25, 14] = convert(T,
        big"1.401348297109657208175105062756204410558450173139305083488955470547181175786794270770")
    A[25, 15] = convert(T,
        big"-5.529311014394990236290103060057643364212760557776581564009075823477850141179520237020")
    A[25, 16] = convert(T,
        big"-.8531382355080633493095468949747849061889275080395525195574976043372338240387571536437")
    A[25, 17] = convert(T,
        big".1035757803736101404118046071677727955182939144585001755737485164728465862790827590424")
    A[25, 18] = convert(T,
        big"-.1404744169506009411425469012021325348706659237000349571965457593133214062801025106415")
    A[25, 19] = convert(T,
        big"-.4185875118565068688740737594265962072264614476042481510800164169962235305164482745364")
    A[25, 20] = convert(T,
        big"-.2294995447689948495828902337105554470738235696665067006625099245761293614167793479615")
    A[25, 21] = convert(T,
        big"-.3486007176283295632068544216296575692746899473678474657537568981834986072873525336473")
    A[25, 22] = convert(T,
        big"-.2916666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    A[25, 23] = convert(T,
        big"-.4212962962962962962962962962962962962962962962962962962962962962962962962962962962963")
    A[25, 24] = convert(T, big"-.7875")
    α[1] = convert(T,
        big".2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1")
    α[2] = convert(T, big".234375e-1")
    α[3] = convert(T, big".3125e-1")
    α[4] = convert(T, big"0.")
    α[5] = convert(T,
        big".4166666666666666666666666666666666666666666666666666666666666666666666666666666666667e-1")
    α[6] = convert(T, big"0.")
    α[7] = convert(T, big".5e-1")
    α[8] = convert(T, big".5e-1")
    α[9] = convert(T, big"0.")
    α[10] = convert(T, big".1")
    α[11] = convert(T,
        big".7142857142857142857142857142857142857142857142857142857142857142857142857142857142857e-1")
    α[12] = convert(T, big"0.")
    α[13] = convert(T,
        big".1384130236807829740053502031450331467488136400899412345912671194817223119377730668077")
    α[14] = convert(T,
        big".2158726906049313117089355111406811389654720741957730511230185948039919737765126474781")
    α[15] = convert(T,
        big".2438095238095238095238095238095238095238095238095238095238095238095238095238095238095")
    α[16] = convert(T,
        big".2158726906049313117089355111406811389654720741957730511230185948039919737765126474781")
    α[17] = convert(T,
        big".1384130236807829740053502031450331467488136400899412345912671194817223119377730668077")
    α[18] = convert(T,
        big"-.7142857142857142857142857142857142857142857142857142857142857142857142857142857142857e-1")
    α[19] = convert(T, big"-.1")
    α[20] = convert(T, big"-.5e-1")
    α[21] = convert(T, big"-.5e-1")
    α[22] = convert(T,
        big"-.4166666666666666666666666666666666666666666666666666666666666666666666666666666666667e-1")
    α[23] = convert(T, big"-.3125e-1")
    α[24] = convert(T, big"-.234375e-1")
    α[25] = convert(T,
        big".2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 12, stability_size = 4.0696679617637175))
end

"""
Tableau form of Feagin14
"""
function constructFeagin14(T::Type = BigFloat)
    A = zeros(T, 35, 35)
    c = zeros(T, 35)
    α = zeros(T, 35)
    αEEst = zeros(T, 35)

    c[2] = convert(T,
        big".1111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
    c[3] = convert(T,
        big".5555555555555555555555555555555555555555555555555555555555555555555555555555555555556")
    c[4] = convert(T,
        big".8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[5] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[6] = convert(T, big"1.")
    c[7] = convert(T,
        big".6699869792727729217646837855059985139388452296384603532851421391683474428303956826239")
    c[8] = convert(T,
        big".2970683842138183573895847168082194132233320946989156873791682903324708698499266217383")
    c[9] = convert(T,
        big".7272727272727272727272727272727272727272727272727272727272727272727272727272727272727")
    c[10] = convert(T,
        big".1401527990421887652761874879669467176298064630825329362873230163439023340348096838456")
    c[11] = convert(T,
        big".7007010397701507371510998548307493379414070492655464089692218490447945746638665522966")
    c[12] = convert(T,
        big".3636363636363636363636363636363636363636363636363636363636363636363636363636363636364")
    c[13] = convert(T,
        big".2631578947368421052631578947368421052631578947368421052631578947368421052631578947368")
    c[14] = convert(T,
        big".392172246650270859125196642501208648863714315266128052078483e-1")
    c[15] = convert(T,
        big".8129175029283767629833931592780365061896123726172385507744269795906758195776958783707")
    c[16] = convert(T,
        big".1666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    c[17] = convert(T, big".9")
    c[18] = convert(T,
        big".6412992574519669233127711938966828094810966516150832254029235721305050295351572963693e-1")
    c[19] = convert(T,
        big".2041499092834288489277446343010234050271495052413337516288702042649259099754335560687")
    c[20] = convert(T,
        big".3953503910487605656156713698273243723522272974566594505545766538389345381768585023057")
    c[21] = convert(T,
        big".6046496089512394343843286301726756276477727025433405494454233461610654618231414976943")
    c[22] = convert(T,
        big".7958500907165711510722553656989765949728504947586662483711297957350740900245664439313")
    c[23] = convert(T,
        big".9358700742548033076687228806103317190518903348384916774597076427869494970464842703631")
    c[24] = convert(T,
        big".1666666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    c[25] = convert(T,
        big".8129175029283767629833931592780365061896123726172385507744269795906758195776958783707")
    c[26] = convert(T,
        big".392172246650270859125196642501208648863714315266128052078483e-1")
    c[27] = convert(T,
        big".3636363636363636363636363636363636363636363636363636363636363636363636363636363636364")
    c[28] = convert(T,
        big".7007010397701507371510998548307493379414070492655464089692218490447945746638665522966")
    c[29] = convert(T,
        big".1401527990421887652761874879669467176298064630825329362873230163439023340348096838456")
    c[30] = convert(T,
        big".2970683842138183573895847168082194132233320946989156873791682903324708698499266217383")
    c[31] = convert(T,
        big".6699869792727729217646837855059985139388452296384603532851421391683474428303956826239")
    c[32] = convert(T,
        big".3333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    c[33] = convert(T,
        big".5555555555555555555555555555555555555555555555555555555555555555555555555555555555556")
    c[34] = convert(T,
        big".1111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
    c[35] = convert(T, big"1.")
    A[2, 1] = convert(T,
        big".1111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
    A[3, 1] = convert(T,
        big"-.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[3, 2] = convert(T,
        big"1.388888888888888888888888888888888888888888888888888888888888888888888888888888888889")
    A[4, 1] = convert(T,
        big".2083333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[4, 2] = convert(T, big"0.")
    A[4, 3] = convert(T, big".625")
    A[5, 1] = convert(T,
        big".1933333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[5, 2] = convert(T, big"0.")
    A[5, 3] = convert(T, big".22")
    A[5, 4] = convert(T, big"-.8e-1")
    A[6, 1] = convert(T, big".1")
    A[6, 2] = convert(T, big"0.")
    A[6, 3] = convert(T, big"0.")
    A[6, 4] = convert(T, big".4")
    A[6, 5] = convert(T, big".5")
    A[7, 1] = convert(T,
        big".1034845616366797766729935465119103444997447982019713166066629728281981965079290745983")
    A[7, 2] = convert(T, big"0.")
    A[7, 3] = convert(T, big"0.")
    A[7, 4] = convert(T,
        big".1220688873064072225896440828689620771395927148341621347412746563709055937325311521675")
    A[7, 5] = convert(T,
        big".4825744903312466224751347801256881128659190238501680496794015023696413273862321544150")
    A[7, 6] = convert(T,
        big"-.3814096000156069997308862400056202056641130724784114774219699240039767479629669855696e-1")
    A[8, 1] = convert(T,
        big".1243805266540944128815164208687993162684914663596714231632892354628068537117612942798")
    A[8, 2] = convert(T, big"0.")
    A[8, 3] = convert(T, big"0.")
    A[8, 4] = convert(T, big"0.")
    A[8, 5] = convert(T,
        big".2261202821975843014222386629792029011967523207426331439651447460281196206643404356021")
    A[8, 6] = convert(T,
        big".1378858876180808806076958370164778145309694174914933853635428709475288586061552782365e-1")
    A[8, 7] = convert(T,
        big"-.6722101339966844497493995074143058569500863415253821828561997825320849038679063596730e-1")
    A[9, 1] = convert(T,
        big".9369190656596738155308854560830059338663496952177500856556033862893464429241815101000e-1")
    A[9, 2] = convert(T, big"0.")
    A[9, 3] = convert(T, big"0.")
    A[9, 4] = convert(T, big"0.")
    A[9, 5] = convert(T, big"0.")
    A[9, 6] = convert(T,
        big"-.6134068434505109872294989956416647356209145071288588710070986068372475355320835997035e-2")
    A[9, 7] = convert(T,
        big".2160198256255030637088600976598665734909794332781173201886676706066128640340557614360")
    A[9, 8] = convert(T,
        big".4236950635157619373376190739609767532058674695441235326831157041055522397561196508237")
    A[10, 1] = convert(T,
        big".8384798124090526646169687913728140859805331392249111310693346670107922625197375034871e-1")
    A[10, 2] = convert(T, big"0.")
    A[10, 3] = convert(T, big"0.")
    A[10, 4] = convert(T, big"0.")
    A[10, 5] = convert(T, big"0.")
    A[10, 6] = convert(T,
        big"-.1179493671009738143197550560312957753679619605907361507776128268875265788248790903515e-1")
    A[10, 7] = convert(T,
        big"-.2472990205688126523394738387431945983259928403533401326974984247503501083158412965835")
    A[10, 8] = convert(T,
        big".9780808583677290122593130140812916655037406554767339407565991037499621093437371932341e-1")
    A[10, 9] = convert(T,
        big".2175906892434206313600086517678603183441681200247821768799893467069296630467914197921")
    A[11, 1] = convert(T,
        big".6152553597694282279545623896143147143334239690648211074539397569215087099333654844097e-1")
    A[11, 2] = convert(T, big"0.")
    A[11, 3] = convert(T, big"0.")
    A[11, 4] = convert(T, big"0.")
    A[11, 5] = convert(T, big"0.")
    A[11, 6] = convert(T,
        big".5922327803245033080429900057980465247383895604442571368349896773084347972825775455007e-2")
    A[11, 7] = convert(T,
        big".4703261599638411122172243032058941134553625307461088250108483236601604516650193568134")
    A[11, 8] = convert(T,
        big".2996888638486790008539818370961923991368311216717812791841936858888827504094204242461")
    A[11, 9] = convert(T,
        big"-.2476568775939949146899922763298108258539580692639470955481886317480090967647905771626")
    A[11, 10] = convert(T,
        big".1108950297714376828939998518390617145224451736006787182086245987785252503880550245038")
    A[12, 1] = convert(T,
        big".4197000733627825798617928647872777872134836565431046112459945389674655429048057710370e-1")
    A[12, 2] = convert(T, big"0.")
    A[12, 3] = convert(T, big"0.")
    A[12, 4] = convert(T, big"0.")
    A[12, 5] = convert(T, big"0.")
    A[12, 6] = convert(T,
        big"-.317987696266205093901912847692712407988609169703103952205634e-2")
    A[12, 7] = convert(T,
        big".8063977149061920772608217115203795063935431115674197501197468839656405367779525213500")
    A[12, 8] = convert(T,
        big".9759831264123889790935228506842888513146720480030545503571875185550549213299958241991e-1")
    A[12, 9] = convert(T,
        big".7785755781583989090275124464529272389997634605941819649588520345133050850477185489203")
    A[12, 10] = convert(T,
        big".2048904238315994281894992020981056033120292350814206535748293420400885242747823516625")
    A[12, 11] = convert(T,
        big"-1.562615796274681883070709439505278252114628922364243608928053762634922556160297217820")
    A[13, 1] = convert(T,
        big".4377267822337301635744652424953398116882149670716141232569729223172939742940416733395e-1")
    A[13, 2] = convert(T, big"0.")
    A[13, 3] = convert(T, big"0.")
    A[13, 4] = convert(T, big"0.")
    A[13, 5] = convert(T, big"0.")
    A[13, 6] = convert(T, big"0.")
    A[13, 7] = convert(T, big"0.")
    A[13, 8] = convert(T, big"0.")
    A[13, 9] = convert(T,
        big".6243650275201952087943586285809336252816312169030959172012504609444028241438248581173e-2")
    A[13, 10] = convert(T,
        big".2000430971095773149944351654696478568290662322182649696087680691197048872391143823078")
    A[13, 11] = convert(T,
        big"-.8053283678049830368238571620489029119233928873370293148442058928084075077460302544840e-2")
    A[13, 12] = convert(T,
        big".2115175280673965219157119035233996013168778251575505730512208770404786743066139905871e-1")
    A[14, 1] = convert(T,
        big".2834992503635145630950235919207173122471376548964770977684956012393009143065795513785e-1")
    A[14, 2] = convert(T, big"0.")
    A[14, 3] = convert(T, big"0.")
    A[14, 4] = convert(T, big"0.")
    A[14, 5] = convert(T, big"0.")
    A[14, 6] = convert(T, big"0.")
    A[14, 7] = convert(T, big"0.")
    A[14, 8] = convert(T, big"0.")
    A[14, 9] = convert(T,
        big".2491632048558174075389491488059951494598846535854176800982219995075912885766744587193e-2")
    A[14, 10] = convert(T,
        big".2301387878545931496383998463737427687720871226381422342236583655735620108657836993957e-1")
    A[14, 11] = convert(T,
        big"-.3221559566929770987244760924671208781894636047606204610433085107190031098987004938258e-2")
    A[14, 12] = convert(T,
        big".9884425494476646689463354144878852560408199827860146481292993078049373245839618405001e-2")
    A[14, 13] = convert(T,
        big"-.2130107713288873513843076428759273848866345654295724666320922464722154754985568313136e-1")
    A[15, 1] = convert(T,
        big".3435118942902430010494322347351479430833531749807014262686507474123120416010457867571")
    A[15, 2] = convert(T, big"0.")
    A[15, 3] = convert(T, big"0.")
    A[15, 4] = convert(T, big"0.")
    A[15, 5] = convert(T, big"0.")
    A[15, 6] = convert(T, big"0.")
    A[15, 7] = convert(T, big"0.")
    A[15, 8] = convert(T, big"0.")
    A[15, 9] = convert(T,
        big".2104519120236273856090970119990106557888074052256267000419050051487632641518018732685")
    A[15, 10] = convert(T,
        big"1.034274520572304119364829268288257099386679996983247401666929134177931632176349026735")
    A[15, 11] = convert(T,
        big".6003036458644224870512404482066405749390780924061569454673075686417142117164254262878e-2")
    A[15, 12] = convert(T,
        big".8559381250996195375780121060024077289150626526164160058172684354881277648341960563008")
    A[15, 13] = convert(T,
        big"-.9772350050367668108722648523725256330131076568928396776974412446349105799705851506077")
    A[15, 14] = convert(T,
        big"-.6600269804792946946162250138563276937205739812199748747775581736879654453322759683463")
    A[16, 1] = convert(T,
        big"-.1435740016721680695382063999350763666577559543783998809757153672896315044426183882232e-1")
    A[16, 2] = convert(T, big"0.")
    A[16, 3] = convert(T, big"0.")
    A[16, 4] = convert(T, big"0.")
    A[16, 5] = convert(T, big"0.")
    A[16, 6] = convert(T, big"0.")
    A[16, 7] = convert(T, big"0.")
    A[16, 8] = convert(T, big"0.")
    A[16, 9] = convert(T,
        big"-.3662532700490399702936857968489747917331190817335522078657913621382824038988807796287e-1")
    A[16, 10] = convert(T,
        big".3502549756362136819768494069798465243467890824711035742020654749717518291597210559354e-1")
    A[16, 11] = convert(T,
        big".3609460163621135089317866587583352398236899298642376718895880083960486970547825683491e-1")
    A[16, 12] = convert(T,
        big"-.2652199675536811063515959468346019236496270124574642848667252606942739787160130682669e-1")
    A[16, 13] = convert(T,
        big".4456990113056981196389115375088399081043363230822267716707629092315111479614958673268e-1")
    A[16, 14] = convert(T,
        big".1243430933313582432862255957417864480389734088951067419167759990001419776217292554191")
    A[16, 15] = convert(T,
        big".4138296932394806944035124962043359604261929086744760344472227418812310333088685698274e-2")
    A[17, 1] = convert(T,
        big".3560324044251202909756091163980891762641062223797488026536968101501275805051289760823")
    A[17, 2] = convert(T, big"0.")
    A[17, 3] = convert(T, big"0.")
    A[17, 4] = convert(T, big"0.")
    A[17, 5] = convert(T, big"0.")
    A[17, 6] = convert(T, big"0.")
    A[17, 7] = convert(T, big"0.")
    A[17, 8] = convert(T, big"0.")
    A[17, 9] = convert(T,
        big"-.450192758947562595966821779075956175110645100214763601190349")
    A[17, 10] = convert(T,
        big".430527907083710898626656292808782917793030154094709462877146")
    A[17, 11] = convert(T,
        big".5119730290110222376685569603940716920771257870306513863906805244405042813755411512463")
    A[17, 12] = convert(T,
        big".9083036388864042603901591246381102139974962148199046305445452541870528153933236088278")
    A[17, 13] = convert(T,
        big"-1.239210933719339317573724691515340288544138892486057261860887966510000755220957594942")
    A[17, 14] = convert(T,
        big"-.6490486616717614651416723488790625539054028319671910976544025456235491510559878435372")
    A[17, 15] = convert(T,
        big".2517089045868192922104805299489705414048878529314474912189256354259853776829630937658")
    A[17, 16] = convert(T,
        big".7799064703455863988107567952823344760235405934115501870206452879298798513199886085571")
    A[18, 1] = convert(T,
        big".1309356874065130664068812064188349801274704382131924878449566575565302965696195341197e-1")
    A[18, 2] = convert(T, big"0.")
    A[18, 3] = convert(T, big"0.")
    A[18, 4] = convert(T, big"0.")
    A[18, 5] = convert(T, big"0.")
    A[18, 6] = convert(T, big"0.")
    A[18, 7] = convert(T, big"0.")
    A[18, 8] = convert(T, big"0.")
    A[18, 9] = convert(T, big"0.")
    A[18, 10] = convert(T, big"0.")
    A[18, 11] = convert(T, big"0.")
    A[18, 12] = convert(T, big"0.")
    A[18, 13] = convert(T,
        big"-.9320530679851139459084619627671082378586315096846671421247697017556505173897578610165e-4")
    A[18, 14] = convert(T,
        big".5053743342622993596400904431385907267709423447161223817027456630856526555478831396014e-1")
    A[18, 15] = convert(T,
        big".8044703419444879791095791096101977976413118689308653610493721999399129417586629251430e-6")
    A[18, 16] = convert(T,
        big".5917260294941711905287557427777172598443409719243215281782302034071342229921661278343e-3")
    A[18, 17] = convert(T,
        big"-.4016147221545573370646916849063755877322642479500938046774565993013424294867398455789e-6")
    A[19, 1] = convert(T,
        big".2079264844660530125419445440007656521672552061443734079797586969853055549175505457737e-1")
    A[19, 2] = convert(T, big"0.")
    A[19, 3] = convert(T, big"0.")
    A[19, 4] = convert(T, big"0.")
    A[19, 5] = convert(T, big"0.")
    A[19, 6] = convert(T, big"0.")
    A[19, 7] = convert(T, big"0.")
    A[19, 8] = convert(T, big"0.")
    A[19, 9] = convert(T, big"0.")
    A[19, 10] = convert(T, big"0.")
    A[19, 11] = convert(T, big"0.")
    A[19, 12] = convert(T, big"0.")
    A[19, 13] = convert(T,
        big".5826959188000859151019026978372841089514061030298715701031065480360641416298102920851e-3")
    A[19, 14] = convert(T,
        big"-.8017007323588159390833421865258527466405584659196335246554992680506588169863285718822e-2")
    A[19, 15] = convert(T,
        big".4038476438471369403751708217435605704841172903308955066191655368223862388605213690921e-5")
    A[19, 16] = convert(T,
        big".8546099980555061442250561145675356025101146220336224918025961310211940592009621595606e-1")
    A[19, 17] = convert(T,
        big"-.2044864809358042427067075696910043079044428375526774562331430989116458814609927891477e-5")
    A[19, 18] = convert(T,
        big".1053285788244318933997994029790939973542409042351728431465827473723673651882417656762")
    A[20, 1] = convert(T,
        big"1.401534497957360214154462473557713067184864529175977331289881318884096354294079099114")
    A[20, 2] = convert(T, big"0.")
    A[20, 3] = convert(T, big"0.")
    A[20, 4] = convert(T, big"0.")
    A[20, 5] = convert(T, big"0.")
    A[20, 6] = convert(T, big"0.")
    A[20, 7] = convert(T, big"0.")
    A[20, 8] = convert(T, big"0.")
    A[20, 9] = convert(T, big"0.")
    A[20, 10] = convert(T, big"0.")
    A[20, 11] = convert(T, big"0.")
    A[20, 12] = convert(T, big"0.")
    A[20, 13] = convert(T,
        big"-.2302520009842212616162724103674156212611302982744556219175010157057031125814669239016")
    A[20, 14] = convert(T,
        big"-7.211068404669129056595822371068742471658564935099615697324849532576890894506619405031")
    A[20, 15] = convert(T,
        big".3729015606948363352369953278521323402177595666786623882373057096229137360164435411243e-2")
    A[20, 16] = convert(T,
        big"-4.714154957271250206787781793922247570113233732218200980194845522013711035054762664884")
    A[20, 17] = convert(T,
        big"-.1763676575453492420538419950327976735749038866956001340593194717236122233799126229446e-2")
    A[20, 18] = convert(T,
        big"7.641305480386987655630293108802376511851733678139370059818519661401442202665741111270")
    A[20, 19] = convert(T,
        big"3.506020436597518349898960829497447109682129498933753736341591881470708008233521976557")
    A[21, 1] = convert(T,
        big"11.95146506941206867993723858307164016744736108265535168242754934626543968357331742096")
    A[21, 2] = convert(T, big"0.")
    A[21, 3] = convert(T, big"0.")
    A[21, 4] = convert(T, big"0.")
    A[21, 5] = convert(T, big"0.")
    A[21, 6] = convert(T, big"0.")
    A[21, 7] = convert(T, big"0.")
    A[21, 8] = convert(T, big"0.")
    A[21, 9] = convert(T, big"0.")
    A[21, 10] = convert(T, big"0.")
    A[21, 11] = convert(T, big"0.")
    A[21, 12] = convert(T, big"0.")
    A[21, 13] = convert(T,
        big"7.794809321081759687835167002317643882202842795989809549197917776161588225206322580459")
    A[21, 14] = convert(T,
        big"-56.45013938673257925235609911209042814404681000613405538635967763011214022629172907669")
    A[21, 15] = convert(T,
        big".9123763069306449013445304492902766457096074504036737047499704936582270274950128398912e-1")
    A[21, 16] = convert(T,
        big"-12.73362799254348862019455243091992750381627175299189605168457824373779389828110581300")
    A[21, 17] = convert(T,
        big"-.3968959219047197123135428109397366747123830704331478729319411886202118671113516172493e-1")
    A[21, 18] = convert(T,
        big"54.43921418835708869962257651553077918614383784233053341001985423053366890118247056463")
    A[21, 19] = convert(T,
        big"-3.644116379215692368464069903613506458067214784092667356589342345057374050114156075061")
    A[21, 20] = convert(T,
        big"-.8045032499105099108990307879585794993156949132107878807481027183961246894903442258757")
    A[22, 1] = convert(T,
        big"-148.8094265071004884278388682686476255619306120821485965777899951377767737092911763254")
    A[22, 2] = convert(T, big"0.")
    A[22, 3] = convert(T, big"0.")
    A[22, 4] = convert(T, big"0.")
    A[22, 5] = convert(T, big"0.")
    A[22, 6] = convert(T, big"0.")
    A[22, 7] = convert(T, big"0.")
    A[22, 8] = convert(T, big"0.")
    A[22, 9] = convert(T, big"0.")
    A[22, 10] = convert(T, big"0.")
    A[22, 11] = convert(T, big"0.")
    A[22, 12] = convert(T, big"0.")
    A[22, 13] = convert(T,
        big"-91.72952782912564843579356624023216234952287290363542836291360346578688265538801398361")
    A[22, 14] = convert(T,
        big"707.6561449715983598345757192863357161548211289666495623584804744987957677893379157809")
    A[22, 15] = convert(T,
        big"-1.10563611857482440905296961311590930801338308942637769555540")
    A[22, 16] = convert(T,
        big"176.1345918838113725878598980760556604069995167623016865882869129962911416096097878945")
    A[22, 17] = convert(T,
        big".4913848242148806622688983451644545574168846314027647925019604519368994965045299923826")
    A[22, 18] = convert(T,
        big"-684.2780004498149443582375356108950819560771678936002751371799726829821841834791232605")
    A[22, 19] = convert(T,
        big"27.99106049983982589842243321243804074460025184006686868209688958109916979926727384229")
    A[22, 20] = convert(T,
        big"13.19397100302823334436709643711532384350641596237449753683872220663989495376087330358")
    A[22, 21] = convert(T,
        big"1.251287812839804454501149741480560063172688300773964063605141347518040989702499199856")
    A[23, 1] = convert(T,
        big"-9.673079469481967636441261184332193958399514085718772596349277868068021458303626779169")
    A[23, 2] = convert(T, big"0.")
    A[23, 3] = convert(T, big"0.")
    A[23, 4] = convert(T, big"0.")
    A[23, 5] = convert(T, big"0.")
    A[23, 6] = convert(T, big"0.")
    A[23, 7] = convert(T, big"0.")
    A[23, 8] = convert(T, big"0.")
    A[23, 9] = convert(T, big"0.")
    A[23, 10] = convert(T, big"0.")
    A[23, 11] = convert(T, big"0.")
    A[23, 12] = convert(T, big"0.")
    A[23, 13] = convert(T,
        big"-4.469901508585055314438462277019603604978306814087514357488023393670679083633020106516")
    A[23, 14] = convert(T,
        big"45.51271286909526819682419504000527511789059078173984816890412459840121969200961260987")
    A[23, 15] = convert(T,
        big"-.713085086183826912791492024438246129930559805352394367050813e-1")
    A[23, 16] = convert(T,
        big"11.22736140684127415825906244799393842078268007767944830815221105133516977144595052189")
    A[23, 17] = convert(T,
        big".1262443767176227245162379129091388093617868898191054263714925416869147773104813482457")
    A[23, 18] = convert(T,
        big"-43.54393395494833136058106249072421076238143044676214056937881652359375369765457150165")
    A[23, 19] = convert(T,
        big".7871743075430589783987929949965509020645460914432340378113766124779028133099797867162")
    A[23, 20] = convert(T,
        big".5322646967446842156693007086038866907853957768215038536520118921656033723449302296244")
    A[23, 21] = convert(T,
        big".4224227339963253260102251274713887725750865388096033468497941673910509540050957057177")
    A[23, 22] = convert(T,
        big".8591312495030671073084380314998594434411150562941549563989586466154235621165245563192e-1")
    A[24, 1] = convert(T,
        big"-10.06640324470547024033966069004268914722028247579687652710623604380152449409080444899")
    A[24, 2] = convert(T, big"0.")
    A[24, 3] = convert(T, big"0.")
    A[24, 4] = convert(T, big"0.")
    A[24, 5] = convert(T, big"0.")
    A[24, 6] = convert(T, big"0.")
    A[24, 7] = convert(T, big"0.")
    A[24, 8] = convert(T, big"0.")
    A[24, 9] = convert(T,
        big"-.3662532700490399702936857968489747917331190817335522078657913621382824038988807796287e-1")
    A[24, 10] = convert(T,
        big".3502549756362136819768494069798465243467890824711035742020654749717518291597210559354e-1")
    A[24, 11] = convert(T,
        big".3609460163621135089317866587583352398236899298642376718895880083960486970547825683491e-1")
    A[24, 12] = convert(T,
        big"-.2652199675536811063515959468346019236496270124574642848667252606942739787160130682669e-1")
    A[24, 13] = convert(T,
        big"-6.270889721814641435905531494788716038393561229573960230194057818533161624674313994502")
    A[24, 14] = convert(T,
        big"48.20792374425629890907021030081950639234925931416361161278899187780407980462426656808")
    A[24, 15] = convert(T,
        big"-.694471689136165640882395180583732834557754169149088630301342e-1")
    A[24, 16] = convert(T,
        big"12.68106902048502956983413709136098070661084838114121251454273060707937017246509534894")
    A[24, 17] = convert(T,
        big".119671168968323754838161435501011294100927813964199613229864e-1")
    A[24, 18] = convert(T,
        big"-46.72497649924824080033582682426626955932013216597956070401309263301039263373634230581")
    A[24, 19] = convert(T,
        big"1.330296133266267113147100392982165913990335111912271192356479099067512051132965697343")
    A[24, 20] = convert(T,
        big"1.007667875033982983534389036199266577711627177936617199056121787956529680139072027935")
    A[24, 21] = convert(T,
        big".2095120519336650916641223884754807028927707538644872411247284065032940106679251005781e-1")
    A[24, 22] = convert(T,
        big".2101347063312641773177354243313964074244121884437574908902263894855162847478911411134e-1")
    A[24, 23] = convert(T,
        big".9521960144171217941751015424545759073763602336583562405468424451848266905185171865534e-2")
    A[25, 1] = convert(T,
        big"-409.4780816777437087725890974093703576244243416067520683455326035855162023776088699896")
    A[25, 2] = convert(T, big"0.")
    A[25, 3] = convert(T, big"0.")
    A[25, 4] = convert(T, big"0.")
    A[25, 5] = convert(T, big"0.")
    A[25, 6] = convert(T, big"0.")
    A[25, 7] = convert(T, big"0.")
    A[25, 8] = convert(T, big"0.")
    A[25, 9] = convert(T,
        big".2104519120236273856090970119990106557888074052256267000419050051487632641518018732685")
    A[25, 10] = convert(T,
        big"1.034274520572304119364829268288257099386679996983247401666929134177931632176349026735")
    A[25, 11] = convert(T,
        big".6003036458644224870512404482066405749390780924061569454673075686417142117164254262878e-2")
    A[25, 12] = convert(T,
        big".8559381250996195375780121060024077289150626526164160058172684354881277648341960563008")
    A[25, 13] = convert(T,
        big"-250.5169985474478604927776577293161303865840504207820779326393997812026874735614210230")
    A[25, 14] = convert(T,
        big"1946.424666523884277660537503282647585958298508957614274560610260899186136259514015246")
    A[25, 15] = convert(T,
        big"-3.045038821023103655061058090868608827869505440976021016842196622317831446605499698935")
    A[25, 16] = convert(T,
        big"490.6263795282817135212082652991680838415985422740616633051003594128766152337185220086")
    A[25, 17] = convert(T,
        big"1.566475895312709071154840670135974457395956152459667753199388690841173424714434871921")
    A[25, 18] = convert(T,
        big"-1881.974289940111733622172673770358706192159066384530557689275696031792911993357071098")
    A[25, 19] = convert(T,
        big"75.25922247248471752788377136433031498216206189142459440229301807516615379972994062700")
    A[25, 20] = convert(T,
        big"34.57343569803310676224343447365546896967286447935510158001529990937243976348724448442")
    A[25, 21] = convert(T,
        big"3.211476794409689614354173618470737551690229667488916278855754113243135684398993410117")
    A[25, 22] = convert(T,
        big"-.4604080417384143913072014042370588488672450952653828208427296561415079214017074427602")
    A[25, 23] = convert(T,
        big"-.8707183398418105224318841379579862457242520473889365722145748143125162133630944128398e-1")
    A[25, 24] = convert(T,
        big"-7.393518141583030675670169521955210639991857732491329543926346613193825315394087286297")
    A[26, 1] = convert(T,
        big"3.433474758535508789210934962575967811206238910720084588712755786644583035514752699598")
    A[26, 2] = convert(T, big"0.")
    A[26, 3] = convert(T, big"0.")
    A[26, 4] = convert(T, big"0.")
    A[26, 5] = convert(T, big"0.")
    A[26, 6] = convert(T, big"0.")
    A[26, 7] = convert(T, big"0.")
    A[26, 8] = convert(T, big"0.")
    A[26, 9] = convert(T,
        big".2491632048558174075389491488059951494598846535854176800982219995075912885766744587193e-2")
    A[26, 10] = convert(T,
        big".2301387878545931496383998463737427687720871226381422342236583655735620108657836993957e-1")
    A[26, 11] = convert(T,
        big"-.3221559566929770987244760924671208781894636047606204610433085107190031098987004938258e-2")
    A[26, 12] = convert(T,
        big".9884425494476646689463354144878852560408199827860146481292993078049373245839618405001e-2")
    A[26, 13] = convert(T,
        big"2.162527993779225077883078419047573540457592253357327094851479956564246957314476133478")
    A[26, 14] = convert(T,
        big"-16.26998645464574213280656406601394890069875520402288517985775075363232756881970486667")
    A[26, 15] = convert(T,
        big"-.1285345021205245528435834174709350105380290375426545062302651848844352856037884822181")
    A[26, 16] = convert(T,
        big"-8.98915042666504253089307820833379330486511746063552853023189")
    A[26, 17] = convert(T,
        big"-.3485953632320253333870802018510136501924017672505137649688730136175086767654181319387e-2")
    A[26, 18] = convert(T,
        big"15.79361941133398075362351873886955741358533870251397376656158275266140525531011608606")
    A[26, 19] = convert(T,
        big"-.5744033309140950656281654820173358201483836631956754708231458398423255984252281047127")
    A[26, 20] = convert(T,
        big"-.3456020390213932966927224966081249825352372288276553067081833889419898565070467534157")
    A[26, 21] = convert(T,
        big"-.6622414902065850917316199913837577811330679927074186873906450413385445874036001388495e-2")
    A[26, 22] = convert(T,
        big"-.7777881292422041640325464586073643097593472096267591120155367761150273183248441708392e-2")
    A[26, 23] = convert(T,
        big"-.3560841924022749133388272326974373646752408187917065879526063406092336300493607300593e-2")
    A[26, 24] = convert(T,
        big"4.792825064499307996497977496298401894572969341393590555417712618624354747222657791607")
    A[26, 25] = convert(T,
        big".153725464873068577844576387402512082757034273069877432944621")
    A[27, 1] = convert(T,
        big"32.30385208719854423269947344400315350913649750477846297617061421719281146058139852238")
    A[27, 2] = convert(T, big"0.")
    A[27, 3] = convert(T, big"0.")
    A[27, 4] = convert(T, big"0.")
    A[27, 5] = convert(T, big"0.")
    A[27, 6] = convert(T,
        big"-.317987696266205093901912847692712407988609169703103952205634e-2")
    A[27, 7] = convert(T,
        big".8063977149061920772608217115203795063935431115674197501197468839656405367779525213500")
    A[27, 8] = convert(T,
        big".9759831264123889790935228506842888513146720480030545503571875185550549213299958241991e-1")
    A[27, 9] = convert(T,
        big".7785755781583989090275124464529272389997634605941819649588520345133050850477185489203")
    A[27, 10] = convert(T,
        big".2048904238315994281894992020981056033120292350814206535748293420400885242747823516625")
    A[27, 11] = convert(T,
        big"-1.562615796274681883070709439505278252114628922364243608928053762634922556160297217820")
    A[27, 12] = convert(T, big"0.")
    A[27, 13] = convert(T,
        big"16.34298918823105706485042439739271747087533535041545512917666902744198799725970841669")
    A[27, 14] = convert(T,
        big"-154.5445552935436212307301896314710363993166836696091165017078152549564923882084122674")
    A[27, 15] = convert(T,
        big"1.569710887033348726920342834176217614662635935824970859658624964687079589089479471888")
    A[27, 16] = convert(T,
        big"3.276855450872481313214298172699007311655224049747336000450385269517693130775985884604")
    A[27, 17] = convert(T,
        big"-.5034892451936531763480407271997836265340810956916323972462042700071863164675818955838e-1")
    A[27, 18] = convert(T,
        big"153.3211518580416650705937678859146940112243631025945564907021486707139114294996134941")
    A[27, 19] = convert(T,
        big"7.175681863277204958467664848147841435678263080348653386540185145833155908488128910568")
    A[27, 20] = convert(T,
        big"-2.940367486753004819459176598969309892153205943807775979427615740476908865098135595635")
    A[27, 21] = convert(T,
        big"-.6658459460768031444707496760226288702819204931972568878708744783028558369468497032253e-1")
    A[27, 22] = convert(T,
        big"-.4623460549908436612292486685622172611769665140168592842374268449140643068786760618896e-1")
    A[27, 23] = convert(T,
        big"-.2041987335856794015393882286172697788485797748215817776751235910664984352284968100100e-1")
    A[27, 24] = convert(T,
        big"-53.35231064387358505159534411659981079740450904957915977996876390672711239156977103431")
    A[27, 25] = convert(T,
        big"-1.355487147150786549787321867059964040175545016141913251148206738329360142936656282958")
    A[27, 26] = convert(T,
        big"-1.571962758012327518829017351714592491776872191144425834618663282570958684038698495739")
    A[28, 1] = convert(T,
        big"-16.64514674863415128720312944039317587645603711308189782044257016154825923946758475845")
    A[28, 2] = convert(T, big"0.")
    A[28, 3] = convert(T, big"0.")
    A[28, 4] = convert(T, big"0.")
    A[28, 5] = convert(T, big"0.")
    A[28, 6] = convert(T,
        big".5922327803245033080429900057980465247383895604442571368349896773084347972825775455007e-2")
    A[28, 7] = convert(T,
        big".4703261599638411122172243032058941134553625307461088250108483236601604516650193568134")
    A[28, 8] = convert(T,
        big".2996888638486790008539818370961923991368311216717812791841936858888827504094204242461")
    A[28, 9] = convert(T,
        big"-.2476568775939949146899922763298108258539580692639470955481886317480090967647905771626")
    A[28, 10] = convert(T,
        big".1108950297714376828939998518390617145224451736006787182086245987785252503880550245038")
    A[28, 11] = convert(T, big"0.")
    A[28, 12] = convert(T,
        big"-.4917190438462291470706666287041940976780819072106730449888664749836403474888832394921")
    A[28, 13] = convert(T,
        big"-11.47431544272894969683894925643525363508424541308531757856483965863898534849416840511")
    A[28, 14] = convert(T,
        big"80.25931665762302725417024858864844001527933666235899875893849400507278534931158408231")
    A[28, 15] = convert(T,
        big"-.3841323039800428476253125267590291037469268413420882192068133107492120348263618466046")
    A[28, 16] = convert(T,
        big"7.281476674681075834713269509261361157676125818628777243483988994104498714011047355205")
    A[28, 17] = convert(T,
        big"-.1326993846122483795105717081760352748368273416167518843018178653526280269065470590467")
    A[28, 18] = convert(T,
        big"-81.07998325257307266746792897522552400060707166336329885641562357237166810196760593013")
    A[28, 19] = convert(T,
        big"-1.250374928356206395217681856561791199622537474924031863192434629401819729868852090550")
    A[28, 20] = convert(T,
        big"2.592635949695436810237763795043773249942264473592968880837586883560068434349818491911")
    A[28, 21] = convert(T,
        big"-.3014402983464045398301639972605268752644315372756414953420797074457552586137488110716")
    A[28, 22] = convert(T,
        big".2213844607898323374517064515727737916952468390573184143179573617704323166985265217363")
    A[28, 23] = convert(T,
        big".8275772747718929319559898709746931529962764354298098905497078729734353980896315305691e-1")
    A[28, 24] = convert(T,
        big"18.99606620406115204646724500372432639981751614122371589366718674999943569769696943522")
    A[28, 25] = convert(T,
        big".2692319464096396856234680151283341674600519103489128451211866688910668614577677735665")
    A[28, 26] = convert(T,
        big"1.626748274470665374629893649296289339881250292841836802790201430504847697803528636395")
    A[28, 27] = convert(T,
        big".4917190438462291470706666287041940976780819072106730449888664749836403474888832394921")
    A[29, 1] = convert(T,
        big".8384798124090526646169687913728140859805331392249111310693346670107922625197375034871e-1")
    A[29, 2] = convert(T, big"0.")
    A[29, 3] = convert(T, big"0.")
    A[29, 4] = convert(T, big"0.")
    A[29, 5] = convert(T, big"0.")
    A[29, 6] = convert(T,
        big"-.1179493671009738143197550560312957753679619605907361507776128268875265788248790903515e-1")
    A[29, 7] = convert(T,
        big"-.2472990205688126523394738387431945983259928403533401326974984247503501083158412965835")
    A[29, 8] = convert(T,
        big".9780808583677290122593130140812916655037406554767339407565991037499621093437371932341e-1")
    A[29, 9] = convert(T,
        big".2175906892434206313600086517678603183441681200247821768799893467069296630467914197921")
    A[29, 10] = convert(T, big"0.")
    A[29, 11] = convert(T,
        big".1375856067633252248656596321967877466474472229750848659754400903987833771639575727867")
    A[29, 12] = convert(T,
        big".4398702297150466850587900923415450260461038902942613590425808839943205635447284745074e-1")
    A[29, 13] = convert(T, big"0.")
    A[29, 14] = convert(T,
        big"-.5137008137681933419570044566186303037387573636419640300869712169933398305905931343468")
    A[29, 15] = convert(T,
        big".8263556911513155086442113083991534587014231586161685769224194977471882335420141183213")
    A[29, 16] = convert(T,
        big"25.70181397198118326258738829725199395111365563419600781824702737091645129169813134401")
    A[29, 17] = convert(T, big"0.")
    A[29, 18] = convert(T, big"0.")
    A[29, 19] = convert(T, big"0.")
    A[29, 20] = convert(T, big"0.")
    A[29, 21] = convert(T, big"0.")
    A[29, 22] = convert(T, big"0.")
    A[29, 23] = convert(T, big"0.")
    A[29, 24] = convert(T,
        big"-25.70181397198118326258738829725199395111365563419600781824702737091645129169813134401")
    A[29, 25] = convert(T,
        big"-.8263556911513155086442113083991534587014231586161685769224194977471882335420141183213")
    A[29, 26] = convert(T,
        big".5137008137681933419570044566186303037387573636419640300869712169933398305905931343468")
    A[29, 27] = convert(T,
        big"-.4398702297150466850587900923415450260461038902942613590425808839943205635447284745074e-1")
    A[29, 28] = convert(T,
        big"-.1375856067633252248656596321967877466474472229750848659754400903987833771639575727867")
    A[30, 1] = convert(T,
        big".1243805266540944128815164208687993162684914663596714231632892354628068537117612942798")
    A[30, 2] = convert(T, big"0.")
    A[30, 3] = convert(T, big"0.")
    A[30, 4] = convert(T, big"0.")
    A[30, 5] = convert(T,
        big".2261202821975843014222386629792029011967523207426331439651447460281196206643404356021")
    A[30, 6] = convert(T,
        big".1378858876180808806076958370164778145309694174914933853635428709475288586061552782365e-1")
    A[30, 7] = convert(T,
        big"-.6722101339966844497493995074143058569500863415253821828561997825320849038679063596730e-1")
    A[30, 8] = convert(T, big"0.")
    A[30, 9] = convert(T, big"0.")
    A[30, 10] = convert(T,
        big"-.8562389750854283547553497698795017721121215974115638028550665385850612741040225222977")
    A[30, 11] = convert(T,
        big"-1.963375228668589089282628500280938139881804405182674045535756631526916950083353845169")
    A[30, 12] = convert(T,
        big"-.2323328227241194012372462573089218472501081992304199949782180319905262045718872259601")
    A[30, 13] = convert(T, big"0.")
    A[30, 14] = convert(T,
        big"4.306607190864533494616689368765629477724325620534780926267640393608500758570100495873")
    A[30, 15] = convert(T,
        big"-2.927229632494654826597879112023904466876873949506336126307786635262992367484998786517")
    A[30, 16] = convert(T,
        big"-82.31316663978589444544923341054587077357619664281386893950601309356417181948645997040")
    A[30, 17] = convert(T, big"0.")
    A[30, 18] = convert(T, big"0.")
    A[30, 19] = convert(T, big"0.")
    A[30, 20] = convert(T, big"0.")
    A[30, 21] = convert(T, big"0.")
    A[30, 22] = convert(T, big"0.")
    A[30, 23] = convert(T, big"0.")
    A[30, 24] = convert(T,
        big"82.31316663978589444544923341054587077357619664281386893950601309356417181948645997040")
    A[30, 25] = convert(T,
        big"2.927229632494654826597879112023904466876873949506336126307786635262992367484998786517")
    A[30, 26] = convert(T,
        big"-4.306607190864533494616689368765629477724325620534780926267640393608500758570100495873")
    A[30, 27] = convert(T,
        big".2323328227241194012372462573089218472501081992304199949782180319905262045718872259601")
    A[30, 28] = convert(T,
        big"1.963375228668589089282628500280938139881804405182674045535756631526916950083353845169")
    A[30, 29] = convert(T,
        big".8562389750854283547553497698795017721121215974115638028550665385850612741040225222977")
    A[31, 1] = convert(T,
        big".1034845616366797766729935465119103444997447982019713166066629728281981965079290745983")
    A[31, 2] = convert(T, big"0.")
    A[31, 3] = convert(T, big"0.")
    A[31, 4] = convert(T,
        big".1220688873064072225896440828689620771395927148341621347412746563709055937325311521675")
    A[31, 5] = convert(T,
        big".4825744903312466224751347801256881128659190238501680496794015023696413273862321544150")
    A[31, 6] = convert(T,
        big"-.3814096000156069997308862400056202056641130724784114774219699240039767479629669855696e-1")
    A[31, 7] = convert(T, big"0.")
    A[31, 8] = convert(T,
        big"-.5504995253108023241383885070205081774114143110000375617128363206424473498745141065969")
    A[31, 9] = convert(T, big"0.")
    A[31, 10] = convert(T,
        big"-.7119158115851892278876482620437943875782918824067455704957652139710574799878630163853")
    A[31, 11] = convert(T,
        big"-.5841296056715513404329887301584808720953353296452275957070524410065417676683463009109")
    A[31, 12] = convert(T, big"0.")
    A[31, 13] = convert(T, big"0.")
    A[31, 14] = convert(T,
        big"2.110463081258649321287173000466227503003750542789369878507182287710881470618943318741")
    A[31, 15] = convert(T,
        big"-.8374947367395721355257420230010379926952601753351235177405529298334532793741463162845e-1")
    A[31, 16] = convert(T,
        big"5.100214990723209140752959690433441131075450608628042491597346388445135412965217165555")
    A[31, 17] = convert(T, big"0.")
    A[31, 18] = convert(T, big"0.")
    A[31, 19] = convert(T, big"0.")
    A[31, 20] = convert(T, big"0.")
    A[31, 21] = convert(T, big"0.")
    A[31, 22] = convert(T, big"0.")
    A[31, 23] = convert(T, big"0.")
    A[31, 24] = convert(T,
        big"-5.100214990723209140752959690433441131075450608628042491597346388445135412965217165555")
    A[31, 25] = convert(T,
        big".8374947367395721355257420230010379926952601753351235177405529298334532793741463162845e-1")
    A[31, 26] = convert(T,
        big"-2.110463081258649321287173000466227503003750542789369878507182287710881470618943318741")
    A[31, 27] = convert(T, big"0.")
    A[31, 28] = convert(T,
        big".5841296056715513404329887301584808720953353296452275957070524410065417676683463009109")
    A[31, 29] = convert(T,
        big".7119158115851892278876482620437943875782918824067455704957652139710574799878630163853")
    A[31, 30] = convert(T,
        big".5504995253108023241383885070205081774114143110000375617128363206424473498745141065969")
    A[32, 1] = convert(T,
        big".1933333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[32, 2] = convert(T, big"0.")
    A[32, 3] = convert(T, big".22")
    A[32, 4] = convert(T, big"-.8e-1")
    A[32, 5] = convert(T, big"0.")
    A[32, 6] = convert(T, big"0.")
    A[32, 7] = convert(T,
        big".1099934255807247039194624048650683408451190582958464264636524271459687549994002654752")
    A[32, 8] = convert(T,
        big"-.2542970480762701613840685069971531221418356269767039208462421656164179875269042982442")
    A[32, 9] = convert(T, big"0.")
    A[32, 10] = convert(T,
        big".8655707771166942543437703438210982818328474012330118593467368132762510892051242759318")
    A[32, 11] = convert(T,
        big"3.324164491140930831067995527865720183368600929369864071601998386039920635781409865040")
    A[32, 12] = convert(T, big"0.")
    A[32, 13] = convert(T, big"0.")
    A[32, 14] = convert(T,
        big"-12.01022233159779338823523851486618412603019426339968151272769528462035002110216728101")
    A[32, 15] = convert(T,
        big".4766014662424932394304427768620618996029637820035802094825720242694315551196576125507")
    A[32, 16] = convert(T,
        big"-29.02430112210363905258026232136540995962512213324709106915239870601916450708546744075")
    A[32, 17] = convert(T, big"0.")
    A[32, 18] = convert(T, big"0.")
    A[32, 19] = convert(T, big"0.")
    A[32, 20] = convert(T, big"0.")
    A[32, 21] = convert(T, big"0.")
    A[32, 22] = convert(T, big"0.")
    A[32, 23] = convert(T, big"0.")
    A[32, 24] = convert(T,
        big"29.02430112210363905258026232136540995962512213324709106915239870601916450708546744075")
    A[32, 25] = convert(T,
        big"-.4766014662424932394304427768620618996029637820035802094825720242694315551196576125507")
    A[32, 26] = convert(T,
        big"12.01022233159779338823523851486618412603019426339968151272769528462035002110216728101")
    A[32, 27] = convert(T, big"0.")
    A[32, 28] = convert(T,
        big"-3.324164491140930831067995527865720183368600929369864071601998386039920635781409865040")
    A[32, 29] = convert(T,
        big"-.8655707771166942543437703438210982818328474012330118593467368132762510892051242759318")
    A[32, 30] = convert(T,
        big".2542970480762701613840685069971531221418356269767039208462421656164179875269042982442")
    A[32, 31] = convert(T,
        big"-.1099934255807247039194624048650683408451190582958464264636524271459687549994002654752")
    A[33, 1] = convert(T,
        big"-.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333")
    A[33, 2] = convert(T,
        big"1.388888888888888888888888888888888888888888888888888888888888888888888888888888888889")
    A[33, 3] = convert(T, big"0.")
    A[33, 4] = convert(T, big"0.")
    A[33, 5] = convert(T, big"-.75")
    A[33, 6] = convert(T, big"0.")
    A[33, 7] = convert(T,
        big"-.4925295437180263044226820491140213202002146815806577847190740839644346370048749342561")
    A[33, 8] = convert(T, big"0.")
    A[33, 9] = convert(T, big"0.")
    A[33, 10] = convert(T, big"0.")
    A[33, 11] = convert(T, big"0.")
    A[33, 12] = convert(T, big"0.")
    A[33, 13] = convert(T, big"0.")
    A[33, 14] = convert(T, big"0.")
    A[33, 15] = convert(T, big"0.")
    A[33, 16] = convert(T, big"0.")
    A[33, 17] = convert(T, big"0.")
    A[33, 18] = convert(T, big"0.")
    A[33, 19] = convert(T, big"0.")
    A[33, 20] = convert(T, big"0.")
    A[33, 21] = convert(T, big"0.")
    A[33, 22] = convert(T, big"0.")
    A[33, 23] = convert(T, big"0.")
    A[33, 24] = convert(T, big"0.")
    A[33, 25] = convert(T, big"0.")
    A[33, 26] = convert(T, big"0.")
    A[33, 27] = convert(T, big"0.")
    A[33, 28] = convert(T, big"0.")
    A[33, 29] = convert(T, big"0.")
    A[33, 30] = convert(T, big"0.")
    A[33, 31] = convert(T,
        big".4925295437180263044226820491140213202002146815806577847190740839644346370048749342561")
    A[33, 32] = convert(T, big".75")
    A[34, 1] = convert(T,
        big".1111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
    A[34, 2] = convert(T, big"0.")
    A[34, 3] = convert(T,
        big"-.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222")
    A[34, 4] = convert(T, big"0.")
    A[34, 5] = convert(T, big"0.")
    A[34, 6] = convert(T, big"0.")
    A[34, 7] = convert(T, big"0.")
    A[34, 8] = convert(T, big"0.")
    A[34, 9] = convert(T, big"0.")
    A[34, 10] = convert(T, big"0.")
    A[34, 11] = convert(T, big"0.")
    A[34, 12] = convert(T, big"0.")
    A[34, 13] = convert(T, big"0.")
    A[34, 14] = convert(T, big"0.")
    A[34, 15] = convert(T, big"0.")
    A[34, 16] = convert(T, big"0.")
    A[34, 17] = convert(T, big"0.")
    A[34, 18] = convert(T, big"0.")
    A[34, 19] = convert(T, big"0.")
    A[34, 20] = convert(T, big"0.")
    A[34, 21] = convert(T, big"0.")
    A[34, 22] = convert(T, big"0.")
    A[34, 23] = convert(T, big"0.")
    A[34, 24] = convert(T, big"0.")
    A[34, 25] = convert(T, big"0.")
    A[34, 26] = convert(T, big"0.")
    A[34, 27] = convert(T, big"0.")
    A[34, 28] = convert(T, big"0.")
    A[34, 29] = convert(T, big"0.")
    A[34, 30] = convert(T, big"0.")
    A[34, 31] = convert(T, big"0.")
    A[34, 32] = convert(T, big"0.")
    A[34, 33] = convert(T,
        big".2222222222222222222222222222222222222222222222222222222222222222222222222222222222222")
    A[35, 1] = convert(T,
        big".2858351403889715587960888421638364148529275378945964668924322897553490152559792262023")
    A[35, 2] = convert(T,
        big".2916666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    A[35, 3] = convert(T, big".21875")
    A[35, 4] = convert(T, big"0.")
    A[35, 5] = convert(T, big".1640625")
    A[35, 6] = convert(T, big"0.")
    A[35, 7] = convert(T,
        big".2181943549455566583271882415813521070932888243221879411415164327116967439531911272777")
    A[35, 8] = convert(T,
        big".1803928984786977668636352219467754377196200536418492285624347210514163759703679527180")
    A[35, 9] = convert(T, big"0.")
    A[35, 10] = convert(T,
        big".2057138394048450188591207551229295422775700949828089053939914789386228504942804843989")
    A[35, 11] = convert(T,
        big".2427157915817702399702829279594465157627459713866705419485763522859549196625913978401")
    A[35, 12] = convert(T,
        big".2464657808136293058336092911818914077992281038693057051370210135284213379790417930740")
    A[35, 13] = convert(T,
        big"-3.449919407908908249798341546016226620603704606149316442883265523381128452524989278943")
    A[35, 14] = convert(T,
        big".2288755621600360817607290607384585842942203725527402184592948392511281334278617959957")
    A[35, 15] = convert(T,
        big".2832905997021514153215274190567333359784365954938557898314048426595070708424182066065")
    A[35, 16] = convert(T,
        big"3.210851258377666409601314905442367870055573203322387098512984999880577120008173123283")
    A[35, 17] = convert(T,
        big"-.2235387773648456999202337562141625079641252300836740320899016275445898395177373582441")
    A[35, 18] = convert(T,
        big"-.7071211572044190735187272862074872121300912319552061607910521928571247612111795934106")
    A[35, 19] = convert(T,
        big"3.211233451502870804081747292028565008932600344430223743249588034157195885590228893622")
    A[35, 20] = convert(T,
        big"1.409543483096697660304144743011231757690459455735489863573218752821178310978199657967")
    A[35, 21] = convert(T,
        big"-.1513620534437426131216022767425181110909630262036760559494590353712667648924754181285")
    A[35, 22] = convert(T,
        big".3723505745270142764547240802146199843971210282021482987373568243836683323798121465643")
    A[35, 23] = convert(T,
        big".2529787464063613367221999077621412859157757281294143192610824780367182739421617243696")
    A[35, 24] = convert(T,
        big"-3.210851258377666409601314905442367870055573203322387098512984999880577120008173123283")
    A[35, 25] = convert(T,
        big"-.2832905997021514153215274190567333359784365954938557898314048426595070708424182066065")
    A[35, 26] = convert(T,
        big"-.2288755621600360817607290607384585842942203725527402184592948392511281334278617959957")
    A[35, 27] = convert(T,
        big"-.2464657808136293058336092911818914077992281038693057051370210135284213379790417930740")
    A[35, 28] = convert(T,
        big"-.2427157915817702399702829279594465157627459713866705419485763522859549196625913978401")
    A[35, 29] = convert(T,
        big"-.2057138394048450188591207551229295422775700949828089053939914789386228504942804843989")
    A[35, 30] = convert(T,
        big"-.1803928984786977668636352219467754377196200536418492285624347210514163759703679527180")
    A[35, 31] = convert(T,
        big"-.2181943549455566583271882415813521070932888243221879411415164327116967439531911272777")
    A[35, 32] = convert(T, big"-.1640625")
    A[35, 33] = convert(T, big"-.21875")
    A[35, 34] = convert(T,
        big"-.2916666666666666666666666666666666666666666666666666666666666666666666666666666666667")
    α[1] = convert(T,
        big".1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1")
    α[2] = convert(T, big".5859375e-2")
    α[3] = convert(T, big".1171875e-1")
    α[4] = convert(T, big"0.")
    α[5] = convert(T, big".17578125e-1")
    α[6] = convert(T, big"0.")
    α[7] = convert(T, big".234375e-1")
    α[8] = convert(T, big".29296875e-1")
    α[9] = convert(T, big"0.")
    α[10] = convert(T, big".3515625e-1")
    α[11] = convert(T, big".41015625e-1")
    α[12] = convert(T, big".46875e-1")
    α[13] = convert(T, big"0.")
    α[14] = convert(T, big".52734375e-1")
    α[15] = convert(T, big".5859375e-1")
    α[16] = convert(T, big".64453125e-1")
    α[17] = convert(T, big"0.")
    α[18] = convert(T,
        big".1053521135717530196914960328878781622276730830805238840416702908213176249782427570033")
    α[19] = convert(T,
        big".1705613462417521823821203385538740858875554878027908047375010369442754416180982144816")
    α[20] = convert(T,
        big".2062293973293519407835264857011048947419142862595424540779715293772640762608018856579")
    α[21] = convert(T,
        big".2062293973293519407835264857011048947419142862595424540779715293772640762608018856579")
    α[22] = convert(T,
        big".1705613462417521823821203385538740858875554878027908047375010369442754416180982144816")
    α[23] = convert(T,
        big".1053521135717530196914960328878781622276730830805238840416702908213176249782427570033")
    α[24] = convert(T, big"-.64453125e-1")
    α[25] = convert(T, big"-.5859375e-1")
    α[26] = convert(T, big"-.52734375e-1")
    α[27] = convert(T, big"-.46875e-1")
    α[28] = convert(T, big"-.41015625e-1")
    α[29] = convert(T, big"-.3515625e-1")
    α[30] = convert(T, big"-.29296875e-1")
    α[31] = convert(T, big"-.234375e-1")
    α[32] = convert(T, big"-.17578125e-1")
    α[33] = convert(T, big"-.1171875e-1")
    α[34] = convert(T, big"-.5859375e-2")
    α[35] = convert(T,
        big".1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1")

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 14, stability_size = 1.873815353775539))
end

function constructButcher7(T::Type = Float64)
    t = sqrt(convert(T, 21))
    c2 = (7 + t) / 42
    c6 = (7 - t) / 14
    c3 = 2c2
    c = [0; c2; c3; 3c2; 1 / 2; c6; 1 / 2; 3c2; 1]
    A = zeros(T, 9, 9)
    A[2, 1] = c2
    A[3, 2] = c3
    A[4, 1:3] = [(7 + t) / 56; 0; (21 + 3t) / 56]
    A[5, 1:4] = [(8 - t) / 16; 0; (-21 + 6t) / 16; (21 - 5t) / 16]
    A[6, 1:4] = [(-1687 + 374t) / 196; 0; (969 - 210t) / 28; (-381 + 83t) / 14]
    A[6, 5] = (84 - 20t) / 49
    A[7, 1:4] = [(583 - 131t) / 128; 0; (-2373 + 501 * t) / 128; (4221 - 914 * t) / 288]
    A[7, 5:6] = [(-9 + 4t) / 18; (189 + 35t) / 576]
    A[8, 1:4] = [(-623 + 169t) / 392; 0; (435 - 81t) / 56; (-1437 + 307t) / 252]
    A[8, 5:7] = [(-2028 - 1468t) / 7497; (-21 - 4 * t) / 126; (384 + 80t) / 833]
    A[9, 1:4] = [(579 - 131t) / 24; 0; (-791 + 167 * t) / 8; (8099 - 1765t) / 108]
    A[9, 5:8] = [(-1976 + 784t) / 459; (70 + 7t) / 54; (160 - 80t) / 153; (49 - 7t) / 18]
    α = [1 // 20; 0; 0; 0; 0; 49 // 180; 16 // 45; 49 // 180; 1 // 20]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 7, stability_size = 2.6662179256363756))
end

function constructDverk(T::Type = Float64)
    c = [0; 1 // 6; 4 // 15; 2 // 3; 5 // 6; 1; 1 // 15; 1]
    A = [0 0 0 0 0 0 0 0
         1//6 0 0 0 0 0 0 0
         4//75 16//75 0 0 0 0 0 0
         5//6 -8//3 5//2 0 0 0 0 0
         -165//64 55//6 -425//64 85//96 0 0 0 0
         12//5 -8 4015//612 -11//36 88//255 0 0 0
         -8263//15000 124//75 -643//680 -81//250 2484//10625 0 0 0
         3501//1720 -300//43 297275//52632 -319//2322 24068//84065 0 3850//26703 0]
    α = [3 // 40; 0; 875 // 2244; 23 // 72; 264 // 1955; 0; 125 // 11592; 43 // 616]
    αEEst = [13 // 160; 0; 2375 // 5984; 5 // 16; 12 // 85; 3 // 44; 0; 0]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, αEEst = αEEst, adaptiveorder = 5,
        stability_size = 4.064777441244461))
end

"""
EXPLICIT RUNGE-KUTTA METHODS WITH
ESTIMATES OF THE LOCAL TRUNCATION ERROR
"""
function constructClassicVerner6(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0
         1//18 0 0 0 0 0 0 0
         -1//12 1//4 0 0 0 0 0 0
         -2//81 4//27 8//81 0 0 0 0 0
         40//33 -4//11 -56//11 54//11 0 0 0 0
         -369//73 72//73 5380//219 -12285//584 2695//1752 0 0 0
         -8716//891 656//297 39520//891 -416//11 52//27 0 0 0
         3015//256 -9//4 -4219//78 5985//128 -539//384 0 693//3328 0]
    c = [0; 1 // 18; 1 // 6; 2 // 9; 2 // 3; 1; 8 // 9; 1]
    αEEst = [3 // 80; 0; 4 // 25; 243 // 1120; 77 // 160; 73 // 700; 0; 0]
    α = [57 // 640; 0; -16 // 65; 1377 // 2240; 121 // 320; 0; 891 // 8320; 2 // 35]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 6, αEEst = αEEst, adaptiveorder = 5,
        stability_size = 4.457242953908324))
end

"""
EXPLICIT RUNGE-KUTTA METHODS WITH
ESTIMATES OF THE LOCAL TRUNCATION ERROR
"""
function constructClassicVerner7(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0 0 0
         1//12 0 0 0 0 0 0 0 0 0
         0 1//6 0 0 0 0 0 0 0 0
         1//16 0 3//16 0 0 0 0 0 0 0
         21//16 0 -81//16 9//2 0 0 0 0 0 0
         1344688//250563 0 -1709184//83521 1365632//83521 -78208//250563 0 0 0 0 0
         -559//384 0 6 -204//47 14//39 -4913//78208 0 0 0 0
         -625//224 0 12 -456//47 48//91 14739//136864 6//7 0 0 0
         -12253//99144 0 16//27 16//459 29072//161109 -2023//75816 112//12393 0 0 0
         30517//2512 0 -7296//157 268728//7379 2472//2041 -3522621//10743824 132//157 0 -12393//4396 0]
    c = [0; 1 // 12; 1 // 6; 1 // 4; 3 // 4; 16 // 17; 1 // 2; 1; 2 // 3; 1]
    αEEst = [7 // 90; 0; 0; 16 // 45; 16 // 45; 0; 2 // 15; 7 // 90; 0; 0]
    α = [2881 // 40320; 0; 0; 1216 // 2961; -2624 // 4095; 24137569 // 57482880; -4 // 21;
         0; 4131 // 3920; -157 // 1260]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 7, αEEst = αEEst, adaptiveorder = 6,
        stability_size = 5.062347991796731))
end

"""
EXPLICIT RUNGE-KUTTA METHODS WITH
ESTIMATES OF THE LOCAL TRUNCATION ERROR
"""
function constructClassicVerner8(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0 0 0 0 0 0
         1//4 0 0 0 0 0 0 0 0 0 0 0 0
         5//72 1//72 0 0 0 0 0 0 0 0 0 0 0
         1//32 0 3//32 0 0 0 0 0 0 0 0 0 0
         106//125 0 -408//125 352//125 0 0 0 0 0 0 0 0 0
         1//48 0 0 8//33 125//528 0 0 0 0 0 0 0 0
         -1263//2401 0 0 39936//26411 -64125//26411 5520//2401 0 0 0 0 0 0 0
         37//392 0 0 0 1625//9408 -2//15 61//6720 0 0 0 0 0 0
         17176//25515 0 0 -47104//25515 1325//504 -41792//25515 20237//145800 4312//6075 0 0 0 0 0
         -23834//180075 0 0 -77824//1980825 -636635//633864 254048//300125 -183//7000 8//11 -324//3773 0 0 0 0
         12733//7600 0 0 -20032//5225 456485//80256 -42599//7125 339227//912000 -1029//4180 1701//1408 5145//2432 0 0 0
         -27061//204120 0 0 40448//280665 -1353775//1197504 17662//25515 -71687//1166400 98//225 1//16 3773//11664 0 0 0
         11203//8680 0 0 -38144//11935 2354425//458304 -84046//16275 673309//1636800 4704//8525 9477//10912 -1029//992 0 729//341 0]
    c = [0; 1 // 4; 1 // 12; 1 // 8; 2 // 5; 1 // 2; 6 // 7; 1 // 7; 2 // 3; 2 // 7; 1;
         1 // 3; 1]
    αEEst = [13 // 288; 0; 0; 0; 0; 32 // 125; 31213 // 144000; 2401 // 12375;
             1701 // 14080; 2401 // 19200; 19 // 450; 0; 0]
    α = [31 // 720; 0; 0; 0; 0; 16 // 75; 16807 // 79200; 16807 // 79200; 243 // 1760; 0; 0;
         243 // 1760; 31 // 720]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.076743619212314))
end

"""
constructRKF8()

Constructs the tableau object for the Runge-Kutta-Fehlberg Order 7/8 method.
"""
function constructRKF8(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0 0 0 0 0 0
         2//27 0 0 0 0 0 0 0 0 0 0 0 0
         1//36 1//12 0 0 0 0 0 0 0 0 0 0 0
         1//24 0 1//8 0 0 0 0 0 0 0 0 0 0
         5//12 0 -25//16 25//16 0 0 0 0 0 0 0 0 0
         1//20 0 0 1//4 1//5 0 0 0 0 0 0 0 0
         -25//108 0 0 125//108 -65//27 125//54 0 0 0 0 0 0 0
         31//300 0 0 0 61//225 -2//9 13//900 0 0 0 0 0 0
         2 0 0 -53//6 704//45 -107//9 67//90 3 0 0 0 0 0
         -91//108 0 0 23//108 -976//135 311//54 -19//60 17//6 -1//12 0 0 0 0
         2383//4100 0 0 -341//164 4496//1025 -301//82 2133//4100 45//82 45//164 18//41 0 0 0
         3//205 0 0 0 0 -6//41 -3//205 -3//41 3//41 6//41 0 0 0
         -1777//4100 0 0 -341//164 4496//1025 -289//82 2193//4100 51//82 33//164 12//41 0 1 0]
    αEEst = [41 // 840; 0; 0; 0; 0; 34 // 105; 9 // 35; 9 // 35; 9 // 280; 9 // 280;
             41 // 840; 0; 0]
    α = [0; 0; 0; 0; 0; 34 // 105; 9 // 35; 9 // 35; 9 // 280; 9 // 280; 0; 41 // 840;
         41 // 840]
    c = [0; 2 // 27; 1 // 9; 1 // 6; 5 // 12; 1 // 2; 5 // 6; 1 // 6; 2 // 3; 1 // 3; 1; 0;
         1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.007588848940572))
end

"""
constructDormandPrice8_64bit()

Constructs the tableau object for the Dormand-Prince Order 6/8 method with the
approximated coefficients from the paper. This works until below 64-bit precision.

"""
function constructDormandPrince8_64bit(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0 0 0 0 0 0
         1//18 0 0 0 0 0 0 0 0 0 0 0 0
         1//48 1//16 0 0 0 0 0 0 0 0 0 0 0
         1//32 0 3//32 0 0 0 0 0 0 0 0 0 0
         5//16 0 -75//64 75//64 0 0 0 0 0 0 0 0 0
         3//80 0 0 3//16 3//20 0 0 0 0 0 0 0 0
         29443841//614563906 0 0 77736538//692538347 -28693883//1125000000 23124283//1800000000 0 0 0 0 0 0 0
         16016141//946692911 0 0 61564180//158732637 22789713//633445777 545815736//2771057229 -180193667//1043307555 0 0 0 0 0 0
         39632708//573591083 0 0 -433636366//683701615 -421739975//2616292301 100302831//723423059 790204164//839813087 800635310//3783071287 0 0 0 0 0
         246121993//1340847787 0 0 -37695042795//15268766246 -309121744//1061227803 -12992083//490766935 6005943493//2108947869 393006217//1396673457 123872331//1001029789 0 0 0 0
         -1028468189//846180014 0 0 8478235783//508512852 1311729495//1432422823 -10304129995//1701304382 -48777925059//3047939560 15336726248//1032824649 -45442868181//3398467696 3065993473//597172653 0 0 0
         185892177//718116043 0 0 -3185094517//667107341 -477755414//1098053517 -703635378//230739211 5731566787//1027545527 5232866602//850066563 -4093664535//808688257 3962137247//1805957418 65686358//487910083 0 0
         403863854//491063109 0 0 -5068492393//434740067 -411421997//543043805 652783627//914296604 11173962825//925320556 -13158990841//6184727034 3936647629//1978049680 -160528059//685178525 248638103//1413531060 0 0]
    α = [14005451 // 335480064; 0; 0; 0; 0; -59238493 // 1068277825; 181606767 // 758867731;
         561292985 // 797845732; -1041891430 // 1371343529; 760417239 // 1151165299;
         118820643 // 751138087; -528747749 // 2220607170; 1 // 4]
    αEEst = [13451932 // 455176623; 0; 0; 0; 0; -808719846 // 976000145;
             1757004468 // 5645159321; 656045339 // 265891186; -3867574721 // 1518517206;
             465885868 // 322736535; 53011238 // 667516719; 2 // 45; 0]
    c = [0; 1 // 18; 1 // 12; 1 // 8; 5 // 16; 3 // 8; 59 // 400; 93 // 200;
         5490023248 // 9719169821; 13 // 20; 1201146811 // 1299019798; 1; 1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.166633619968106))
end

"""
constructDormandPrice8()

Constructs the tableau object for the Dormand-Prince Order 6/8 method.
"""
function constructDormandPrince8(T::Type = Float64)
    A = [0 0 0 0 0 0 0 0 0 0 0 0 0
         1//18 0 0 0 0 0 0 0 0 0 0 0 0
         1//48 1//16 0 0 0 0 0 0 0 0 0 0 0
         1//32 0 3//32 0 0 0 0 0 0 0 0 0 0
         5//16 0 -75//64 75//64 0 0 0 0 0 0 0 0 0
         3//80 0 0 3//16 3//20 0 0 0 0 0 0 0 0
         215595617//4500000000 0 0 202047683//1800000000 -28693883//1125000000 23124283//1800000000 0 0 0 0 0 0 0
         14873762658037143//879168438156250000 0 0 3467633544794897//8940695981250000 1474287494383247//40978189914062500 26709270507070017//135600555715625000 -14591655588284//84484570233063 0 0 0 0 0 0
         BigInt(7586331039021946882049083502441337664277676907617750536566352)//BigInt(109794461601491217860220353338581031394059220336451160078730445) 0 0 -BigInt(236057339412812449835946465344221735535939129430991059693568)//BigInt(372184615598275314780407977418918750488336340123563254504171) -BigInt(3299739166368883603096250588167927276977533790499480498577408)//BigInt(20470153857905142312922438758040531276858498706795978997729405) BigInt(4695919603694846215470554638065271273971468502369170235542016)//BigInt(33868800019443053645017125945121606294438606951244256159879561) BigInt(291851811898394201384602939640627532330843113837053004434432000000)//BigInt(310174233778061645620360730197195350622945922304711702829528117367) BigInt(6992959981041103840944260661352231159203510904000000)//BigInt(33042342481018810238716485165383193327572243242031481) 0 0 0 0 0
         99299034813490800741867453179778547//540971123539151162906952826011200000 0 0 -2493835259080554724582//1010153717930905426875 -48550347897506146536052//166675363458599395434375 -24871192635697392099560348960246//939492072180864357472739828818125 478776089216929482237673925052922000//168119099731629344552415590032785027 6560308981643238155096750//23314158982833116227901307 BigInt(1586281686644478270321241459439899956623408540189275177)//BigInt(12818966182821619734532382093543907143647820508227904000) 0 0 0 0
         -BigInt(102116003386322998978127600084904875522141269364791505043913504184525097818434721165778087547359160299919872547571820573487921693)//BigInt(84016717385376362440519288454722754561118206109968455863629915569413007015484884082989277327146750610032897620427741658059440288) 0 0 BigInt(338590872606752219742507143357021902717271169524361004010718467428498066558752974165816979255870352236800)//BigInt(20308212073515087965058545521329962060416948491603802421256875704911573108931922671691153944392874968051) BigInt(68189290605616416787948548385820859588684790288743680764422049661712817412461535969920258787664375619072)//BigInt(74463444269555322538548000244876527554862144469213942211275210918009101399417049796200897796107208216187) -BigInt(1734282043732424474072631498514610096486835338935534079058145376760622893524503274680375038942168945756187943481380463560951840)//BigInt(286345537377499805912462279621622489249909215975695809863482134802066603511244489020404711919081949534640172065152437496912477) -BigInt(3399549280223124443696423490103003766707892326374755946138975000967466690241111348721006509128775254952212682658842765965521154240000000)//BigInt(212424385105117691648087703103838079790425456287912424851546922389328485500145214289225448961304538830766072442444722564103495915888123) BigInt(14452808190943733856347403293564049428070036006455540637351575894308889412108389906599600485253194980566957563315340127500000)//BigInt(973298753951638431793701721528200883789914680313298926814615071301495341142665245758696799918623095581715765886887649741383) -BigInt(847205714160239289113307424793539077951658318917591980262304042838612275700008766016957700930195545053374220841398660187944621107065829310608865394026418258355)//BigInt(63358704383980726998416112830322706485300332630289060627019459285960825979588560697460438306253611095891491565590971432387489415884103732012574255897878321824) BigInt(115188988949323598098458035263894669359112068207548636038131244599058496172710646622536373145562218909633738697549245770000)//BigInt(22435701423704647109276644681016984863989966659062291511947760084943925084166270812354794844590216383205333034660617626349) 0 0 0
         BigInt(21969012306961489525323859125985377266525845354279828748)//BigInt(84868015648089839210997460517819380601933600521692915045) 0 0 -2291872762438069505504//480025046760766258851 -3829018311866050387904//8800459190614048078935 -BigInt(607977714773374460437401016185253441418120832060126402968)//BigInt(199370728929424959394190105343852509479613745231838418799) BigInt(5302029233035772894614097632213626682295966947853615180783170000000)//BigInt(950538766256052885387161080614691196420735587733978871061913292363) BigInt(102968047255116137164987219663037502898143843145000000)//BigInt(16726911019578511096352500731821705820659977305290973) -BigInt(111383789341965407321602142444917514115800834690201329379027449761759895100011973929185171163615)//BigInt(22003454775272439861723739055800175619777853128055268766511800511549546753240522083740083243539) BigInt(44737471541467333111555512048686345065750)//BigInt(20391511842264262870398286145868618178341) BigInt(596546910748352988538198147432444829112451075399436970876618894337461087953328002664759407401623072330633057948252)//BigInt(4431076125983762085449284205348478790535717302043416234911901479328512794465980800998816354448181196721636373483787) 0 0
         BigInt(1066221205855832326088695778460159015192405644968016897066521076847764032613686056268693633)//BigInt(1296431693610525557488309197474904206216262654240544950471874305723890174339356551609704000) 0 0 -BigInt(1335791413506612664643690684478806471077526746614666064)//BigInt(114574907798601779179110271814903983120429559544320175) -BigInt(1591415543044168099882026495959288688569084060473110176)//BigInt(2100539976307699284950354983273239690541208591645869875) BigInt(33975758488532631832742416857645572913178866704247539610423012370193845167470455176890924)//BigInt(47586856225469573819304596274208152402640120925455970356063642741972959597009066064956075) BigInt(12176653428667113090492984656207574633063967759246601254930448409444470870786024235115138527800000)//BigInt(1008353786145118968620988891518234034224047994442049071310258686840184337101721351612973016221399) -BigInt(339784374935367314296824613776444883113869450234942131172912300100535979345925250000)//BigInt(159698690787587746004588725210359673189662237866695585709500421500486548151424426361) BigInt(4955095692700499418628052380948016677978733013841365878109775677669056866398110949788869771135857671298802131693154421086808143)//BigInt(2489789885462873158531234022579722982784822257458164105126884288597324542930882581099522281388970940826324647386340365850671680) -BigInt(563115171027780776675066866318087406247194110301648522108648094708415)//BigInt(2403532595444498372383116767918060257292523183751650851596520916634577) BigInt(147332487580158450887955957061658718012538967463083369806963200702426559434915876714751833908862217396388157664714990174448521780809)//BigInt(837599084085749358149340415048050308970085851893614803629073546048735327947816070400330404870816820234727495143522673498826476267825) 0 0]
    α = [212810988215683677989664967567559 // 5097575504458999984164528930580800; 0; 0; 0;
         0; -570667999368605802515460802224128 // 10291145812277763122885317774476825;
         BigInt(3970894643399159150754126826496000000000000) //
         BigInt(16592904867230933191457493387696939021741363)
         BigInt(177094288219480472437690862000000000000) //
         BigInt(251729356670100506734814442705774463449)
         -BigInt(66822609448295850920212176513645119787713273203022994500406050793972052314809461629969645683) //
         BigInt(87952305220338336969447643899150816363456821562985998778022435070001091778042097545895594560)
         314652731163869955629145958568800000 // 476340207420551356675670184044905167;
         BigInt(177014954088789647707522848990757432519504314686067075784476503038212450536095365316360385634933688213244039743969578872631174179769) //
         BigInt(1119019983628991838522384101261104859676427163726922121733732080377576616485631933067985100908132443862205090961383250990215178108200)
         -454665916000392064556420344242099 // 1909482158429176288068071462671400; 1 // 4]
    αEEst = [7136040226482108704342809557217 // 241464102842794736092004001974880; 0; 0; 0;
             0; -15349154422148033115423212285265536 // 18524062462099973621193571994058285;
             BigInt(45434521806506196832804182374790400000000) //
             BigInt(145978635195580057402851847985603569106229)
             BigInt(365696286946774693155766999232150000000) //
             BigInt(148214481030059176862554298041717674741)
             -BigInt(836336669851503831866889530158468123932231502753408325817124013619515886965077571) //
             BigInt(328368994730082689886153304749497093954319862912916225944630536728837081959128864)
             294694385044387823293019951454286000 // 204145803180236295718144364590673643;
             BigInt(1759482754698187564675489259591170188433054767657805212470918093603353527288272972728828708146708084742711724049636) //
             BigInt(22155380629918810427246421026742393952678586510217081174559507396642563972329904004994081772240905983608181867418935)
             2 // 45; 0]
    c = [0; 1 // 18; 1 // 12; 1 // 8; 5 // 16; 3 // 8; 59 // 400; 93 // 200;
         5490023248 // 9719169821; 13 // 20; 30992876149296355 // 33518267164510641; 1; 1]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    αEEst = map(T, αEEst)
    return (ExplicitRKTableau(A, c, α, 8, αEEst = αEEst, adaptiveorder = 7,
        stability_size = 5.166633619968108))
end

"""
  Runge Kutta Oliver method, 6 stages, 5th order
  https://www.sciencedirect.com/science/article/pii/S0096300319301511
"""

function constructRKO65(T::Type = Float64)
    A = zeros(T, 6, 6)
    c = zeros(T, 6)
    α = zeros(T, 6)

    c[1] = 2 // 3
    c[2] = 1 // 6
    c[3] = 3 // 4
    c[4] = 1 // 1
    c[5] = 4 // 5
    c[6] = 1 // 1

    α[2] = 54 // 133
    α[3] = 32 // 21
    α[4] = 1 // 18
    α[5] = -125 // 114
    α[6] = 1 // 9

    A[2, 1] = 1 // 6
    A[3, 1] = -15 // 8
    A[4, 1] = -9 // 1
    A[5, 1] = -3 // 1

    A[3, 2] = 21 // 8
    A[4, 2] = 75 // 7
    A[5, 2] = 34257 // 8750
    A[6, 2] = 123 // 380 # changed w.r.t. the paper

    A[4, 3] = -5 // 7
    A[5, 3] = -114 // 875
    A[6, 3] = 5 // 2

    A[5, 4] = 19 // 1250
    A[6, 4] = 3 // 20

    A[6, 5] = -75 // 38

    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    return (ExplicitRKTableau(A, c, α, 5, stability_size = 4.165854606804698))
end
