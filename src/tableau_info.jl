"""
`Base.length(tab::ODERKTableau)`

Defines the length of a Runge-Kutta method to be the number of stages.
"""
Base.length(tab::ODERKTableau) = tab.stages

"""
`stability_region(z,tab::ODERKTableau)`

Calculates the stability function from the tableau at `z`. Stable if <1.

```math
r(z) = 1 + z bᵀ(I - zA)⁻¹ e
```
where e denotes a vector of ones.
"""
function stability_region(z, tab::ODERKTableau; embedded = false)
    A, c = tab.A, tab.c
    b = embedded ? tab.αEEst : tab.α
    e = ones(eltype(A), length(b))
    stages = (I - z * A) \ e
    1 + z * (transpose(b) * stages)
end

"""
`stability_region(tab::ODERKTableau; initial_guess=-3.0)`

Calculates the length of the stability region in the real axis.
"""
function stability_region(tab::ODERKTableau; initial_guess = -3.0, kw...)
    residual! = function (resid, x)
        resid[1] = abs(stability_region(x[1], tab)) - 1
    end
    sol = nlsolve(residual!, [initial_guess]; kw...)
    sol.zero[1]
end

function RootedTrees.residual_order_condition(tab::ODERKTableau, order::Int,
        reducer = nothing, mapper = x -> x^2;
        embedded = false)
    A, c = tab.A, tab.c
    b = embedded ? tab.αEEst : tab.α
    if reducer === nothing
        resid = map(RootedTreeIterator(order)) do t
            residual_order_condition(t, A, b, c)
        end
    else
        resid = mapreduce(reducer, RootedTreeIterator(order)) do t
            mapper(residual_order_condition(t, RungeKuttaMethod(A, b, c)))
        end
    end
    return resid
end

isfsal(tab::ExplicitRKTableau) = tab.fsal
isfsal(::ImplicitRKTableau) = nothing
function check_tableau(tab; tol = 10eps(1.0))
    order = all(i -> residual_order_condition(tab, i, +, abs) < tol, 1:(tab.order))
    if !order
        error("Tableau's order is not correct.")
    end
    if tab.adaptiveorder != 0
        embedded_order = all(
            i -> residual_order_condition(tab, i, +, abs;
                embedded = true) < tol,
            tab.adaptiveorder)
        if !embedded_order
            error("Tableau's embedded order is not correct.")
        end
    end
    return true
end
