"""
`Base.length(tab::ODERKTableau)`

Defines the length of a Runge-Kutta method to be the number of stages.
"""
Base.length(tab::ODERKTableau) = tab.stages

"""
`stability_region(z,tab::ODERKTableau)`

Calculates the stability function from the tableau at `z`. Stable if <1.

```math
r(z) = \\frac{\\det(I-zA+zeb^T)}{\\det(I-zA)}
```
"""
stability_region(z,tab::ODERKTableau) = det(eye(tab.stages)- z*tab.A + z*ones(tab.stages)*tab.α')/det(eye(tab.stages)-z*tab.A)
