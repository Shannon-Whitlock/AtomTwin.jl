module Visualization

export animate

"""
    animate(args...; kwargs...)

Animate atomic configurations.

This method requires GLMakie.
"""
function animate(args...; kwargs...)
    error(
        "AtomTwin.Visualization.animate requires GLMakie.\n" *
        "Install it with: pkg> add GLMakie and add `using GLMakie` in your script before calling `using AtomTwin.Visualization: animate`
"
    )
end

end