# Visualization

AtomTwin provides animation support for 2D atomic trajectories via a
GLMakie extension. The extension is loaded automatically when GLMakie is
available in the environment.

## Setup

Install GLMakie and load it before `using AtomTwin`:

```julia
# In your Julia environment:
# pkg> add GLMakie

using GLMakie
using AtomTwin
using AtomTwin.Visualization: animate
```

## Animation

`animate` reads position data recorded by motion detectors and produces an
interactive 2D scatter animation showing how atoms move through the tweezer
array over time. It can also save a GIF.

The typical workflow is:

1. Add `MotionDetectorSpec` entries to your system before running.
2. Call `play` to get the output.
3. Call `animate` on the output, specifying display options per detector group.

```julia
# Add a motion detector for a group of atoms
add_detector!(system, MotionDetectorSpec(tweezers; name = "row1"))

out = play(system, seq; initial_state = initial_state)

# Animate interactively
animate(out;
    options    = Dict("row1" => (color = :blue, markersize = 10)),
    limits     = ((-50, 50), (-20, 20)),   # μm
    unit       = 1e-6,                      # convert m → μm
    sleep_time = 0.001,
)

# Or save to a GIF
animate(out;
    options   = Dict("row1" => (color = :blue,)),
    gifname   = "sorting.gif",
    framerate = 30,
)
```

See the [Atom sorting example](atom_sorting.md) for a complete worked example
including motion detectors and visualization.

```@docs
AtomTwin.Visualization.animate
```
