module AtomTwinPlotsExt

using AtomTwin
using AtomTwin.Visualization
using Plots

@recipe function f(curve::PolarizabilityCurve)
    λs_main = collect(curve.λ_main)

    # Compute ylabel
    ylabel_txt = curve.unit == :au ? "Polarizability α (a.u.)" : "Δν / I  [Hz/(W/cm²)]"

    # Main plot settings
    xlabel --> "Wavelength λ (nm)"
    ylabel --> ylabel_txt
    legend --> :bottomright
    linewidth --> 2
    framestyle --> :box
    minorgrid --> true
    size --> (900, 500)
    left_margin --> 5Plots.mm
    bottom_margin --> 5Plots.mm
    right_margin --> 3Plots.mm
    top_margin --> 3Plots.mm

    # Plot main curves on subplot 1
    for model in curve.models
        if curve.unit == :au
            vals_main = [polarizability_au(model, λ) for λ in λs_main]
        else
            vals_main = [light_shift_coeff_Hz_per_Wcm2(model, λ) for λ in λs_main]
        end

        @series begin
            subplot --> 1
            ylims --> curve.ylim_main
            label --> model.state
            λs_main, vals_main
        end
    end

    # Inset zoom (if λ_inset is provided)
    if curve.λ_inset !== nothing
        λs_inset = collect(curve.λ_inset)
        
        inset_subplots --> (1, bbox(curve.inset_position...))

        for model in curve.models
            if curve.unit == :au
                vals_inset = [polarizability_au(model, λ) for λ in λs_inset]
            else
                vals_inset = [light_shift_coeff_Hz_per_Wcm2(model, λ) for λ in λs_inset]
            end

            @series begin
                subplot --> 2
                ylims --> curve.ylim_inset
                label --> ""
                guidefontsize --> 10
                legend --> false
                framestyle --> :box
                linewidth --> 2
                λs_inset, vals_inset
            end
        end
    end
end



end #module
