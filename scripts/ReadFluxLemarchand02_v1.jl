using CSV, DataFrames, GLMakie, Interpolations, SpecialFunctions, LaTeXStrings
Makie.inline!(false)

#---------------------------------------------------------#

function MakeFigure() # Main function

    printfig = false

    # Read data into a dataframe
    df_OceanCrust = (CSV.read("data/OceanicProductionLemarchand02.csv",DataFrame)); df_OceanCrust[:,1] .-= 140.0; df_OceanCrust[:,1] .= reverse(df_OceanCrust[:,1])
    df_CrustSedim = (CSV.read("data/CrustalSedFluxLemarchand02.csv",   DataFrame)); df_CrustSedim[:,1] .-= 140.0; df_CrustSedim[:,1] .= reverse(df_CrustSedim[:,1])
    df_Carbonates = (CSV.read("data/CarbonatesLemarchand02.csv",       DataFrame)); df_Carbonates[:,1] .-= 140.0; df_Carbonates[:,1] .= reverse(df_Carbonates[:,1])
    df_Runoff     = (CSV.read("data/RunoffLemarchand02.csv",           DataFrame)); df_Runoff[:,1]     .-= 140.0; df_Runoff[:,1]     .= reverse(df_Runoff[:,1])
   
    @show df_OceanCrust[1,1]
    @show df_OceanCrust[end,1]

    # Time
    t = collect(LinRange(-140., 0.,  100))
    t_years   = t.*1e6 

    # # Interpolate data on regular grid
    itp_CrustSedim = linear_interpolation(df_CrustSedim[:,1], df_CrustSedim[:,2], extrapolation_bc=Line())
    itp_OceanCrust = linear_interpolation(df_OceanCrust[:,1], df_OceanCrust[:,2], extrapolation_bc=Line())
    itp_Carbonates = linear_interpolation(df_Carbonates[:,1], df_Carbonates[:,2], extrapolation_bc=Line())
    itp_Runoff     = linear_interpolation(df_Runoff[:,1],     df_Runoff[:,2],     extrapolation_bc=Line())
    CrustSedim     = itp_CrustSedim.(t)
    OceanCrust     = itp_OceanCrust.(t)
    Carbonates     = itp_Carbonates.(t)
    Runoff         = itp_Runoff.(t)

    # Make 6 fluxes out of 4 digitized dataframes
    # Ocean crust: 1) alteration of oceanic crust, 2) hydrothermal input 3) accretionary prisms, see Lemarchand et al., (2002)
    # Outflux
    qAlterOcean     = OceanCrust.*27*1e10 # g Boron/y
    qCarbonates     = Carbonates.*6*1e10  # g Boron/y
    qCrustSedim     = CrustSedim.*13*1e10 # g Boron/y 
    # Influx
    qHydroThermal   = OceanCrust.*4*1e10  # g Boron/y
    qAccrePrism     = OceanCrust.*2*1e10  # g Boron/y
    qRunOff         = Runoff.*38*1e10     # g Boron/y

    # Initial condition
    Bref = 1.39e24                 # g - modulates the amplitude of the signal, not the shape
    B0   = 0.5*Bref/1e6            # ppm -> g
    Δt   = (t_years[2]-t_years[1]) # years
    nt   = length(t)
    B    = B0 * ones(nt)           # ppm initial boron concentration

    # Integrate
    for it=2:nt
        ∑sources = (qHydroThermal[it] + qAccrePrism[it] + qRunOff[it])      # g/year # Lemarchand et al., 2002
        ∑sinks   = (qAlterOcean[it]   + qCarbonates[it] + qCrustSedim[it])  # g/year
        B[it]    = B[it-1] + Δt*( ∑sources - ∑sinks )
    end

    # Boron balance
    ΔB  = (qHydroThermal + qAccrePrism + qRunOff) - (qAlterOcean + qCarbonates + qCrustSedim)
    ΔB  = diff(B, dims=1) ./ B0

    # Figure
    f = Figure(resolution = (1000,1200), fontsize=25, aspect = 2.0)
    ax1 = Axis(f[1, 1], title = L"$$OceanCrust, Lemarchand et al., 2002", xlabel = L"$t$ [Ma]", ylabel = L"$$OceanCrust [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)
    # Sampled data
    scatter!(ax1, df_OceanCrust[1:1:end,1], df_OceanCrust[1:1:end,2], marker=:circle, label="OceanCrust raw")
    # Interpolated data
    lines!(ax1, t, OceanCrust, label="OceanCrust")

    ax1 = Axis(f[2, 1], title = L"$$CrustSedim, Lemarchand et al., 2002", xlabel = L"$t$ [Ma]", ylabel = L"$$CrustSedim [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)
    scatter!(ax1, df_CrustSedim[1:1:end,1], df_CrustSedim[1:1:end,2], marker=:circle, label="CrustSedim raw" )
    lines!(ax1, t, CrustSedim, label="CrustSedim" )

    ax1 = Axis(f[1, 2], title = L"$$Carbonates, Lemarchand et al., 2002", xlabel = L"$t$ [Ma]", ylabel = L"$$Carbonates [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)
    scatter!(ax1, df_Carbonates[1:1:end,1], df_Carbonates[1:1:end,2], marker=:circle, label="Carbonates raw" )
    lines!(ax1, t, Carbonates, label="Carbonates" )

    ax1 = Axis(f[2, 2], title = L"$$Runoff, Lemarchand et al., 2002", xlabel = L"$t$ [Ma]", ylabel = L"$$Runoff [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)
    scatter!(ax1, df_Runoff[1:1:end,1], df_Runoff[1:1:end,2], marker=:circle, label="Runoff raw" )
    lines!(ax1, t, Runoff, label="Runoff" )

    # f[3, 1] = Legend(f, ax1, "Legend", framevisible = false)
    ax2 = Axis(f[3, 1], title = L"$$Boron concentration", xlabel = L"$t$ [Ga]", ylabel = L"$$B [ppm]",  xreversed = false, xgridvisible = true, ygridvisible = false)
    lines!(ax2, t, B./(Bref/1e6), label="Boron concentration" )

    # f[3, 2] = Legend(f, ax1, "Legend", framevisible = false)
    ax2 = Axis(f[3, 2], title = L"$$Boron balance %", xlabel = L"$t$ [Ga]", ylabel = L"$$ ΔB [%]",  xreversed = false, xgridvisible = true, ygridvisible = false)
    lines!(ax2, t[1:end-1], ΔB.*100, label="Boron balance %" )

    # Display figure
    if printfig
        save("figures/Lemarchand02.png", f, px_per_unit = 300)
    else
        DataInspector(f)
        display(f)
    end

end

#---------------------------------------------------------#

MakeFigure() # Call the main function
