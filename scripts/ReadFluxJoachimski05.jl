using CSV, DataFrames, GLMakie, Interpolations, SpecialFunctions
Makie.inline!(false)

#---------------------------------------------------------#

function MakeFigure() # Main function

    printfig = false

    # Read data into a dataframe
    df_CrustSedim = (CSV.read("data/CrustalSedimentFluxJoachimski05.csv", DataFrame))
    df_OceanCrust = (CSV.read("data/OceanCrustProductionJoachimski05.csv",      DataFrame))
    df_CrustSedim[:,2] .*= 2.
    # # Interpolate data on regular grid
    t = collect(LinRange(-550., -150.,  100))
    t_years   = t.*1e6 
    t_years .-= t_years[1]
    itp_CrustSedim = linear_interpolation(df_CrustSedim[:,1], df_CrustSedim[:,2], extrapolation_bc=Line())
    itp_OceanCrust = linear_interpolation(df_OceanCrust[:,1], df_OceanCrust[:,2], extrapolation_bc=Line())
    CrustSedim     = itp_CrustSedim.(t)
    OceanCrust     = itp_OceanCrust.(t)

    Bref = 1.39e24           # g
    @show B0 = 3.5*Bref/1e6 # g
    @show Δt   = (t_years[2]-t_years[1])   # years

    nt = length(t)
    B  = B0 * ones(nt)    # ppm initial boron concentration
 

    # B  = 1e10.*((CrustSedim*38 .+ 4 .+ 2 ) .- (OceanCrust*27 .+ 13 .+ 6) ).*t_years .+ B0

    # for it=2:nt
    #     B[it] = 1e10*((CrustSedim[it]*38 + 4 + 2 ) - (OceanCrust[it]*27 + 13 + 6) )*t_years[it] .+ B0
    # end

    for it=2:nt
        sources = (CrustSedim[it]*38 + 6 + 0*30 )*1e10  # g/year # Lemarchand et al., 2002
        sinks   = (OceanCrust[it]*27 + 13 + 6)*1e10  # g/year
        B[it]   = B[it-1]  + Δt*( sources - sinks )
    end

    # Figure
    f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
    ax1 = Axis(f[1, 1], title = L"$$Joachimski et al., 2005", xlabel = L"$t$ [Ga]", ylabel = L"$$Continental crust volume [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)
    # Sampled data
    scatter!(ax1, df_CrustSedim[1:20:end,1], df_CrustSedim[1:20:end,2], marker=:circle, label="CrustSedim raw" )
    scatter!(ax1, df_OceanCrust[1:20:end,1], df_OceanCrust[1:20:end,2], marker=:circle, label="OceanCrust raw")
    # Interpolated data
    lines!(ax1, t, CrustSedim, label="CrustSedim" )
    lines!(ax1, t, OceanCrust, label="OceanCrust")
    # lines!(ax1, t, CrustSedim.-OceanCrust, label="OceanCrust")
    f[1, 2] = Legend(f, ax1, "Legend", framevisible = false)
    ax2 = Axis(f[2, 1], title = L"$$Joachimski et al., 2005", xlabel = L"$t$ [Ga]", ylabel = L"$$B [ppm]",  xreversed = false, xgridvisible = false, ygridvisible = false)
    lines!(ax2, t, B./(Bref/1e6), label="Boron concentration" )
    # # Display figure
    # if printfig
    #     save("figures/ContinentalCrust.png", f, px_per_unit = 300)
    # else
    # hidespines!(ax2)
    # hidedecorations!(ax2)
        DataInspector(f)
        display(f)
    # end

end

#---------------------------------------------------------#

MakeFigure() # Call the main function
