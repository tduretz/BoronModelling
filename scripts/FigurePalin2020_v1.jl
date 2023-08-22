using CSV, DataFrames, GLMakie, Interpolations, SpecialFunctions
Makie.inline!(false)

#---------------------------------------------------------#

function MakeFigure() # Main function

    printfig = false

    # Read data into a dataframe
    df_Condie10   =  reverse(CSV.read("data/CrustGrowthCondieAlster2010.csv", DataFrame))
    df_Spencer17  =  reverse(CSV.read("data/CrustGrowthSpencer2017.csv",      DataFrame))
    df_Kalderon21 =  reverse(CSV.read("data/Global_Outgassing_Kalderon.csv",  DataFrame))
    df_Palin20    =  reverse(CSV.read("data/PalinFreeBoard.csv",              DataFrame))

    # Interpolate data on regular grid
    t = collect(LinRange(4.0, 0.0,  1000))
    t_years   = t.*1e9 
    itp_Condie10  = linear_interpolation(df_Condie10[:,1],  df_Condie10[:,2]./100,    extrapolation_bc=Line())
    itp_Spencer17 = linear_interpolation(df_Spencer17[:,1], df_Spencer17[:,2]./100,   extrapolation_bc=Line())
    itp_Kalderon21= linear_interpolation(df_Kalderon21[:,1], df_Kalderon21[:,2]./2.3, extrapolation_bc=Line())
    itp_Palin20   = linear_interpolation(df_Palin20[:,1], df_Palin20[:,2],            extrapolation_bc=Line())

    Condie10      = itp_Condie10.(t)
    Spencer17     = itp_Spencer17.(t)
    Kalderon21    = itp_Kalderon21.(t)
    Palin20       = itp_Palin20.(t)
    ConstFlux     = ones(length(t))
    LinDecay      = LinRange(1.042, 1.042, length(t))

    # Synthetic functions
    log_type1 = 1.0 .- 0.6*log.(t .+ 1)
    erf_type1 = 1.0 .- 1*erf.(t./2)
    exp_type1 = 2.0 .- (exp.(t./8 ))

    # Make 6 fluxes out of 4 digitized dataframes
    # Ocean crust: 1) alteration of oceanic crust, 2) hydrothermal input 3) accretionary prisms, see Lemarchand et al., (2002)

    # Outflux
    qAlterOcean     =  LinDecay.*27*1e10 # g Boron/y
    # qCarbonates     = ConstFlux.* 6*1e10 # g Boron/y
    qCrustSedim     = Spencer17.*13*1e10 # g Boron/y 

    # Influx
    qHydroThermal   =  Kalderon21.* 4*1e10 # g Boron/y
    # qAccrePrism     = OceanCrust.* 2*1e10 # g Boron/y
    qRunOff         = Spencer17.*    38*1e10 # g Boron/y

    # Initial condition
    Bref = 1.39e24                 # g - modulates the amplitude of the signal, not the shape
    B0   = 0.0*Bref/1e6           # ppm -> g
    Δt   = (t_years[2]-t_years[1]) # years
    nt   = length(t)
    B    = B0 * ones(nt)           # ppm initial boron concentration

    # Integrate
    for it=2:nt
        ∑sources = (qRunOff[it] + qHydroThermal[it]) #+ qAccrePrism[it] )    # g/year # Lemarchand et al., 2002
        ∑sinks   = (qAlterOcean[it]  + qCrustSedim[it]) # + qCarbonates[it]) # g/year
        B[it]    = B[it-1] + Δt*( ∑sources - ∑sinks )
    end

    # Figure
    f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
    ax1 = Axis(f[1, 1], title = L"$$Palin et al., 2020", xlabel = L"$t$ [Ga]", ylabel = L"$$Continental crust volume [%]",  xreversed = false)
    # Sampled data
    # scatter!(ax1, df_Condie10[:,1],  df_Condie10[:,2]./100,  marker=:circle, label="Condie10 raw" )
    # scatter!(ax1, df_Spencer17[:,1], df_Spencer17[:,2]./100, marker=:circle, label="Spencer17 raw")
    # scatter!(ax1, df_Kalderon21[:,1], df_Kalderon21[:,2]./2.3, marker=:circle, label="Kalderon21 raw")
    
    # Interpolated data
    lines!(ax1, t, Condie10,    label="Condie10" )
    lines!(ax1, t, Spencer17,   label="Spencer17")
    # lines!(ax1, t, LinDecay,    label="LinDecay")
    lines!(ax1, t, Palin20,     label="Palin20")
    # lines!(ax1, t, Kalderon21, label="Kalderon21")

    # Synthetic data
    # lines!(ax1, t, log_type1, label="log" )
    # lines!(ax1, t, erf_type1, label="erf" )
    # lines!(ax1, t, exp_type1, label="exp" )
    f[1, 2] = Legend(f, ax1, "Legend", framevisible = false)

    ax2 = Axis(f[2, 1], title = L"$$Boron concentration", xlabel = L"$t$ [Ga]", ylabel = L"$$B [ppm]",  xreversed = false, xgridvisible = true, ygridvisible = false)
    lines!(ax2, t, B./(Bref/1e6), label="Boron concentration" )
    # ylims!(ax2, 3.0, 6.0)

    # Display figure
    if printfig
        save("figures/ContinentalCrustPalin20.eps", f, px_per_unit = 300)
    else
        DataInspector(f)
        display(f)
    end

end

#---------------------------------------------------------#

MakeFigure() # Call the main function