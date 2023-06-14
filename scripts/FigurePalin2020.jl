using CSV, DataFrames, GLMakie, Interpolations, SpecialFunctions
Makie.inline!(false)

#---------------------------------------------------------#

function MakeFigure() # Main function

    printfig = false

    # Read data into a dataframe
    df_Condie10  = reverse(CSV.read("data/CrustGrowthCondieAlster2010.csv", DataFrame))
    df_Spencer17 = reverse(CSV.read("data/CrustGrowthSpencer2017.csv",      DataFrame))

    # Interpolate data on regular grid
    t = LinRange(0., 4.5,  1000)
    itp_Condie10  = linear_interpolation(df_Condie10[:,1],  df_Condie10[:,2],  extrapolation_bc=Line())
    itp_Spencer17 = linear_interpolation(df_Spencer17[:,1], df_Spencer17[:,2], extrapolation_bc=Line())
    Condie10      = itp_Condie10.(t)
    Spencer17     = itp_Spencer17.(t)

    # Synthetic functions
    log_type1 = 100.0 .- 60*log.(t .+ 1)
    erf_type1 = 100.0 .- 100*erf.(t./2)
    exp_type1 = 100.0 .- (exp.(t ) .- 1.)

    # Figure
    f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
    ax = Axis(f[1, 1], title = L"$$Palin et al., 2020", xlabel = L"$t$ [Ga]", ylabel = L"$$Continental crust volume [%]",  xreversed = false)
    # Sampled data
    scatter!(ax, df_Condie10[:,1],  df_Condie10[:,2],  marker=:circle, label="Condie10 raw" )
    scatter!(ax, df_Spencer17[:,1], df_Spencer17[:,2], marker=:circle, label="Spencer17 raw")
    # Interpolated data
    lines!(ax, t, Condie10, label="Condie10" )
    lines!(ax, t, Spencer17, label="Spencer17")
    # Synthetic data
    lines!(ax, t, log_type1, label="log" )
    lines!(ax, t, erf_type1, label="erf" )
    lines!(ax, t, exp_type1, label="exp" )
    f[1, 2] = Legend(f, ax, "Legend", framevisible = false)
    # Display figure
    if printfig
        save("figures/ContinentalCrust.png", f, px_per_unit = 300)
    else
        DataInspector(f)
        display(f)
    end

end

#---------------------------------------------------------#

MakeFigure() # Call the main function