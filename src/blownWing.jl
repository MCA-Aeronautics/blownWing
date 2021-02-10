module blownWing

    using Revise; revise()
    #Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
    import VLMMCA.generatePanels

    #Pkg.develop(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/Code/NonlinearLiftingLine/"))
    import NonlinearLiftingLine.NLL

    #Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))

    #Pkg.develop(PackageSpec(url="https://github.com/byuflowlab/CCBlade.jl"))
    #Pkg.add("CCBlade")
    using CCBlade

    #Pkg.add("FLOWMath")
    using FLOWMath

    #Pkg.add("CSV")
    using CSV

    #Pkg.add("Tables")
    using Tables

    #Pkg.add("PyPlot")
    using PyPlot

    # Eduardo's stuff
    import AirfoilPrep
    ap = AirfoilPrep

    import AirfoilDatabase
    adb = AirfoilDatabase

    import JuliaDB
    jdb = JuliaDB

    #Pkg.add("Statistics")
    #Pkg.add("IterativeSolvers")
    import Statistics
    import IterativeSolvers

    #Pkg.add("DataFrames")
    using DataFrames

    function createRBFS(database_path)

        # Some of Eduardo's code
        # database_path = "/Users/markanderson/Box/FLOW-MCA/Code/NonlinearLiftingLine/airfoil-data/eduardo";
        db  = jdb.loadtable(joinpath(database_path, "index.csv"))

        # Singular kernel
        zeta_sing(r) = r==0 ? 1 : 0

        # erf Gaussian kernel
        zeta_gauserf(r) = 1/(2*pi)^(3/2) * exp(-r^2/2)

        # Gaussian kernel
        zeta_gaus(r) = 3/(4*pi)*exp(-r^3)

        # Winckelmans algebraic kernel
        zeta_wnklmns(r) = 7.5 / (4*pi) / (r^2 + 1)^3.5

        #--MAGIC CODE--#

        # Generate RBF interpolation function
        rbf_axes = [:alpha, :re, :ma, :ncrit]    # Axes of the RBF


        # Read all files in database
        Xps = Dict{Symbol, Array{Array{Float64, 1}, 1}}()
        vals = Dict{Symbol, Array{Float64, 1}}()
        rbfs = Dict{Symbol, Function}()

        # This section just extracts all of the data from the files and feeds it into the generate_RBF() function

        # for file in [:clfile, :cdfile, :cmfile]  # Iterate over each file
        for file in [:clfile]
            
            Xp = Array{Float64, 1}[]
            val = Float64[]
            
            # Convert file symbol to column header
            colheader = Symbol(adb.FIELD2_HEADER[file])
            
            Xmin = [Inf for xi in 1:length(rbf_axes)]
            Xmax = [-Inf for xi in 1:length(rbf_axes)]
            valmin, valmax = Inf, -Inf
            
            for (rowi, row) in enumerate(db) # Iterate over each row

                filename = row[colheader]   # File to read

                # Read data in file
                # println(joinpath(database_path, adb.DIRS[file], filename))
                # You may need to run `using CSV; using DataFrames;` in the REPL before this next line will work
                data = CSV.read(joinpath(database_path, adb.DIRS[file], filename),DataFrame)
                
                for drow in eachrow(data) # Iterate over rows in the data
                    this_Xp = Float64[ax != :alpha ? row[Symbol(adb.FIELD2_HEADER[ax])] : drow[1] for ax in rbf_axes]
                    this_val = drow[2]
                    push!(Xp, this_Xp)
                    push!(val, this_val)
                    
                    for xi in 1:length(rbf_axes)
                        this_Xp[xi] < Xmin[xi] ? Xmin[xi]=this_Xp[xi] : 
                        this_Xp[xi] > Xmax[xi] ? Xmax[xi]=this_Xp[xi] :
                                                nothing
                    end
                    this_val < valmin ? valmin = this_val : this_val > valmax ? valmax = this_val : nothing
                end
                
            end
            
            
            println("Generating $file RBF function with $(length(Xp)) data points...")
            
            # Scale each variable in the range 0 to 1
            X_scaling = [x == 0 ? 1 : x for x in Xmax .- Xmin]
            Xp_scaled = [(X .- Xmin) ./ X_scaling for X in Xp]
            val_scaled = (val .- valmin) ./ (valmax - valmin)
            
            # Generate RBF interpolation function
            rbf, A = generate_RBF(Xp, val; zeta=zeta_gaus, sigmas=0.1)
            
            # Generate scaled RBF interpolation function
        #     rbf_scaled, A = generate_RBF(Xp_scaled, val_scaled; zeta=zeta_gaus, sigmas=1.50)
            rbf_scaled, A = generate_RBF(Xp_scaled, val_scaled; zeta=zeta_gaus, sigmas=5.0)
            rbf(X) = valmin + (valmax - valmin)*rbf_scaled( (X.-Xmin)./ X_scaling )
            
            
            Xps[file] = Xp
            vals[file] = val
            rbfs[file] = rbf
        end

        return rbfs

    end

    # Returns a radial basis function interpolation of a field
    # with values `val` at positions `Xp`. `zeta` is the chosen
    # basis function
    function generate_RBF(Xp, val; zeta=zeta_gaus, sigmas=0.1)
        
        # ERROR CASES
        if size(Xp,1)!=size(val,1)
            error("size(Xp,1)!=size(val,1)")
        end
        
        Np = size(Xp, 1)                     # Number of data points
        
        if size(sigmas)==()
            sgms = sigmas*ones(Np)           # Spreading of every basis function
        else
            sgms = sigmas
        end
            
        # j-th scaled basis evaluated at X
        zetasgm(j, X) = zeta(Statistics.norm(X-Xp[j])/sgms[j])/sgms[j]^3
        
        # Matrix with basis functions evaluated at every point
        # Z[i,j] corresponds to the j-th basis evaluated at i-th point
        Z = [zetasgm(j, Xi) for Xi in Xp, j in 1:Np]
        
        # Solves for the alpha coefficients of every basis
        # A = Z\val
        # A = LinearAlgebra.pinv(Z)*val
        A = IterativeSolvers.gmres(Z, val)
        
        # Generates RBF interpolation function
        rbf(X) = sum([A[j]*zetasgm(j, X) for j in 1:Np])
        
        return rbf, A
    end
    

    function solveBlownWing(panels,
                   airfoil,
                   airfoilName,
                   freestream,
                   rbfs)

        CL, CDi, cl, spanLocations = NLL(panels,airfoil,airfoilName,freestream,rbfs);

        return CL, CDi, cl, spanLocations

    end # solve

    function generateWingGeometry(coordinates,numPanelsSpan)

        # "coordinates" is used to define the starboard side of the wing, starting with the root point
        # on the leading edge and moving clockwise from there

        firstCoordinate = coordinates[1,:]
        secondCoordinate = coordinates[2,:]
        thirdCoordinate = coordinates[3,:]
        fourthCoordinate = coordinates[4,:]
        
        wingGeometry = generatePanels(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,numPanelsSpan,1)

        return wingGeometry

    end # generateWingGeometry

    function generatePropellerWake(wakeData,propDiameter,propPosition,wingGeometry)

        radialPosition = wakeData[:,1]
        axialVelocity = wakeData[:,2]
        swirlVelocity = wakeData[:,3]

        panels = wingGeometry
        numPanels = length(panels[:,1])

        propellerWake = zeros(numPanels,3)

        # Determine whether each panel is within the wake. If it is, then add the axial and swirl components to the freestream
        for i = 1:numPanels

            # determine distance to center of Propeller
            distanceFromProp = abs(panels[i,2]) - propPosition

            if abs(distanceFromProp) <= propDiameter/2 # See if it's within the wake
                axial = akima(radialPosition,axialVelocity,abs(distanceFromProp))
                swirl = akima(radialPosition,swirlVelocity,abs(distanceFromProp))

                # If you are to the outside of the propeller, reverse the swirl because the propeller is
                # swinging down
                if distanceFromProp > 0
                    swirl = -swirl
                end

                propellerWake[i,:] = [axial, 0, swirl]

            end

        end

        return propellerWake

    end # generatePropeller

    function calculatePropellerProperties(airfoilData,Rhub,Rtip,numBlades,r,chord,theta,J,rho,Vinf,Omega)

        # # Temporary hard-coding
        # # Defining large-scale properties
        # Rhub = 0.01
        # Rtip = 0.236 / 2
        # B = 4 # Four blades
        # turbine = false

        # # Defining blade properties
        # r = [0.148, 0.254237, 0.381356, 0.508475, 0.635593, 0.762712, 0.889831, 0.99].*Rtip;
        # chord = [9.88, 11.88, 15.59, 18.81, 19.55, 18.32, 13.96, 0.01].*Rtip*1e-3;
        # theta = [35.0, 32.5, 26.5, 23.5, 19, 16.5, 14.0, 10.0].*pi/180;

        # af = af_from_files("Propeller Data/E212 Propeller Data.dat")
        # airfoils = fill(af,length(r))

        # rotor = Rotor(Rhub, Rtip, B, turbine)
        # sections = Section.(r, chord, theta, airfoils);

        # # Defining operating point
        # J = 0.85
        # rho = 1.225
        # Vinf = 10
        # rev_per_second = (J * 2 * Rtip) / Vinf
        # Omega = 8000 * pi / 30

        # op = simple_op.(Vinf, Omega, r, rho);

        # # Getting and plotting the outputs
        # outputs = solve.(Ref(rotor), sections ,op);

        # # End of temporary hard-coding

        af = af_from_files(airfoilData)
        airfoils = fill(af,length(r))

        rotor = Rotor(Rhub, Rtip, numBlades, false)
        sections = Section.(r, chord, theta, airfoils);

        op = simple_op.(Vinf, Omega, r, rho);

        # Getting and plotting the outputs
        outputs = solve.(Ref(rotor), sections ,op)

        data = cat(r,outputs.u,outputs.v,dims=2)
        newData = zeros(length(data[:,1]) + 2,3)

        for i = 1:(length(data[:,1]) + 2)
            if i == 1
                newData[i,:] = [0,0,0]
            elseif i == length(data[:,1]) + 2
                newData[i,:] = [Rtip, 0, 0]
            else
            newData[i,:] = data[i-1,:] 
            end
        end

        data = Tables.table(newData)

        CSV.write("/Users/markanderson/Box/FLOW-MCA/Code/blownWing/Propeller Data/E212 Wake.csv",data)

        # Plotting the wake velocity
        figure(10)
        plot(r/Rtip, outputs.u/Vinf)
        plot(r/Rtip, outputs.v/Vinf)
        title("Wake Velocity")
        xlabel("r/Rtip")
        ylabel("(normalized) induced velocity at rotor disk")
        legend(["axial velocity", "swirl velocity"]);
        
        T, Q = thrusttorque(rotor, sections, outputs)
        
        eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor)
        println("Efficiency = ", eff)
        println("Thrust Coefficient = ", CT)



    end # calculatePropellerProperties

end # module
