module blownWing

    using Pkg
    Pkg.add("Revise"); using Revise; revise()
    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
    import VLMMCA.generatePanels

    Pkg.develop(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/Code/NonlinearLiftingLine/"))
    import NonlinearLiftingLine.NLL

    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))

    Pkg.develop(PackageSpec(url="https://github.com/byuflowlab/CCBlade.jl"))
    Pkg.add("CCBlade")
    using CCBlade

    Pkg.add("FLOWMath")
    using FLOWMath

    Pkg.add("CSV")
    using CSV

    Pkg.add("Tables")
    using Tables

    Pkg.add("PyPlot")
    using PyPlot
    

    function solveBlownWing(panels,
                   airfoil,
                   airfoilName,
                   freestream)

        CL, CDi, cl, spanLocations = NLL(panels,airfoil,airfoilName,freestream);

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
