module blownWing

    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/NonlinearLiftingLine"))
    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
    Pkg.develop(PackageSpec(url="https://github.com/byuflowlab/CCBlade.jl"))

    import VLMMCA.generatePanels
    

    function solve()

    end

    function generateWingGeometry(coordinates,numPanelsSpan)

        # "coordinates" is used to define the starboard side of the wing, starting with the root point
        # on the leading edge and moving clockwise from there

        firstCoordinate = coordinates[1,:]
        secondCoordinate = coordinates[2,:]
        thirdCoordinate = coordinates[3,:]
        fourthCoordinate = coordinates[4,:]
        
        wingGeometry = VLMMCA.generatePanels(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,numPanelsSpan,1)

        return wingGeometry

    end

    function generatePropeller()

    end

end # module
