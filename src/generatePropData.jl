# Installing Packages
using Pkg
Pkg.add("CCBlade")
Pkg.add("PyPlot")
Pkg.add("CSV")
Pkg.add("Tables")
using CCBlade
using PyPlot
using CSV
using Tables

# Rhub and Rtip must be in meters

function generatePropData(propName,Rhub = 0.0254*0.5,Rtip = 0.0254*3.8,rho = 1.225,Vinf = 10.0,Omega = 8000*pi/30,numBlades = 2)

    turbine = false

    # Defining blade properties
    bladeProperties = CSV.read(string("Propeller Data/",propName,"/Blade Properties.csv"))
    r = convert(Array,bladeProperties[1:end,1]).*Rtip
    chord = convert(Array,bladeProperties[1:end,2]).*Rtip
    theta = convert(Array,bladeProperties[1:end,3]).*pi/180

    af = af_from_files(string("Propeller Data/",propName,"/CL vs Alpha.dat"))   
    airfoils = fill(af,length(r))

    rotor = Rotor(Rhub, Rtip, numBlades, turbine)
    sections = Section.(r, chord, theta, airfoils);

    op = simple_op.(Vinf, Omega, r, rho);

    # Getting and plotting the outputs
    outputs = solve.(Ref(rotor), sections ,op);

    data = cat(dims=2,outputs.u,outputs.v)

    data = Tables.table(data)

    touch("Propeller Data/Example/Wake Properties.csv")

    CSV.write("Propeller Data/Example/Wake Properties.csv",data)

end # generatePropData