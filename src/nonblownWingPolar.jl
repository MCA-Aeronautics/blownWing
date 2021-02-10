using Pkg
# Pkg.add("Revise")
# Pkg.add("CSV")
# Pkg.add("PyPlot")
using Revise
using CSV
using PyPlot
using Tables
using Plots
revise()

# Developing the blownWing package and extracting the functions
# Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/NonlinearLiftingLine"))
import blownWing.generateWingGeometry
import blownWing.generatePropellerWake
import blownWing.solveBlownWing
import blownWing.calculatePropellerProperties
import blownWing.createRBFS

# Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
import VLMMCA.VLM

# Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
import makeAirfoil.naca

#-- Defining the wing geometry --#

# the coordinates are defined starting with the leading edge root, and then clockwise from there,
# units are in meters

coordinates = [0     0     0;
               0     0.500 0;
               0.083 0.500 0;
               0.083 0     0];

numPanelsSpan = 50;

wingGeometry = generateWingGeometry(coordinates, numPanelsSpan);
numPanels = length(wingGeometry[:,1])

# Generating the interpolation functions
# This next line can be commented out after it has run once in the current session
rbfs = createRBFS("/Users/markanderson/.julia/dev/blownWing/airfoil-data/NACA23018/");

airfoilName = "NACA23018";
airfoil = "Filler";

angles = 0:5:40 # angles of attack for polar
CL_VLM = zeros(length(angles),1)
CL_VLM_Waked = similar(CL_VLM)
CL = similar(CL_VLM)
CL_Waked = similar(CL_VLM)
for i = 1:length(angles)

    #-- Defining the freestream --#
    Vinf = 10; # m/s
    alpha = angles[i] # degrees
    sideslipAngle = 0; # degrees

    freestream = zeros(numPanels,3)
    for j = 1:numPanels
        freestream[j,:] = [Vinf*cosd(alpha),0,Vinf*sind(alpha)]
    end

    #-- Running the Linear solver --#

    CL_VLM[i], CDi_near_VLM, cl_VLM, cd_near_VLM, spanLocations, GammaValues_VLM = VLM(wingGeometry,freestream)

    #-- Running the Nonlinear solver --#

    println("Angle of Attack: ",alpha," Degrees")

    CL[i], CDi, cl, spanLocations = solveBlownWing(wingGeometry,airfoil,airfoilName,freestream,rbfs)

    #-- Saving the results to a CSV File --#

    filename = string("/Users/markanderson/.julia/dev/blownWing/results/",angles[i],".csv")
    
    data = cat(dims=2,spanLocations,cl,cl_VLM);

    data = Tables.table(data)

    #touch(filename)

    CSV.write(filename,data)

end

# figure()
# plot(angles, CL_VLM, label = "VLM No Props", color = "green", linestyle = "--")
# plot(angles, CL, label = "Strip Theory No Props", color = "orange", linestyle = "--")
# xlabel("Angle of Attack (Degrees)")
# ylabel("CL")
# title("Veldhuis Wing Simulation")
# legend()

#-- Saving the results to a CSV File --#

filename = string("/Users/markanderson/.julia/dev/blownWing/results/Polar.csv")
    
data = cat(dims=2,angles,CL,CL_VLM);

data = Tables.table(data)

#touch(filename)

CSV.write(filename,data)

Plots.plot(angles,CL,title = "CL vs AOA",legend = false,ylabel = "Lift Coefficient",xlabel = "Angle of Attack (Degrees)")