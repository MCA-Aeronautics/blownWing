using Pkg
Pkg.add("Revise")
Pkg.add("CSV")
Pkg.add("PyPlot")
using Revise
using CSV
using PyPlot
revise()

# Developing the blownWing package and extracting the functions
Pkg.develop(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/Code/blownWing"))
import blownWing.generateWingGeometry
import blownWing.generatePropellerWake
import blownWing.solveBlownWing
import blownWing.calculatePropellerProperties

Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
import makeAirfoil.naca

#-- Defining the wing geometry --#

# the coordinates are defined starting with the leading edge root, and then clockwise from there,
# units are in meters

coordinates = [0     0     0;
               0     0.640 0;
               0.240 0.640 0;
               0.240 0     0];

numPanelsSpan = 100;

wingGeometry = generateWingGeometry(coordinates, numPanelsSpan);
numPanels = length(wingGeometry[:,1])

#-- Defining the airfoil --#
data = CSV.read("/Users/markanderson/Box/FLOW-MCA/Code/blownWing/airfoil-data/NACA642-015A Coordinates.csv")
xcoords = convert(Array,data[1:end,1])
ycoords = convert(Array,data[1:end,2])
airfoil = cat(xcoords,ycoords,dims=2)
airfoilName = "NACA642-015A"

#-- Defining the propeller --#
airfoilData = "/Users/markanderson/Box/FLOW-MCA/Code/blownWing/Propeller Data/E212 Propeller Data.dat"
Rhub = 0.01
Rtip = 0.236 / 2
numBlades = 4 # Four blades

# Blade-wise properties
r = [0.148, 0.254237, 0.381356, 0.508475, 0.635593, 0.762712, 0.889831, 0.99].*Rtip;
chord = [9.88, 11.88, 15.59, 18.81, 19.55, 18.32, 13.96, 0.01].*Rtip*1e-3;
theta = [35.0, 32.5, 26.5, 23.5, 19, 16.5, 14.0, 10.0].*pi/180;

# Defining operating point
J = 0.85
rho = 1.225
Vinf = 50
Omega = 2*pi*Vinf / (J* 2 * Rtip)
println(Omega)

calculatePropellerProperties(airfoilData,Rhub,Rtip,numBlades,r,chord,theta,J,rho,Vinf,Omega)

#-- Defining the freestream --#
#Vinf = 10; # m/s
alpha = 4; # degrees
sideslipAngle = 0; # degrees

freestream = zeros(numPanels,3)
for i = 1:numPanels
    freestream[i,:] = [Vinf*cosd(alpha),0,Vinf*sind(alpha)]
end

#-- Defining the propeller wake --#

wakeData = CSV.read("Propeller Data/E212 Wake.csv");
propDiameter = 0.236; # meters
propPosition = 0.300; # meters

propellerWake = generatePropellerWake(wakeData,propDiameter,propPosition,wingGeometry);

wakedFreestream = freestream .+ propellerWake

#-- Running the solver --#

CL, CDi, cl, spanLocations = solveBlownWing(wingGeometry,airfoil,airfoilName,wakedFreestream)

#-- Validating the results against wind tunnel data --#

# Get the Veldhius data
data = CSV.read("/Users/markanderson/Box/FLOW-MCA/Code/blownWing/Validation Data/Veldhuis Propeller Data.csv")
spancoords = convert(Array,data[1:end,1])
unnormalizedLift = convert(Array,data[1:end,2])

figure()
plot(spanLocations./maximum(spanLocations),cl,label="Strip Theory", color = "orange", linestyle = "--", linewidth = 3)
plot(spancoords./maximum(spanLocations),unnormalizedLift,label="Veldhius Data", color = "black",marker = "o")
xlim(0,1)
title("Veldhuis CL Comparison")
legend()
xlabel("Span Location (2y/b)")
ylabel("cl")