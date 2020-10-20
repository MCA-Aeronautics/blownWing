using Pkg
Pkg.add("Revise")
Pkg.add("CSV")
Pkg.add("PyPlot")
using Revise
using CSV
using PyPlot
revise()

# Developing the blownWing package and extracting the functions
Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/NonlinearLiftingLine"))
import blownWing.generateWingGeometry
import blownWing.generatePropellerWake
import blownWing.solveBlownWing
import blownWing.calculatePropellerProperties
import blownWing.createRBFS

Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
import VLMMCA.VLM

Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
import makeAirfoil.naca

#-- Defining the wing geometry --#

# the coordinates are defined starting with the leading edge root, and then clockwise from there,
# units are in meters

coordinates = [0     0     0;
               0     0.640 0;
               0.240 0.640 0;
               0.240 0     0];

numPanelsSpan = 200;

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
r_coarse = [0.148, 0.254237, 0.381356, 0.508475, 0.635593, 0.762712, 0.889831, 0.99]*Rtip;
chord_coarse = [9.88, 11.88, 15.59, 18.81, 19.55, 18.32, 13.96, 0.01].*1e-3;
theta_coarse = [35.0, 32.5, 26.5, 23.5, 19, 16.5, 14.0, 10.0].*pi/180;

#-- Interpolating the blade-wise properties --#
r = range(r_coarse[1],r_coarse[end],length = 100)
chord = zeros(length(r),1)
theta = similar(chord)

for i = 1:length(r)
    
    chord[i] = akima(r_coarse,chord_coarse,r[i])
    theta[i] = akima(r_coarse,theta_coarse,r[i])
    
end

# Defining operating point
J = 0.85/2;
J = 0.25
rho = 1.225;
D = 2*Rtip;
Vinf = 10;
n = Vinf / (J*D);
Omega = 2*pi*n;

calculatePropellerProperties(airfoilData,Rhub,Rtip,numBlades,r,chord,theta,J,rho,Vinf,Omega)

# Generating the interpolation functions
rbfs = createRBFS("/Users/markanderson/Box/FLOW-MCA/Code/NonlinearLiftingLine/airfoil-data/eduardo";)

angles = 0:2:30 # angles of attack for polar
CL_VLM = zeros(length(angles),1)
CL_VLM_Waked = similar(CL_VLM)
CL = similar(CL_VLM)
CL_Waked = similar(CL_VLM)
for i = 1:length(angles)

    #-- Defining the freestream --#
    #Vinf = 10; # m/s
    alpha = angles[i] # degrees
    sideslipAngle = 0; # degrees

    freestream = zeros(numPanels,3)
    for j = 1:numPanels
        freestream[j,:] = [Vinf*cosd(alpha),0,Vinf*sind(alpha)]
    end

    #-- Defining the propeller wake --#

    wakeData = CSV.read("Propeller Data/E212 Wake.csv");
    propDiameter = 0.236; # meters
    propPosition = 0.300; # meters

    propellerWake = generatePropellerWake(wakeData,propDiameter,propPosition,wingGeometry);

    # # Tilting the propellerWake according to the angle of attack (not necessary, the wing is flat in the 3D space)
    # for j = 1:length(propellerWake[:,1])
    #     propellerWake[j,:] = [propellerWake[j,1] * cosd(alpha), 0, propellerWake[j,3] * sind(alpha)]
    # end

    wakedFreestream = freestream .+ propellerWake;

    #-- Running the Linear solver --#

    CL_VLM[i], CDi_near_VLM, cl_VLM, cd_near_VLM, spanLocations, GammaValues_VLM = VLM(wingGeometry,freestream)
    CL_VLM_Waked[i], CDi_near_VLM, cl_VLM, cd_near_VLM, spanLocations, GammaValues_VLM = VLM(wingGeometry,wakedFreestream)

    #-- Running the Nonlinear solver --#

    CL[i], CDi, cl, spanLocations = solveBlownWing(wingGeometry,airfoil,airfoilName,freestream,rbfs)
    CL_Waked[i], CDi, cl, spanLocations = solveBlownWing(wingGeometry,airfoil,airfoilName,wakedFreestream,rbfs)

end

figure()
plot(angles, CL_VLM, label = "VLM No Props", color = "green", linestyle = "--")
plot(angles, CL_VLM_Waked, label = "VLM With Props", color = "green", linestyle = "-")
plot(angles, CL, label = "Strip Theory No Props", color = "orange", linestyle = "--")
plot(angles, CL_Waked, label = "Strip Theory With Props", color = "orange", linestyle = "-")
xlabel("Angle of Attack (Degrees)")
ylabel("CL")
title("Veldhuis Wing Simulation")
legend()
