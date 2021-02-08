# using Pkg

# Importing the necessary julia packages
# using Revise
import CSV

using Revise
using Plots
Revise.revise()

# Importing the necessary FLOW packages
import AirfoilPrep
ap = AirfoilPrep

import AirfoilDatabase
adb = AirfoilDatabase

# macro for bash commands?
# Defining the database
database_path = "airfoil-data/NACA23018"
adb.new_database(database_path)

# # Defining the airfoil and flow conditions
# airfoil_file = "airfoil-data/NACA4415.csv"
# Re = 10e6
# Ma = 0.22
# alphas = [i for i in -30:1.0:30];

# # Reading in the airfoil
# x, y = ap.readcontour(airfoil_file; header_len = 1, delim = ",");

# # Running XFOIL to get the polar
# polar = ap.runXFOIL(x, y, Re; alphas = alphas, verbose = false, Mach = Ma, iter = 100);

# # Adding the polar to the database
# adb.new_entry(polar; database_path=database_path, airfoilname = "NACA 4415");

airfoil_file = "airfoil-data/NACA23018.csv"

alphas = [i for i in -30:1.0:30]
Res = 4e6:1e5:8e6
Ma = 0
ncrit = 9

# Read the airfoil contour
x, y = ap.readcontour(airfoil_file; header_len=1, delim=",")

#------RUN SWEEP-------------------------
for Re in Res
    
    println("Sweep at Re = $Re...")
    
    # Run XFOIL to create polar
    polar = ap.runXFOIL(x, y, Re; alphas = alphas, Mach = Ma, ncrit = ncrit, verbose = false, iter = 100,
                        clmaxstop = true, clminstop = true)

    # Viterna extrapolation
    println("Extrapolating...")
    polar = ap.extrapolate(polar,0.0;AR = 10) # CDmax = 0 as a dummy guess because CDmax is calculated from AR
    println("Done extrapolating!")
    
    # Add the newly-created polar to the database
    adb.new_entry(polar; database_path = database_path, airfoilname = "NACA 4415", warn = false)
    
end

# Defining the things we're looking for
pathToData = "airfoil-data/NACA23018"
airfoilName = "NACA23018"
resultType = "Cl"
ReynoldsNumber = 5000000
machNumber = 0
ncrit = 9

# Creating the file name
filename = string(pathToData,"/",resultType,"/",
                  airfoilName,"-",
                  resultType,"-",
                  "re",ReynoldsNumber,"-",
                  "ma",machNumber,"p0-",
                  "ncrit",ncrit,"p0-0.csv")

data = CSV.read(filename)

data = Matrix(data)

plotly()
Plots.plot(data[:,1],data[:,2],
            title = "Extrapolated Airfoil Data",
            legend = false,
            xlim = [-30,30])