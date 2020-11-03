using Pkg

# Importing the necessary julia packages
using Revise
import CSV

# Making the environment here the same as on Eduardo's computer
Pkg.instantiate()

# Importing the necessary FLOW packages
import AirfoilPrep
ap = AirfoilPrep

import AirfoilDatabase
adb = AirfoilDatabase

# macro for bash commands?
# Defining the database
database_path = "airfoil-data/NACA4415"
adb.new_database(database_path)

# Defining the airfoil and flow conditions
airfoil_file = "airfoil-data/NACA4415.csv"
Re = 10e6
Ma = 0.22
alphas = [i for i in -30:1.0:30];

# Reading in the airfoil
x, y = ap.readcontour(airfoil_file; header_len = 1, delim = ",");

# Running XFOIL to get the polar
polar = ap.runXFOIL(x, y, Re; alphas = alphas, verbose = false, Mach = Ma, iter = 100);

# Adding the polar to the database
adb.new_entry(polar; database_path=database_path, airfoilname = "NACA 4415");

airfoil_file = "airfoil-data/NACA4415.csv"

alphas = [i for i in -30:1.0:30]
Res = 1e5:1e5:5e6
Ma = 0
ncrit = 9

# Read the airfoil contour
x, y = ap.readcontour(airfoil_file; header_len=1, delim=",")

#------RUN SWEEP-------------------------
for Re in Res
    
    println("Sweep at Re = $Re...")
    
    # Run XFOIL to create polar
    polar = ap.runXFOIL(x, y, Re; alphas = alphas, Mach = Ma, ncrit = ncrit, verbose = false, iter = 100)
    
    # Add the newly-created polar to the database
    adb.new_entry(polar; database_path = database_path, airfoilname = "NACA 4415", warn = false)
    
end

# Defining the things we're looking for
pathToData = "../airfoil-data/NACA4415"
airfoilName = "NACA4415"
resultType = "Cl"
ReynoldsNumber = 100000
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