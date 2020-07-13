using Pkg
Pkg.add("Revise")
using Revise
revise()

Pkg.develop(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/Code/blownWing"))


#-- Defining the wing geometry --#

# the coordinates are defined starting with the leading edge root, and then clockwise from there,
# units are in meters

coordinates = [0     0     0;
               0     0.640 0;
               0.240 0.640 0;
               0.240 0     0];

numPanelsSpan = 50;

blownWing.generateWingGeometry(coordinates, numPanelsSpan)