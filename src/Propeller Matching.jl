using Pkg
Pkg.add("Revise")
Pkg.add("CSV")
Pkg.add("PyPlot")
Pkg.add("CCBlade")
using Revise
using CSV
using PyPlot
using CCBlade
revise()

#-- Defining the propeller --#
airfoilData = "/Users/markanderson/Box/FLOW-MCA/Code/blownWing/Propeller Data/E212 Propeller Data.dat"
Rhub = 0.01
Rtip = 0.236 / 2
numBlades = 4 # Four blades

# Blade-wise properties
r = [0.148, 0.254237, 0.381356, 0.508475, 0.635593, 0.762712, 0.889831, 0.99]*Rtip;
chord = [9.88, 11.88, 15.59, 18.81, 19.55, 18.32, 13.96, 0.01].*1e-3*Rtip;
theta = [35.0, 32.5, 26.5, 23.5, 19, 16.5, 14.0, 10.0].*pi/180;

# Defining operating point
J = 0.85
rho = 1.225
D = 2*Rtip
Vinf = 50
n = Vinf / (J*D)
Omega = 2*pi*n

#-- Standard CCBlade Code --#

af = af_from_files(airfoilData)
airfoils = fill(af,length(r))

rotor = Rotor(Rhub, Rtip, numBlades, false)
sections = Section.(r, chord, theta, airfoils);

op = simple_op.(Vinf, Omega, r, rho);

# Getting and plotting the outputs
outputs = solve.(Ref(rotor), sections ,op)

# Plotting the wake velocity
figure(10)
plot(r/Rtip, outputs.u/Vinf)
plot(r/Rtip, outputs.v/Vinf)
title("Wake Velocity")
xlabel("r/Rtip")
ylabel("(normalized) induced velocity at rotor disk")
legend(["axial velocity", "swirl velocity"]);

#-- Calculating the Thrust Coefficient --#

T, Q = thrusttorque(rotor, sections, outputs)

eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor)
println("Efficiency = ", eff)
println("Thrust Coefficient = ", CT)
