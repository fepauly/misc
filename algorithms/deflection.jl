#= Solve a numerical problem in Julia

deflection curve -> beam with area force

=#

using LinearAlgebra, Plots

l = 1. # length of beam
q0 = 200. # area force
E = 2-1e11 # Young’s modulus
J = 5.4e-7 # geometrical moment of inertia

n = 1000 # number of points
dx = l/(n-1)

# Initialize
x = [i*dx for i in 0:n-1]
d = [-2. for i in 0:n-1]
a = [1. for i in 0:n-2]
b = [1. for i in 0:n-2]

r_func(x) = (q0 * x * l * dx^2)/(6*E*J) * ((x/l)^2-1)
r = [r_func(i) for i in x]

# Border Conditions
a[1] = 0; 
b[n-1] = 0
d[1] = 1; d[n] = 1
r[1] = 0; r[n] = 0

# Create system of equations
tridiag = Tridiagonal(b, d, a)
# display(tridiag)

# Solve it
w = tridiag\r

# Simple plotting of the result
plot(x, w, label="Biegelinie")
