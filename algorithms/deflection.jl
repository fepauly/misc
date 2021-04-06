#= Solve a numerical problem in Julia

Biegelinie -> Balken mit Streckenlast

=#

using LinearAlgebra, Plots

l = 1. # Länge des Balken
q0 = 200. # Streckenlast
E = 2-1e11 # E-Modul
J = 5.4e-7 # Flächenträgheitsmoment

n = 1000 # Anzahl Aufpunkte
dx = l/(n-1)

# Initialisierung der GLS
x = [i*dx for i in 0:n-1]
d = [-2. for i in 0:n-1]
a = [1. for i in 0:n-2]
b = [1. for i in 0:n-2]

r_func(x) = (q0 * x * l * dx^2)/(6*E*J) * ((x/l)^2-1)
r = [r_func(i) for i in x]

# Randbedinungen
a[1] = 0; 
b[n-1] = 0
d[1] = 1; d[n] = 1
r[1] = 0; r[n] = 0

tridiag = Tridiagonal(b, d, a)
# display(tridiag)

w = tridiag\r

plot(x, w, label="Biegelinie")
