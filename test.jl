import Elliptic

include("EllipticIntegral.jl")

k = 0.5
m = k^2
n = 0.5

print(Elliptic.K(m))
print("\n")
print(elliptic_K(m))
print("\n")

print(Elliptic.E(m))
print("\n")
print(elliptic_E(m))
print("\n")

print(Elliptic.Π(n, π/2, m))
print("\n")
print(elliptic_Π(n, m))