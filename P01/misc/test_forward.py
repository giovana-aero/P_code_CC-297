x = [-2.078306,-1.612645,-1.240116]
phi = [-2.078306,-1.612645,-1.240116]

# forward
i = 0
# a = 1
# b = -3/2
# c = -1/2
# d = 1
# dphi = (a*phi[i+1] + b*phi[i])/(x[i+1] - x[i]) + (c*phi[i+2] + phi[i+1])/(x[i+2] - x[i+1])

# dphi = ((-3*phi[i] + 3*phi[i+1]) + (phi[i+1] - phi[i+2]))/((x[i+1] - x[i]) + (x[i+2] - x[i+1]))

dphi = (-3*phi[i] + 4*phi[i+1] - phi[i+2])/(-3*x[i] + 4*x[i+1] - x[i+2])

print(dphi)