import numpy as np
import math

# Pre-allocating values
N = 10  # Number of iterations 
u0 = np.matrix([1.5, 0.5])  # Initial guess
u = np.zeros([N, 2])
resid = np.zeros([N,1])
u[0,:] = u0

def f(dat):
    f1 = float(dat[0]**2 + math.sin(dat[1]) - 3)   # f_1 function
    f2 = float(dat[1]**3 - 2*dat[0]/dat[1] + 10)   # f_2 function

    return np.matrix([[f1], [f2]])

def print_results(Uk, R):
    f = open('q3_results',"w")  # Filename

    output = ''
    for i in range(len(u)-1):
        output += str.format('{0:.1f}',i) + r'& ('+str.format('{0:.8f}',Uk[i,0])+', '+str.format('{0:.8f}',Uk[i,1])+ r') &'+str.format('{0:.5e}',R[i,0])+ r'\\' # Output results to LaTeX environment
    f.write(output)
    f.close()

def final_vals(u):
    f = open('q3_final_vals',"w")  # Filename

    idx = len(u) - 1
    f.write('(' + str.format('{0:.5f}',u[idx,0])+', '+str.format('{0:.5f}',u[idx,1]) + ')')
    f.close()

resid[0] = float(np.linalg.norm(np.transpose(f(u[0,:]))))
for i in range(N-1):
    jacobian = np.matrix([[2*u[i,0] , math.cos(u[i,1])], [-2/u[i,1], 3*u[i,1]**2 + 2*u[i,0]/(u[i,1]**2)]]) # Compute partial R/partial U @(U_k)
    deltaux  = -np.linalg.inv(jacobian) * f(u[i,:]) # Compute the delta_Ux

    u[i+1,:] = u[i,:] + np.array([deltaux[0,0], deltaux[1,0]])  # Compute the next step
    resid[i+1] = np.linalg.norm(np.transpose(f(u[i+1,:])))    # Compute residual

print(resid)
print_results(u, resid) # Output results to LaTeX table
final_vals(u)