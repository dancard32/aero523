import numpy as np
import math

def write_file(v, filename):
    f = open(filename,"w")

    output = ''
    for i in range(len(v)):
        output += ' ' + str.format('{0:.5f}', v[i]) + ' ' + r'\\' # Output results to LaTeX environment
    f.write(output)
    f.close()

def write_w(filename, projection):
    f = open(filename, "w")
    
    output = ''
    for i in range(len(projection)):
        output += str.format('{0:.5f}', projection[i]) + r'\\ ' # Output results to LaTeX environment
    f.write(output)
    f.close()

# Pre-allocate matrices
u = np.array([[1, 0, 3], [2, -1, 0], [0, 0, 1], [0, 3, 0], [0, 0, -2]])
v = np.zeros(u.shape)

# Apply Gram-Schmidt Algorithm
for i in range(min(u.shape)):
    v[:,i] = u[:,i]
    for j in range(i):
        if j >= 0:
            v[:,i] = v[:,i] - np.dot(v[:,j], v[:,i])*v[:,j]     # Apply inner product
    v[:,i] = v[:,i]/math.sqrt(np.dot(v[:,i], v[:,i]))    # Normalize

# Output results to LaTeX
write_file(v[:,0], "v1")
write_file(v[:,1], "v2")
write_file(v[:,2], "v3")

w = np.linspace(1, 5, len(v), dtype = int, endpoint=True)
w_parallel = np.zeros(5)
for i in range(min(u.shape)):
    print(i)
    v_norm = v[:,i] / np.linalg.norm(v[:,i])
    w_parallel += np.dot(w,v_norm)*v_norm

w_perp = w - w_parallel
write_w("w_parallel", w_parallel)
write_w("w_perp", w_perp)
