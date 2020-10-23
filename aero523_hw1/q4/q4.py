import matplotlib.pyplot as plt
import numpy as np

def read_txt(letter):
    f = open('%s.txt' % letter, 'r')
    an = f.readlines()

    B = []
    for each in an:
        temp = []
        splitnums = each.strip().split(' ')
        for num in splitnums:
            temp.append(float(num))
        B.append(temp)
    f.close()
    return B

V = np.array(read_txt('V'))
E = np.array(read_txt('E')).astype(int) - 1


old_loop = [E[0, 0]]
node1 = E[0, 0]
node2 = 0

k = 0
in_loop = True
while in_loop:
    for j in range(len(E)):
        node2 = E[j, 1]
        if old_loop[0] == node2 and node1 != E[j, 1]:
            in_loop = False
        if node1 == node2:
            node1 = E[j, 0]
            old_loop.append(node2)
            old_loop.append(node1)
            break



f = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
for i in range(len(old_loop)):
    node1 = V[E[old_loop[i], 0], :]
    node2 = V[E[old_loop[i], 1], :]

    plt.plot(np.array([node1[0], node2[0]]), np.array([node1[1], node2[1]]))

plt.xlabel(r'X-Axis', fontsize=16)
plt.ylabel(r'Y-Axis', fontsize=16)
plt.show()
