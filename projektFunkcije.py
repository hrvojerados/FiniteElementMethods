import numpy as np
#Nx i Ny je broj čvorova po horiznotalnoj liniji, a Ny po vertikalnoj
#grid je brElem x 3 matrica
import numpy as np

def generateRectangleMesh(grid, Nx, Ny, h):
    nodeIndices = np.arange(Nx * Ny).reshape(Ny, Nx)
    
    edgeNodes = []
    edgeNodes.extend(nodeIndices[0, :])        
    edgeNodes.extend(nodeIndices[-1, :])      
    edgeNodes.extend(nodeIndices[1:-1, 0])    
    edgeNodes.extend(nodeIndices[1:-1, -1])   
    edgeNodes = list(set(edgeNodes))
    #print(edgeNodes)

    innerNodes = list(set(np.arange(Nx * Ny)) - set(edgeNodes))
    newOrder = innerNodes + edgeNodes
    indexMap = {old: new for new, old in enumerate(newOrder)}
    reverseIndexMap = {new: old for old, new in indexMap.items()}  
    #print(indexMap)
    elemNum = 0
    for yIndex in range(0, Ny - 1):
        pos = ((yIndex + 1) // 2) * 2 * Nx
        goingUp = (yIndex % 2 == 0)
        for xIndex in range(0, Nx - 1):
            nextPos = pos + (1 - Nx + 2 * Nx * goingUp)
            grid[elemNum][0] = indexMap[pos]
            grid[elemNum][1] = indexMap[nextPos]
            grid[elemNum][2] = indexMap[pos + 1]
            elemNum += 1
            grid[elemNum][0] = indexMap[nextPos]
            grid[elemNum][1] = indexMap[pos]
            grid[elemNum][2] = indexMap[nextPos - 1]
            elemNum += 1
            pos = nextPos
            goingUp = not goingUp
    
    return reverseIndexMap

def nodeToCoords(node, reverseIndexMap, Nx, h):
    originalNode = reverseIndexMap[node]
    return ((originalNode % Nx) * h, (originalNode // Nx) * h)

#A je incijalizirana kao nul matrica
def createStifnesMat(A, grid, Nx, Ny, h, reverseIndexMap):
    elNum = (Nx - 1) * (Ny - 1) * 2
    gradPhi0 = np.array([-1, -1])
    gradPhi1 = np.array([1, 0])
    gradPhi2 = np.array([0, 1])
    limit = (Nx - 2) * (Ny - 2)
    for i in range(0, elNum):
        node0 = int(grid[i][0])
        node1 = int(grid[i][1])
        node2 = int(grid[i][2])
        xy0 = nodeToCoords(node0, reverseIndexMap, Nx, h)
        print(node0)
        x0 = xy0[0]
        print(x0)
        y0 = xy0[1]
        print(y0)

        xy1 = nodeToCoords(node1, reverseIndexMap, Nx, h)
        print(node1)
        x1 = xy1[0]
        print(x1)
        y1 = xy1[1]
        print(y1)

        xy2 = nodeToCoords(node2, reverseIndexMap, Nx, h)
        print(node2)
        x2 = xy2[0]
        print(x2)
        y2 = xy2[1]
        print(y2)
        det = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        Jc = np.array([[y2 - y0, y0 - y1], [x0 - x2, x1 - x0]])
        print(Jc)
        #0 i 1
        isEdge0 = node0 >= limit 
        isEdge1 = node1 >= limit
        isEdge2 = node2 >= limit 
        if ((not isEdge0) and (not isEdge1)):
            A[node0][node1] += np.dot(np.matmul(gradPhi0, Jc), np.matmul(gradPhi1, Jc)) / (2 * det) 
            A[node1][node0] += np.dot(np.matmul(gradPhi0, Jc), np.matmul(gradPhi1, Jc)) / (2 * det)

        #1 i 2
        if ((not isEdge1) and (not isEdge2)):
            A[node1][node2] += np.dot(np.matmul(gradPhi1, Jc), np.matmul(gradPhi2, Jc)) / (2 * det) 
            A[node2][node1] += np.dot(np.matmul(gradPhi1, Jc), np.matmul(gradPhi2, Jc)) / (2 * det)
        #0 i 2
        if ((not isEdge0) and (not isEdge2)):
            A[node0][node2] += np.dot(np.matmul(gradPhi0, Jc), np.matmul(gradPhi2, Jc)) / (2 * det) 
            A[node2][node0] += np.dot(np.matmul(gradPhi0, Jc), np.matmul(gradPhi2, Jc)) / (2 * det)
        #0 i 0
        if (not isEdge0):
            print(Jc)
            print(gradPhi0)
            print(np.matmul(gradPhi0, Jc))
            print(np.dot(np.matmul(gradPhi0, Jc), np.matmul(gradPhi0, Jc)) / (2 * det))
            A[node0][node0] += np.dot(np.matmul(gradPhi0, Jc), np.matmul(gradPhi0, Jc)) / (2 * det) 
            print(A[node0][node0])
        #1 i 1
        if (not isEdge1):
            print(Jc)
            print(gradPhi1)
            print(np.matmul(gradPhi1, Jc))
            print(np.dot(np.matmul(gradPhi1, Jc), np.matmul(gradPhi1, Jc)) / (2 * det))
            A[node1][node1] += np.dot(np.matmul(gradPhi1, Jc), np.matmul(gradPhi1, Jc)) / (2 * det) 
            print(A[node1][node1])
        #2 i 2
        if (not isEdge2):
            print(Jc)
            print(gradPhi2)
            print(np.matmul(gradPhi2, Jc))
            print(np.dot(np.matmul(gradPhi2, Jc), np.matmul(gradPhi2, Jc)) / (2* det))
            A[node2][node2] += (np.dot(np.matmul(gradPhi2, Jc), np.matmul(gradPhi2, Jc))) / (2 * det) 
            print(A[node2][node2])
#b je inicijaliziran nulama
def createFreeTerm(b, grid, Nx, Ny, h, f, reverseIndexMap):
    elNum = (Nx - 1) * (Ny - 1) * 2
    limit = (Nx - 2) * (Ny - 2)
    for i in range(0, elNum):
        node0 = int(grid[i][0])
        node1 = int(grid[i][1])
        node2 = int(grid[i][2])
        xy0 = nodeToCoords(node0, reverseIndexMap, Nx, h)
        x0 = xy0[0]
        y0 = xy0[1]

        xy1 = nodeToCoords(node1, reverseIndexMap, Nx, h)
        x1 = xy1[0]
        y1 = xy1[1]

        xy2 = nodeToCoords(node2, reverseIndexMap, Nx, h)
        x2 = xy2[0]
        y2 = xy2[1]


        det = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        #lambda za transformaciju 
        transformX = lambda x, y: x0 + (x1 - x0) * x + (x2 - x0) * y
        transformY = lambda x, y: y0 + (y1 - y0) * x + (y2 - y0) * y
        gaussNodes = np.array([(0, 0, 1 / 3), (1, 0, 1 / 3), (0, 1, 1 / 3)])

        #0
        if (node0 >= limit):
            continue
        for xw, yw, w in gaussNodes: #provjeri ovo
            xHat = transformX(xw, yw)
            yHat = transformY(xw, yw)
            b[node0] += w * f(xHat, yHat) * (1 - xHat - yHat) * det
        #1
        if (node1 >= limit):
            continue
        for xw, yw, w in gaussNodes: #provjeri ovo
            xHat = transformX(xw, yw)
            yHat = transformY(xw, yw)
            b[node1] += w * f(xHat, yHat) * xHat * det
        #2
        if (node2 >= limit):
            continue
        for xw, yw, w in gaussNodes: #provjeri ovo
            xHat = transformX(xw, yw)
            yHat = transformY(xw, yw)
            b[node2] += w * f(xHat, yHat) * yHat * det

def f(x, y):
    #return x + y
    #return 0
    return 2 * (np.pi**2) * np.sin(np.pi * x) * np.sin(np.pi * y)

    

x0 = float(input())
y0 = float(input()) 
Nx = int(input())
h = x0 / (Nx - 1)
Ny = int((y0 / h) + 1)
nodeNum = Nx * Ny
elementNumber = (Nx - 1) * (Ny - 1) * 2
grid = np.zeros((elementNumber, 3))
reverseIndexMap = generateRectangleMesh(grid, Nx, Ny, h)
baseFooNum = (Nx - 2) * (Ny - 2)
#print(grid)
A = np.zeros((baseFooNum, baseFooNum))
createStifnesMat(A, grid, Nx, Ny, h, reverseIndexMap)
print(A)
b = np.zeros(baseFooNum)
b = np.transpose(b)

createFreeTerm(b, grid, Nx, Ny, h, f, reverseIndexMap)
print(b)
from numpy.linalg import solve
u = solve(A, b)
print(u)
full_u = np.zeros(Nx * Ny)
for idx, value in enumerate(u):
    full_u[idx] = value
for i in range (nodeNum - 2 * Nx - 2 * Ny + 4,nodeNum):
    full_u[i] = 0


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri

coordsX = []
coordsY = []
for i in range(0, nodeNum):
    x, y = nodeToCoords(i, reverseIndexMap, Nx, h)
    coordsX.append(x)
    coordsY.append(y)
x_coords = np.array(coordsX)
y_coords = np.array(coordsY)

triangles = grid.astype(int)
triangulation = tri.Triangulation(x_coords, y_coords, triangles)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_trisurf(triangulation, full_u, cmap='viridis', edgecolor='none')

fig.colorbar(surf, ax=ax, label='Solution')

ax.set_title('FEM Solution - 3D Plot')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Solution Value')
plt.show()
