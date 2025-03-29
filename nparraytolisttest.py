import numpy as np
A = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        A[i][j]=i+j
print(A)
A=A.tolist()
print(A)
