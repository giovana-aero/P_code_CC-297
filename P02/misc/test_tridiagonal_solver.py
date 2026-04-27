import numpy as np

n = 15

a = np.random.randint(1,10,size=(n-1))
b = np.random.randint(1,10,size=(n))
c = np.random.randint(1,10,size=(n-1))
f = np.random.randint(1,10,size=(n))

print(a)
print(b)
print(c)
print(f)
print('')

A = np.zeros((n,n))
for i in range(n):
  A[i,i] = b[i]

  if i > 0:
    A[i,i-1] = a[i-1]

  if i < n-1:
    A[i,i+1] = c[i]

print(A)
print('')

s = np.linalg.solve(A,f)
print(s)