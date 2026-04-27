import numpy as np

n = 10

a = np.random.randint(10,size=(n))
b = np.random.randint(10,size=(n))
c = np.random.randint(10,size=(n))
f = np.random.randint(10,size=(n))

print(a)
print(b)
print(c)
print(f)
print('')

A = np.zeros((n,n))
for i in range(n):
  A[i,i] = b[i]

  if i > 0:
    A[i,i-1] = a[i]
  else:
    A[i,-1] = a[i]

  if i < n-1:
    A[i,i+1] = c[i]
  else:
    A[i,0] = c[i]

print(A)
print('')

s = np.linalg.solve(A,f)
print(s)