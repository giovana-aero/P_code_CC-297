import numpy as np

# note: the nomenclature is based on the upstream case

def main():
  # # x, original
  # x1 = -2.078306
  # x_le = 0
  # x_lep1 = 0.05
  # n = 11

  # x, downstream of the airfoil
  x1 = 3.078306
  x_le = 1
  x_lep1 = 0.95
  n = 11

  # # x, 1/2
  # x1 = -2.078306
  # x_le = 0
  # x_lep1 = 0.025
  # n = 21

  # # x, 1/4
  # x1 = -2.078306
  # x_le = 0
  # x_lep1 = 0.025/2
  # n = 41

  # # x, 1/8
  # x1 = -2.078306
  # x_le = 0
  # x_lep1 = 0.05/8
  # n = 81

  # # x, 1/.5
  # x1 = -2.078306
  # x_le = 0
  # x_lep1 = 0.1
  # n = 6

  # # y, original
  # x1 = 2.103306
  # x_le = 0.025000
  # x_lep1 = -0.025000
  # n = 11

  # # y, 1/2
  # x1 = 2.103306
  # x_le = 0.025000/2
  # x_lep1 = -0.025000/2
  # n = 23

  # # y, 1/4
  # x1 = 2.103306
  # x_le = 0.025000/4
  # x_lep1 = -0.025000/4
  # n = 47

  # # y, 1/8
  # x1 = 2.103306
  # x_le = 0.025000/8
  # x_lep1 = -0.025000/8
  # n = 95

  # # y, 1/.5
  # x1 = 2.103306
  # x_le = 0.025000*2
  # x_lep1 = -0.025000*2
  # n = 5

  execute(x1,x_le,x_lep1,n)

def execute(x1,x_le,x_lep1,n):
  p = np.ones(n)

  p[-1] = (x_le - x1)/(x_le - x_lep1)

  sf = np.roots(p)

  print(sf)

if __name__ == "__main__":
  main()