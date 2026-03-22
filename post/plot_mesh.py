from numpy import loadtxt,meshgrid
import matplotlib.pyplot as plt

def main():
  address = "../exercise_laplace2d/results/"
  casename = "laplace2d"

  execute(address,casename)

def execute(address,casename):
  base_size = 10

  plt.rcParams['text.usetex'] = True

  if address[-1] != '/':
    address += '/'

  mesh_x = loadtxt(address + casename + '_mesh_x.msh')
  mesh_y = loadtxt(address + casename + '_mesh_y.msh')

  X,Y = meshgrid(mesh_x,mesh_y)

  if mesh_x[-1] > mesh_y[-1]:
    fig_height = base_size*mesh_y[-1]/mesh_x[-1]
    fig_width = base_size
  else:
    fig_width = base_size*mesh_x[-1]/mesh_y[-1]
    fig_height = base_size

  plt.figure(figsize=(fig_width,fig_height))
  plt.plot(X,Y,marker='.',color='k',linestyle='none')
  plt.grid(True)
  plt.xlabel('x')
  plt.ylabel('y')
  plt.show()

if __name__ == "__main__":
  main()

