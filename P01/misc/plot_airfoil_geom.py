import numpy as np
import matplotlib.pyplot as plt

def main():
  t = 0.05
  ILE = 11
  ITE = 31
  mesh_x_file = "../results/bi_air_mesh_x.msh"
  mesh_y_file = "../results/bi_air_mesh_y.msh"

  execute(t,ILE,ITE,mesh_x_file,mesh_y_file)
  
def execute(t,ILE,ITE,mesh_x_file,mesh_y_file):
  base_size = 10
  xlims = [-0.25,1.25]
  ylims = [-0.1,0]

  plt.rcParams['text.usetex'] = True

  mesh_x = np.loadtxt(mesh_x_file)
  mesh_y = np.loadtxt(mesh_y_file)

  X,Y = np.meshgrid(mesh_x,mesh_y)

  if mesh_x[-1] > mesh_y[-1]:
    fig_height = base_size*mesh_y[-1]/mesh_x[-1]
    fig_width = base_size
    ylims[1] = np.diff(xlims)*mesh_y[-1]/mesh_x[-1]
  else:
    fig_width = base_size*mesh_x[-1]/mesh_y[-1]
    fig_height = base_size
    ylims[1] = np.diff(xlims)*mesh_x[-1]/mesh_y[-1]

  y_af = bi_air_shape(mesh_x[ILE-1:ITE],t)

  plt.figure(figsize=(fig_width,fig_height))
  plt.plot(X,Y,marker='.',color='k',linestyle='none')
  plt.grid(True)
  plt.xlabel('$x$')
  plt.ylabel('$y$')

  plt.plot(mesh_x[ILE-1:ITE],y_af,'m')
  plt.plot(mesh_x[ILE-1:ITE],-y_af,'m')
  plt.plot(mesh_x[ILE-1:ITE],np.zeros(len(mesh_x[ILE-1:ITE])),'--k')

  plt.xlim(xlims)
  plt.ylim(ylims)
  
  plt.show()

def bi_air_shape(x,t):
  return 2*t*x*(1 - x)

if __name__ == "__main__":
  main()