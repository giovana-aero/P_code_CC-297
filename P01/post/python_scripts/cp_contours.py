# import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
from subprocess import check_output
import re

def main():

  address = "./results/"
  casename = "bi_air"

  execute(address,casename)

def execute(address,casename):
  base_size = 10

  plt.rcParams['text.usetex'] = True

  if address[-1] != '/':
    address += '/'

  mesh_x = loadtxt(address + casename + '_mesh_x.msh')
  mesh_y = loadtxt(address + casename + '_mesh_y.msh')
  cp = loadtxt(address + casename + '_cp.dat')
  
  if mesh_x[-1] > mesh_y[-1]:
    fig_height = base_size*mesh_y[-1]/mesh_x[-1]
    fig_width = base_size
  else:
    fig_width = base_size*mesh_x[-1]/mesh_y[-1]
    fig_height = base_size

  plt.figure(figsize=(fig_width,fig_height))
  plt.contourf(mesh_x,mesh_y,-cp,cmap="cool")
  plt.grid(True)
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("cp")

  plt.show()

if __name__ == "__main__":
  main()
