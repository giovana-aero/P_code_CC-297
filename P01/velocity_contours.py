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
  u = loadtxt(address + casename + '_u.dat')
  v = loadtxt(address + casename + '_v.dat')
  Ve = loadtxt(address + casename + '_Ve.dat')
  
  casename += '_iter_'

  file_list = str(check_output(['ls',address]))
  all_phi = re.findall(r'(' + casename + r'\d{10}.dat)',file_list)

  # qtimes = int(all_phi[1][len(casename):-4])
  # last_save = int(all_phi[-1][len(casename):-4])
  # numfiles = len(all_phi)

  phi = []
  for save_file in all_phi:
    phi.append(loadtxt(address + save_file))

  if mesh_x[-1] > mesh_y[-1]:
    fig_height = base_size*mesh_y[-1]/mesh_x[-1]
    fig_width = base_size
  else:
    fig_width = base_size*mesh_x[-1]/mesh_y[-1]
    fig_height = base_size

  plt.figure(figsize=(fig_width,fig_height))
  plt.contourf(mesh_x,mesh_y,phi[-1],cmap="cool")
  plt.grid(True)
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("phi")

  plt.figure(figsize=(fig_width,fig_height))
  plt.contourf(mesh_x,mesh_y,u,cmap="cool")
  plt.grid(True)
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("u")

  plt.figure(figsize=(fig_width,fig_height))
  plt.contourf(mesh_x,mesh_y,v,cmap="cool")
  plt.grid(True)
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("v")

  plt.figure(figsize=(fig_width,fig_height))
  plt.contourf(mesh_x,mesh_y,Ve,cmap="cool")
  plt.grid(True)
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("Ve")

  plt.show()

if __name__ == "__main__":
  main()
