# import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
from subprocess import check_output
import re

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
  
  casename += '_iter_'

  file_list = str(check_output(['ls',address]))
  all_results = re.findall(r'(' + casename + r'\d{10}.dat)',file_list)


  # qtimes = int(all_results[1][len(casename):-4])
  # last_save = int(all_results[-1][len(casename):-4])
  # numfiles = len(all_results)

  results = []
  for save_file in all_results:
    results.append(loadtxt(address + save_file))

  if mesh_x[-1] > mesh_y[-1]:
    fig_height = base_size*mesh_y[-1]/mesh_x[-1]
    fig_width = base_size
  else:
    fig_width = base_size*mesh_x[-1]/mesh_y[-1]
    fig_height = base_size

  plt.figure(figsize=(fig_width,fig_height))
  plt.contourf(mesh_x,mesh_y,results[-1],cmap="cool")
  plt.grid(True)
  plt.colorbar()
  plt.xlabel("x")
  plt.ylabel("y")
  # plt.colormap("cool")
  plt.show()
  
  # op = 1
  # while op != -2:
  #   check_output(["clear"])
  #   print("qtimes = " + str(qtimes))
  #   print("Options:")
  #   print("0 (0000000000) to " + str(numfiles) + "(" + str(last_save) +");")
  #   print("-1 to animate evolution;")
  #   print("-2 to quit")
  #   op = int(input("Selection: "))

  #   if op > -1:
  #     plt.contourf()



    
  



  


if __name__ == "__main__":
  main()
