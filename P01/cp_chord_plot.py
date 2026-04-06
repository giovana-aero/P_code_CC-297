from numpy import loadtxt
import matplotlib.pyplot as plt

def main():

  address = "./results/"
  casename = "bi_air"

  ILE = 10
  ITE = 30

  execute(address,casename,ILE,ITE)

def execute(address,casename,ILE,ITE):
  base_size = 10

  plt.rcParams['text.usetex'] = True

  if address[-1] != '/':
    address += '/'

  mesh_x = loadtxt(address + casename + '_mesh_x.msh')
  cp = loadtxt(address + casename + '_cp_chord.dat')

  chord = mesh_x[ILE:ITE+1]

  # plt.figure(figsize=(fig_width,fig_height))
  plt.figure()
  plt.plot(chord,-cp)
  plt.grid(True)
  plt.xlabel("x")
  plt.ylabel("-cp")
  plt.title("cp chord")

  plt.show()

if __name__ == "__main__":
  main()