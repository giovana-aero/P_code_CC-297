import numpy as np
import matplotlib.pyplot as plt

def main():

  address_list = ["../exercise_laplace2d/results"]
  casename = "laplace2d"

  execute(address_list,casename)

def execute(address_list,casename):

  plt.rcParams['text.usetex'] = True

  for address in address_list:

    if address[-1] != '/':
      address += '/'

    filename = address + casename + ".log"

    data = np.loadtxt(filename)

    max_iter = data.shape[0]

    fig,ax = plt.subplots(1,1)
    ax.semilogy(np.arange(1,max_iter+1),data)
  
  plt.grid(True)
  plt.xlabel("$n$")
  # plt.ylabel("$$")
  plt.title("finalizar isto!!!")
  plt.show()

if __name__ == "__main__":
  main()
