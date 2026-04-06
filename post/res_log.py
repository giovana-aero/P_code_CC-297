import numpy as np
import matplotlib.pyplot as plt

def main():

  # address_list = ["../exercise_laplace2d/results"]
  # casename = "laplace2d"

  # address_list = ["../P01/results","../P01/results"]
  address_list = ["../P01/results"]
  casename = "bi_air"

  execute(address_list,casename)

def execute(address_list,casename):

  plt.rcParams['text.usetex'] = True

  fig,ax = plt.subplots(1,1)

  for address in address_list:

    if address[-1] != '/':
      address += '/'

    filename = address + casename + ".log"

    data = np.loadtxt(filename)

    max_iter = data.shape[0]

    ax.semilogy(np.arange(1,max_iter+1),data,label=address)
  
  plt.grid(True)
  plt.xlabel("$n$")
  # plt.ylabel("$$")
  plt.title("finalizar isto!!!")
  plt.legend()
  plt.show()

if __name__ == "__main__":
  main()
