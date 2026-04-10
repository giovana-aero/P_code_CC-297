import numpy as np
import matplotlib.pyplot as plt

def main():

  # address_list = ["../exercise_laplace2d/results"]
  # casename = "laplace2d"

  # address_list = ["../P01/results",]
  # address_list = ["../P01/results","../P01/results"]
  
  # # r optimization, sor
  # address_list = ["../P01/r_optimization/sor/gs/",
  #                 "../P01/r_optimization/sor/r08/",
  #                 "../P01/r_optimization/sor/r10/",
  #                 "../P01/r_optimization/sor/r12/",
  #                 "../P01/r_optimization/sor/r14/",
  #                 "../P01/r_optimization/sor/r16/",
  #                 "../P01/r_optimization/sor/r18/",
  #                 "../P01/r_optimization/sor/r182/",
  #                 "../P01/r_optimization/sor/r184/",
  #                 "../P01/r_optimization/sor/r186/",
  #                 "../P01/r_optimization/sor/r188/",
  #                 "../P01/r_optimization/sor/r20/",
  #                 "../P01/r_optimization/sor/r22/"]
  # legends = ["$GS$","$r = 0.8$","$r = 1.0$","$r = 1.2$","$r = 1.4$","$r = 1.6$",
  #            "$r = 1.8$","$r = 1.82$","$r = 1.84$","$r = 1.86$","$r = 1.88$",
  #            "$r = 2.0$","$r = 2.2$"]
  # linestyles = ['-','-','--','-','-','-','-','-','-','-','-','-','-'] 

  # r optimization, slor
  address_list = ["../P01/r_optimization/slor/lgs/",
                  "../P01/r_optimization/slor/r08/",
                  "../P01/r_optimization/slor/r10/",
                  "../P01/r_optimization/slor/r12/",
                  "../P01/r_optimization/slor/r14/",
                  "../P01/r_optimization/slor/r16/",
                  "../P01/r_optimization/slor/r18/",
                  "../P01/r_optimization/slor/r182/",
                  "../P01/r_optimization/slor/r184/",
                  "../P01/r_optimization/slor/r186/",
                  "../P01/r_optimization/slor/r188/",
                  "../P01/r_optimization/slor/r20/",
                  "../P01/r_optimization/slor/r22/"]
  legends =["$LGS$","$r = 0.8$","$r = 1.0$","$r = 1.2$","$r = 1.4$","$r = 1.6$",
            "$r = 1.8$","$r = 1.82$","$r = 1.84$","$r = 1.86$","$r = 1.88$",
            "$r = 2.0$","$r = 2.2$"]
  linestyles = ['-','-','--','-','-','-','-','-','-','-','-','-','-'] 




  casename = "bi_air"

  conv_ref = 1e-6

  execute(address_list,casename,legends,linestyles,conv_ref)

def execute(address_list,casename,legends,linestyles,conv_ref):

  # plt.rcParams['text.usetex'] = True

  fig,ax = plt.subplots(1,1)

  i = 0
  for address in address_list:

    if address[-1] != '/':
      address += '/'

    filename = address + casename + ".log"

    data = np.loadtxt(filename)

    max_iter = data.shape[0]

    ax.semilogy(np.arange(1,max_iter+1),data,linestyles[i],label=legends[i])#,linewidth=2)

    i += 1
  
  plt.grid(True)
  plt.xlabel("$n$")
  # plt.ylabel("$$")
  plt.title("finalizar isto!!!")
  plt.legend()
  # plt.plot([0,max_iter],[conv_ref,conv_ref],'m')
  plt.show()

if __name__ == "__main__":
  main()
