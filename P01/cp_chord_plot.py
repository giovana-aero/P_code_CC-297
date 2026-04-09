from numpy import loadtxt
import matplotlib.pyplot as plt

def main():

  casename = "bi_air"

  address_list = ["./mesh_conv/mtype1/",
                  "./mesh_conv/mtype2/",
                  "./mesh_conv/mtype3/",
                  "./mesh_conv/mtype4/"]
  ILE = [10,20,40,5]
  ITE = [30,60,120,15]

  execute(address_list,casename,ILE,ITE)

def execute(address_list,casename,ILE,ITE):
  base_size = 10

  plt.rcParams['text.usetex'] = True

  i = 0
  for address in address_list:
  
    if address[-1] != '/':
      address += '/'

    mesh_x = loadtxt(address + casename + '_mesh_x.msh')
    cp = loadtxt(address + casename + '_cp_chord.dat')

    chord = mesh_x[ILE[i]:ITE[i]+1]

    # plt.figure(figsize=(fig_width,fig_height))
    # plt.figure()
    plt.plot(chord,-cp,label=address)
    plt.grid(True)
    plt.xlabel("x")
    plt.ylabel("-cp")
    plt.title("cp chord")

    i += 1

  plt.legend()
  plt.show()

if __name__ == "__main__":
  main()