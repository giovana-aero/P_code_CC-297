#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/num_methods.h"

void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,char *fname_msh_x,
                               char *fname_msh_y){
  double (*x)[n] = calloc(m,sizeof *x);
  double (*y)[n] = calloc(m,sizeof *y);

  read_2d_array_from_file(m,n,x,fname_msh_x);
  read_2d_array_from_file(m,n,y,fname_msh_y);

  


  free(x);
  free(y);
}

// void mesh_config(int c_type,control_prmtrs c_prmtrs,msh_prmtrs msh){
//   if(load_mesh){


//   }
//   else{
//     set_control_prmtrs(c_type,&c_prmtrs,&msh);
//     config_msh.casename = malloc(sizeof(char)*200);
//     sprintf(config_msh.casename,"%s",output_file_msh);
//     // switch(mesh_type){
//     //   case 1:
//     evaluate_delta_form_eom(&config_msh,&msh,&c_prmtrs,init_only);
//     //     break;
      
//     //   case 2:
//     //     evaluate_delta_form_ecm(&config,&msh,&c_prmtrs,init_only);
//     //     break;
      
//     //   default:
//     //     puts("invalid mesh_type");
//     //     return 1;
//     // }
//   }
// }

// void save_parameters(){

// }