// #include<stdlib.h>
// #include"../include/num_methods.h"
// #include"../include/eom_lib.h"

// // void evaluate_delta_form_fullp(sim_prmtrs *config_msh,msh_prmtrs *msh,
// //                             control_prmtrs *c_prmtrs,int init_only){
// //   // double (*x)[msh->IMAX] = calloc(msh->JMAX,sizeof *x);
// //   // double (*y)[msh->IMAX] = calloc(msh->JMAX,sizeof *y);

// //   config_msh.save_last_only = 1;

// //   if(init_only == 1){
// //     config_msh.save_i_c == 1
// //   }


// // }

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