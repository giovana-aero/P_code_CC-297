clc,clear

% address = '../results/';
% casename = 'bi_air';

% address = '../final_results_slor/t005_mtype4/';
% address = '../final_results_slor/t010_mtype4/';
% address = '../final_results_slor/naca0005_mtype1/';
% address = '../final_results_slor/naca0005_mtype3/';
address = '../final_results_slor/naca0010_mtype1/';
% address = '../final_results_slor/naca0010_mtype3/';
casename = 'bi_air';

base_size = 1000;
fontsize = 14;
% pba = [2,1];
no_lines = 1;
save_figs = 1;
cmap = slanCM('bupu');
% cmap = slanCM('purples');

iter_list = dir([address,casename,'_iter_*']);

mesh_x = readmatrix([address,casename,'_mesh_x.msh'],'filetype','delimitedtext');
mesh_y = readmatrix([address,casename,'_mesh_y.msh'],'filetype','delimitedtext');
phi = readmatrix([address,iter_list(end).name],'filetype','delimitedtext');
u = readmatrix([address,casename,'_u.dat'],'filetype','delimitedtext');
v = readmatrix([address,casename,'_v.dat'],'filetype','delimitedtext');
Ve = readmatrix([address,casename,'_Ve.dat'],'filetype','delimitedtext');
cp = readmatrix([address,casename,'_cp.dat'],'filetype','delimitedtext');


lx = mesh_x(end) - mesh_x(1);
ly = mesh_y(end) - mesh_y(1);
if lx > ly
  fig_height = base_size*ly/lx;
  fig_width = base_size;
  pba = [fig_width,fig_height]/fig_width;
else
  fig_width = base_size*lx/ly;
  fig_height = base_size;
  pba = [fig_width,fig_height]/fig_height;
end

pos = [fig_width,fig_height];

figure(1),clf
[~,h] = contourf(mesh_x,mesh_y,phi);
xylabel_latex('x','y');
colormap(cmap)
% colorbar()
cb_latex(colorbar(),'\phi',fontsize);
set_fontsize_position(fontsize,pos,pba)
if no_lines,set(h,'LineColor','none'),end

figure(2),clf
[~,h] = contourf(mesh_x,mesh_y,u);
xylabel_latex('x','y');
colormap(cmap)
% colorbar()
cb_latex(colorbar(),'u',fontsize);
set_fontsize_position(fontsize,pos,pba)
if no_lines,set(h,'LineColor','none'),end

figure(3),clf
[~,h] = contourf(mesh_x,mesh_y,v);
xylabel_latex('x','y');
colormap(cmap)
cb_latex(colorbar(),'v',fontsize);
set_fontsize_position(fontsize,pos,pba)
if no_lines,set(h,'LineColor','none'),end

figure(4),clf
[~,h] = contourf(mesh_x,mesh_y,Ve);
xylabel_latex('x','y');
colormap(cmap)
cb_latex(colorbar(),'\vec{U}',fontsize);
set_fontsize_position(fontsize,pos,pba)
if no_lines,set(h,'LineColor','none'),end

figure(5),clf
[~,h] = contourf(mesh_x,mesh_y,cp);
xylabel_latex('x','y');
colormap(cmap)
cb_latex(colorbar(),'Cp',fontsize);
set_fontsize_position(fontsize,pos,pba)
if no_lines,set(h,'LineColor','none'),end

if save_figs
  filenames = {'contours_phi','contours_u','contours_v','contours_Ve',...
               'contours_cp'};
  for i = 1:5
    figure(i)
    set_rndrr_painters();
    save_fig_eps(filenames{i});
  end
end