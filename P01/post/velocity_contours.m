clc,clear

% address = '../results/';
% casename = 'bi_air';

address = '../final_results_slor/t005/';
% address = '../final_results_slor/t010/';
casename = 'bi_air';

base_size = 1000;
fontsize = 14;
pba = [2,1];
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


if mesh_x(end) > mesh_y(end)
  fig_height = base_size*mesh_y(end)/mesh_x(end);
  fig_width = base_size;
else
  fig_width = base_size*mesh_x(end)/mesh_y(end);
  fig_height = base_size;
end

figure(1),clf
contourf(mesh_x,mesh_y,phi)
xylabel_latex('x','y');
colormap(cmap)
% colorbar()
cb_latex(colorbar(),'\phi',fontsize);
set_fontsize_position(fontsize,[fig_width,fig_height],pba)

figure(2),clf
contourf(mesh_x,mesh_y,u)
xylabel_latex('x','y');
colormap(cmap)
% colorbar()
cb_latex(colorbar(),'u',fontsize);
set_fontsize_position(fontsize,[fig_width,fig_height],pba)

figure(3),clf
contourf(mesh_x,mesh_y,v)
xylabel_latex('x','y');
colormap(cmap)
cb_latex(colorbar(),'v',fontsize);
set_fontsize_position(fontsize,[fig_width,fig_height],pba)

figure(4),clf
contourf(mesh_x,mesh_y,Ve)
xylabel_latex('x','y');
colormap(cmap)
cb_latex(colorbar(),'\vec{U}',fontsize);
set_fontsize_position(fontsize,[fig_width,fig_height],pba)

figure(5),clf
contourf(mesh_x,mesh_y,cp)
xylabel_latex('x','y');
colormap(cmap)
cb_latex(colorbar(),'Cp',fontsize);
set_fontsize_position(fontsize,[fig_width,fig_height],pba)

if save_figs
  filenames = {'contours_phi','contours_u','contours_v','contours_Ve',...
               'contours_cp'};
  for i = 1:5
    figure(i)
    set_rndrr_painters();
    save_fig_eps(filenames{i});
  end
end