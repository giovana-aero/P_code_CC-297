clc,clear

address = '../results/';
casename = 'bi_air';
cmap = slanCM('bupu');
save_figs = 1;

t = 0.05;
ILE = 11;
ITE = 31;

fontsize = 14;
scttrsz = 10;
lnwdth = 1.5;
base_size = 1000;


mesh_x = readmatrix([address,casename,'_mesh_x.msh'],'filetype','delimitedtext');
mesh_y = readmatrix([address,casename,'_mesh_y.msh'],'filetype','delimitedtext');

[X,Y] = meshgrid(mesh_x,mesh_y);

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
scatter(X,Y,scttrsz,'k','filled'),grid on
xlim([-2.5,3.5])
ylim([-0.2,2.3])
xylabel_latex('x','y');
set_fontsize_position(fontsize,pos,pba);



x = mesh_x(ILE:ITE);

figure(2),clf
scatter(X,Y,scttrsz,'k','filled'),grid on,hold on
plot(x,bi_air_shape(x,t),'color',cmap(150,:),'linewidth',lnwdth)
plot(x,-bi_air_shape(x,t),'color',cmap(150,:),'linewidth',lnwdth)
xylabel_latex('x','y');
xlim([-.2,1.2])
ylim([-.1,.2])
set_fontsize_position(fontsize,pos,pba);


if save_figs
  figure(1)
  set_rndrr_painters();
  save_fig_eps('mesh1');

  figure(2)
  set_rndrr_painters();
  save_fig_eps('mesh1_bi_air')
end

%% 

function y = bi_air_shape(x,t)
  y = 2*t*x.*(1 - x);
end