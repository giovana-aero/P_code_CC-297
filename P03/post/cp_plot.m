clc,clear

address = '../results/';
casename = 'fullp';

save_fig = 0;
base_size = 800;
% fontsize = base_size*12/500;
fontsize = 22;
xlims = [];
ylims = [];

mesh.x = dlmread([address,casename,'_x_mesh.dat']);
mesh.y = dlmread([address,casename,'_y_mesh.dat']);

cp = dlmread([address,casename,'_cp.dat']);

chord_nump = (size(cp,2) + 1)/2;

figure(1),clf
% plot(mesh.x(1,1:chord_nump),mesh.y(1,1:chord_nump)),hold on,grid on
% plot(mesh.x(1,chord_nump:end),mesh.y(1,chord_nump:end))
plot(mesh.x(1,1:chord_nump),cp(1,1:chord_nump),'k'),hold on,grid on
plot(mesh.x(1,chord_nump:end),cp(1,chord_nump:end),'k')

inv_yaxis();

