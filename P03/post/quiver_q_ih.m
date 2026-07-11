clc,clear

address = '../results/';
casename = 'fullp';

mesh.x = dlmread([address,casename,'_x_mesh_ih.dat']);
mesh.y = dlmread([address,casename,'_y_mesh_ih.dat']);
u = dlmread([address,casename,'_u_ih.dat']);
v = dlmread([address,casename,'_v_ih.dat']);

x = mesh.x;
y = mesh.y;

figure(1),clf,hold on
quiver(mesh.x(:,1:end-1),mesh.y(:,1:end-1),u(:,1:end-1),v(:,1:end-1),.1);
grid on,axis equal
