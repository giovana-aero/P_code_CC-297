clc,clear

address = './';
casename = 'test';

% mesh_x = readmatrix([address,casename,'_mesh_x.msh'],'filetype','delimitedtext');
mesh_x = dlmread('../mesh_x.dat');
mesh_y = dlmread('../mesh_y.dat');

figure(1),clf,hold on
for j = 1:size(mesh_x,1)
% for j =1
  plot(mesh_x(j,:),mesh_y(j,:),'k')
end

for i = 1:size(mesh_x,2)-1
  plot(mesh_x(:,i),mesh_y(:,i),'m')
end

% xlim([-.1,1.1])

axis equal
grid on
xlabel('x')
ylabel('y')

% figure(2),clf
% plot mesh in the transformed domain?
