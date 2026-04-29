function f_plot_mesh(fignum,address,casename,filename)

% mesh_x = readmatrix([address,casename,'_initial_x.msh'],'filetype','delimitedtext');
% mesh_x = readmatrix([address,casename,'_initial_y.msh'],'filetype','delimitedtext');

mesh_x = dlmread([address,casename,filename{1}]);
mesh_y = dlmread([address,casename,filename{2}]);

figure(fignum),clf,hold on
for j = 1:size(mesh_x,1)
% for j = 1
  plot(mesh_x(j,:),mesh_y(j,:),'k')
end

for i = 1:size(mesh_x,2)-1
% for i = [1,2,size(mesh_x,2)-1,size(mesh_x,2)-2]
  plot(mesh_x(:,i),mesh_y(:,i),'m')
end

% xlim([-.1,1.1])

axis equal
grid on
xlabel('x')
ylabel('y')

% figure(2),clf
% plot mesh in the transformed domain?

end