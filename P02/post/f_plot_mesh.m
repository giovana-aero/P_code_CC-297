function f_plot_mesh(fignum,address,casename,filename,linestyle)

if nargin == 4
  linestyle = '-';
end

cmap = slanCM('bupu');

% mesh_x = readmatrix([address,casename,'_initial_x.msh'],'filetype','delimitedtext');
% mesh_x = readmatrix([address,casename,'_initial_y.msh'],'filetype','delimitedtext');

disp(filename{1})
disp([address,casename,filename{1}])
mesh_x = dlmread([address,casename,filename{1}]);
mesh_y = dlmread([address,casename,filename{2}]);

figure(fignum),hold on
for j = 1:size(mesh_x,1)
% for j = 1
  plot(mesh_x(j,:),mesh_y(j,:),['k',linestyle],'linewidth',.8)
end

for i = 1:size(mesh_x,2)
% for i = [1,2,size(mesh_x,2)-1,size(mesh_x,2)-2]
  plot(mesh_x(:,i),mesh_y(:,i),'color',cmap(150,:),'linestyle',linestyle,'linewidth',.8)
end

% xlim([-.1,1.1])

axis equal
grid on
xlabel('x')
ylabel('y')

% figure(2),clf
% plot mesh in the transformed domain?

end