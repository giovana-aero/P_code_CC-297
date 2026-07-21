function f_plot_contours(mesh,data,config)

% contour(mesh.x,mesh.y,data,20,'showtext','on','labelformat','%0.2f')
[C,h] = contour(mesh.x,mesh.y,data,20,'labelformat','%0.2f','labelspacing',100);
axis equal

clabel(C,h,'interpreter','latex','fontsize',config.fontsize)

colormap(slanCM('bupu'))

cb_latex(colorbar(),config.cb_label,config.fontsize);

end