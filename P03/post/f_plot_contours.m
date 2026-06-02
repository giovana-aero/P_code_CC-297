function f_plot_contours(mesh,data,config)

[~,h] = contourf(mesh.x,mesh.y,data);
grid on
set(h,'LineColor','none')
axis equal

colormap(slanCM('bupu'))

cb_latex(colorbar(),config.cb_label,config.fontsize);

end