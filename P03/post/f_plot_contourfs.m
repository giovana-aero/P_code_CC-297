function f_plot_contourfs(mesh,data,config)

[~,h] = contourf(mesh.x,mesh.y,data,20);
set(h,'LineColor','none')
axis equal

colormap(slanCM('bupu'))

cb_latex(colorbar(),config.cb_label,config.fontsize);

end