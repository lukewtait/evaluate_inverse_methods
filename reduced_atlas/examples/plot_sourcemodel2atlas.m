function plot_sourcemodel2atlas_glassbrain(sourcemodel,vol)

clf
ft_plot_mesh(vol.bnd(1),'facecolor','brain','edgealpha',0.05,'facealpha',0.4)
ft_plot_mesh(vol.bnd(2),'facecolor',0.8*ones(1,3),'edgealpha',0.05,'facealpha',0.4)
ft_plot_mesh(vol.bnd(3),'facecolor','skin','edgealpha',0.05,'facealpha',0.4)
hold on

if ~isfield(sourcemodel,'tri')
    scatter3(sourcemodel.pos(:,1),sourcemodel.pos(:,2),sourcemodel.pos(:,3),5,...
        sourcemodel.tissue,'filled')
else
    patch('vertices',sourcemodel.pos,'faces',sourcemodel.tri,'facecolor','flat',...
        'facevertexCdata',sourcemodel.tissue)
end

rng('default') ; colormap(repmat(rand(115,3),2,1)) ;