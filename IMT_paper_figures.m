function IMT_paper_figures(fname)

%% make a bunch of figures for the unthresholded weighmap and the threhsolded bootstrapped weightmap

% run from the task results folder (Box/@Predict_IMT/IAPS or Box/@Predict_IMT/Faces/OF/AllFaces)
% input .nii file from which we want figures (e.g., weightmap_unthresholded.nii)

if nargin < 1, fname = 'weightmap_unthresholded.nii'; end

%%
fstring = fname(1:end-4);

%% set colors
% requires 'RGB XKCD' toolbox
% https://www.mathworks.com/matlabcentral/fileexchange/46872-intuitive-rgb-color-values-from-xkcd
negcm = colormap_tor(rgb('blue'), rgb('light blue'));
poscm = colormap_tor(rgb('orange'), rgb('yellow'));
%  pos_colormap = colormap_tor(minposcolor, maxposcolor);
            %neg_colormap = colormap_tor(minnegcolor, maxnegcolor);
fourcolors = {rgb('blue'), rgb('light blue'), rgb('orange'), rgb('yellow')};
% split color is 
%   {minnegcolor maxnegcolor minposcolor maxposcolor}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load map 
w = fmri_data(fname)

%% full montage
% close all; clear o2; o2 = montage(region(w), 'montagetype', 'full', 'splitcolor', fourcolors)
% g = gcf; delete(g.Children(1)) % remove legend
% print('-dpng', '-r500', [fstring '_montage_full.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% four surfaces - these show some brainstem
ss = get(0, 'ScreenSize');
close all; clear o2;  o2 = fmridisplay; figure('Position', [0 0 ss(3)*.9 ss(3)*.225]); axis image; axis off;
%o2 = montage(o2, 'nofigure', 'onerow')
o2 = surface(o2, 'axes', [-0.25 0.1 .75 .75],'direction', 'hires right', 'orientation', 'medial', 'brainstem')
o2 = surface(o2, 'axes', [0.0 0.1 .75 .75], 'direction', 'hires left', 'orientation', 'medial', 'brainstem')
o2 = surface(o2, 'axes', [0.25 0.1 .75 .75], 'direction', 'hires left', 'brainstem')
o2 = surface(o2, 'axes', [0.5 0.1 .75 .75], 'direction', 'hires right', 'brainstem')
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors)
print('-dpng', '-r500', [fstring '_foursurfaces_onerow_brainstem.png'])

%% four surfaces - no brainstem
ss = get(0, 'ScreenSize');
close all; clear o2;  o2 = fmridisplay; figure('Position', [0 0 ss(3)*.9 ss(3)*.225]); axis image; axis off;
%o2 = montage(o2, 'nofigure', 'onerow')
o2 = surface(o2, 'axes', [-0.25 0.1 .75 .75],'direction', 'hires right', 'orientation', 'medial')
o2 = surface(o2, 'axes', [0.0 0.1 .75 .75], 'direction', 'hires left', 'orientation', 'medial')
o2 = surface(o2, 'axes', [0.25 0.1 .75 .75], 'direction', 'hires left')
o2 = surface(o2, 'axes', [0.5 0.1 .75 .75], 'direction', 'hires right')
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors)
print('-dpng', '-r500', [fstring '_foursurfaces_onerow.png'])


% OLD VERSION - colors don't match montages well
% close all; 
% [all_surf_handles, pcl, ncl] = surface(region(w), 'foursurfaces', 'pos_colormap', poscm, 'neg_colormap', negcm)
% print('-dpng', '-r500', [fstring '_foursurfaces.png'])

%[all_surf_handles, pcl, ncl] = surface(region(w), 'foursurfaces', 'splitcolor', fourcolors)

%% limbic and brainstem (also from inflammation meta-analysis

close all; clear o2;  o2 = fmridisplay; figure('Position', [0 0 700 400]); axis image; axis off;
o2 = surface(o2, 'axes', [-0.1 0.1 .75 .75], 'direction', 'limbic right')
o2 = surface(o2, 'axes', [0.4 0.1 .75 .75], 'direction', 'limbic left')
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors)
print('-dpng', '-r500', [fstring '_limbic.png'])


% % OLD VERSION - colors don't match montages well
% close all; create_figure; p = addbrain('limbic'); delete(p(6)); p(6)=[];
% [all_surf_handles, pcl, ncl] = surface(region(w), 'pos_colormap', poscm, 'neg_colormap', negcm, 'surface_handles', p);
% g = gcf; 
% set(g.Children, 'View', [95 5]); lightRestoreSingle;
% print('-dpng', '-r500', [fstring '_limbic_right.png'])
% set(g.Children, 'View', [-95 5]); lightRestoreSingle;
%print('-dpng', '-r500', [fstring '_limbic_left.png'])

%% montage of selected axial slices
zvals = [-24 -16 -8]; % pick your z-values here
axial_slices = zeros(length(zvals), 3);
axial_slices(:, 3) = zvals; % z values for coordinates

close all; clear o2; o2 = fmridisplay;         
o2 = montage(o2, 'axial', 'onerow', 'wh_slice', axial_slices);

axh = axes('Position', [0.0100 0.1100 0.1220 1.1410]);
o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors);
locs = axial_slices(:, 3);
for i = 1:length(locs)
    %draw_vertical_line(locs(i));
    hh(i) = plot([-105 65], [locs(i) locs(i)], 'b', 'LineWidth', 1);
end
print('-dpng', '-r500', [fstring '_montage_axial_selected_slices' sprintf('_%d', zvals) '.png'])


%% sagittal midline slice
close all; clear o2; o2 = fmridisplay;
axh = axes('Position', [0.1 0.1 .8 .8]);
o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors);
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors);
print('-dpng', '-r500', [fstring '_sagittal_midline_nolines.png']);


%% sagittal midline slice - with horizontal blue lines
close all; clear o2; o2 = fmridisplay;
axh = axes('Position', [0.1 0.1 .8 .8]);
o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors);
o2 = addblobs(o2, region(w), 'splitcolor', fourcolors);
locs = axial_slices(:, 3);
for i = 1:length(locs)
    %draw_vertical_line(locs(i));
    hh(i) = plot([-105 65], [locs(i) locs(i)], 'b', 'LineWidth', 3);
end
print('-dpng', '-r500', [fstring '_sagittal_midline_withlines.png']);


%%
return;



%% big figure from inflammation meta-analysis
% THIS COLORBAR DOES NOT MATCH OTHER FIGURES 
% close all; create_figure('surface');
% [all_surf_handles, pcl, ncl] = surface(region(w), 'cutaway', 'ycut_mm', -50, 'pos_colormap', poscm, 'neg_colormap', negcm, 'existingfig');
% print('-dpng', '-r500', [fstring '_big_cutaway.png'])




%% montage of selected coronal slices
coronal_slices = zeros(3, 3);
coronal_slices(:, 2) = [-8     0     8   ]; % y values for slices

close all; clear o2; o2 = fmridisplay;         
o2 = montage(o2, 'coronal', 'onerow', 'wh_slice', coronal_slices, 'text');




%% sagittal slices for insula
close all; clear o2; o2=fmridisplay;
o2 = montage(o2, 'sag', 'slice_range', [-43 -32], 'spacing', 1)
o2 = addblobs(o2, r)

close all; clear o2; o2=fmridisplay;
o2 = montage(o2, 'sag', 'slice_range', [32 43], 'spacing', 1)
o2 = addblobs(o2, r)



%% other stuff below

close all; clear o2; o2 = montage(w, 'splitcolor', {rgb('blue') rgb('light blue') rgb('orange') rgb('yellow')})


close all; clear o2
[all_surf_handles, pcl, ncl] = surface(region(w), 'hires')



% mediation figure
[all_surf_handles, pcl, ncl] = surface(region(w), 'hires')



close all; clear o2; o2=fmridisplay;
o2 = montage(o2, 'sag', 'slice_range', [-43 -32], 'spacing', 1)
o2 = addblobs(o2, r)

close all; clear o2; o2=fmridisplay;
o2 = montage(o2, 'sag', 'slice_range', [32 43], 'spacing', 1)
o2 = addblobs(o2, r)


[all_surf_handles, pcl, ncl] = surface(region(w), 'hires')





[all_surf_handles, pcl, ncl] = surface(region(w), 'foursurfaces', 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]})


create_figure; p = addbrain('limbic'); delete(p(6)); p(6)=[];
[all_surf_handles, pcl, ncl] = surface(region(w), 'pos_colormap', poscm, 'neg_colormap', negcm, 'surface_handles', p)