% Analysis of orientations from electron backscatter diffraction (EBSD)
% indexing of five phases from an Al-steel joint.
%
% Håkon Wiik Ånes (hakon.w.anes@ntnu.no) and Tina Bergh
% (tina.bergh@ntnu.no)
% 2022-05-31
%
% Requires MTEX and export_fig

clear variables
close all

% MTEX configuration
plotx2east
plotzIntoPlane

% export_fig configuration
res = '-r200';

% Crystal and specimen symmetry
colors = {[0, 143, 159], [117, 216, 185], [127, 198, 130],...
    [246, 122, 114], [253, 225, 4]};
cs = {'notIndexed',...
    crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al'),...
    crystalSymmetry('m-3', [12.56, 12.56, 12.56], 'mineral', 'alpha_alfesi'),...
    crystalSymmetry('m-3m', [2.87, 2.87, 2.87], 'mineral', 'ferrite'),...
    crystalSymmetry('mmm', [7.66, 6.41, 4.23], 'mineral', 'fe2al5'),...
    crystalSymmetry('2/m', [15.01, 8.07, 12.47], [90, 107.720, 90] * degree, 'mineral', 'fe4al13')
};
for i=2:length(cs)
    cs{i}.color = colors{i - 1} / 255;
end

% Directory and file names
% Dataset naming (a-c) = (I-III)
dset = 'a';
dir_data = fullfile('/home/hakon/phd/data/tina', dset, 'kp/merged');
dir_mtex = fullfile(dir_data, 'mtex');
fname = 'xmap_merged.ang';

% Read orientation data
ebsd = EBSD.load(fullfile(dir_data, fname), cs, 'columnNames', ...
    {'phi1', 'Phi', 'phi2', 'x', 'y', 'iq', 'ci', 'Phase',...
    'detector_signal', 'fit', 'mean_intensity'},...
    'radians');
ebsd.prop = rmfield(ebsd.prop, {'detector_signal', 'iq', 'fit'});
rot_tsl2mtex = rotation.byAxisAngle(xvector - yvector, 180 * degree);
ebsd = rotate(ebsd, rot_tsl2mtex, 'keepXY');

% Mineral list
phase_names = ebsd.mineralList;

% Dataset specific parameters
mean_intensity_threshold = struct('a', 0.35, 'b', 0.2, 'c', 0);

%% Threshold data
ebsd2 = ebsd;
condition1 = ebsd2.ci < 0.1;
condition2 = ebsd2.mean_intensity < mean_intensity_threshold.(dset);
ebsd2(condition1 | condition2).phase = 0;

%% Relative phase fractions of intermetallic phases
ids = [2, 4, 5];
imp_frac = zeros(3, 1);
for i=1:length(ids)
    imp_frac(i) = sum(ebsd2.phaseId == ids(i) + 1);
end
for i=1:length(ids)
    phase = ebsd2.CSList(ids(i) + 1);
    fprintf('%s: %.4f\n', phase{1}.mineral, imp_frac(i) / sum(imp_frac))
end

%% Reconstruct grains
[grains, ebsd2.grainId, ebsd2.mis2mean] = calcGrains(ebsd2, 'angle',...
    5 * degree);

%% Phase map with boundaries
figure
[~, mP] = plot(ebsd2, 'micronbar', 'off');
legend('hide')
hold on
for i=2:length(phase_names)
    phase_name = phase_names{i};
    for j=i + 1:length(phase_names)
        other_phase_name = phase_names{j};
        plot(grains.boundary(phase_name, other_phase_name), 'linewidth', 3)
        hold on
    end
end
export_fig(fullfile(dir_mtex, 'maps_phase_no_scalebar.png'), res)

%% Orientation color keys
ckeys = [];
for i=2:length(phase_names)
    ckey = ipfHSVKey(cs{i}.properGroup);
    ckeys = [ckeys ckey];

    % Plot color keys
    figure
    plot(ckey)
    export_fig(fullfile(dir_mtex, ['ipf_key_' phase_names{i} '.png']), res)
end

%% Orientation color maps
directions = [xvector, yvector];
labels = {'x', 'y'};
for i=1:length(directions)
    figure
    for j=2:length(phase_names)
        ckey = ckeys(j - 1);
        ckey.inversePoleFigureDirection = directions(i);
        rgb = ckey.orientation2color(ebsd2(phase_names{j}).orientations);
        plot(ebsd2(phase_names{j}), rgb, 'micronbar', 'off')
        hold on
    end
    export_fig(fullfile(dir_mtex, ['maps_ipf_' labels{i} '.png']), res)
end

%% Orientations of eta (Fe2Al5) grains
figure
plotIPDF(ebsd2('fe2al5').orientations, [xvector, yvector, zvector],...
    'markersize', 10, 'contourf')
export_fig(fullfile(dir_mtex, 'ipf_pole_density_function_eta.png'), res)

%% Misorientations Al - Alpha-AlFeSi (Alpha_C)
gb_al_alfesi = grains.boundary('al', 'alpha_alfesi');
mori_al_alfesi = gb_al_alfesi.misorientation;

%% Plot highlighted phase boundary
ori_colors = cell(5, 1);
figure
for i=2:length(phase_names)
    ckey = ckeys(i - 1);
    ckey.inversePoleFigureDirection = yvector;
    rgb = ckey.orientation2color(ebsd2(phase_names{i}).orientations);
    ori_colors{i - 1} = rgb;
    plot(ebsd2(phase_names{i}), rgb, 'facealpha', 0.5)
    hold on
end
plot(gb_al_alfesi, 'linewidth', 3)
hold off
export_fig(fullfile(dir_mtex, 'maps_al_alfesi_boundary.png'), res)

%% Calculate misorientation clusters
[c_id, center_mori] = calcCluster(mori_al_alfesi);

c_id_unique = unique(c_id);
n_clusters = length(c_id_unique);

%% Specify orientation relationships
vectors_al = [...
    [1, -1,  1;  0, 1,  1],... % 1 - Muggerud
    [2,  0,  2;  1, 1, -1],... % 2
    [1,  0,  1;  0, 2,  0],... % 3
    [1,  0,  1; -1, 1,  1],... % 4
    [1,  1,  1;  2, 0, -2],... % 5
    [1,  1,  1;  2, 0, -2],... % 6
    [1,  0,  5;  0, 2,  0],... % 7
    [1,  1, -1; -1, 0, -1],... % 8 - Mackay
    [1,  0,  0;  0, 0,  1]...  % 9
];
vectors_alfesi = [...
    [1, -1,  1;  5, -2, -7],... % 1 - Muggerud
    [5,  2, -7;  1,  1,  1],... % 2
    [1,  2,  3;  5,  2, -3],... % 3
    [2,  0,  5; -5,  1,  2],... % 4
    [1,  1,  3; -8,  2,  2],... % 5
    [1,  0,  0;  0,  0,  1],... % 6
    [1,  0,  1;  0,  2,  0],... % 7
    [4,  1,  1; -2, -1,  9],... % 8 - Mackay
    [1,  0,  0;  0,  0,  1]...  % 9
];
n_or = length(vectors_al) / 3;
vectors_al = reshape(vectors_al, [2, 3, n_or]);
vectors_alfesi = reshape(vectors_alfesi, [2, 3, n_or]);

or_al_alfesi = [];
for i=1:n_or
    uvw1_i = vectors_al(1, :, i);
    uvw2_i = vectors_alfesi(1, :, i);
    hkl1_i = vectors_al(2, :, i);
    hkl2_i = vectors_alfesi(2, :, i);
    or_i = orientation.map(...
        Miller(uvw1_i(1), uvw1_i(2), uvw1_i(3), 'uvw', cs{2}),...
        Miller(uvw2_i(1), uvw2_i(2), uvw2_i(3), 'uvw', cs{3}),...
        Miller(hkl1_i(1), hkl1_i(2), hkl1_i(3), 'hkl', cs{2}),...
        Miller(hkl2_i(1), hkl2_i(2), hkl2_i(3), 'hkl', cs{3})...
    );
    or_al_alfesi = [or_al_alfesi or_i];
end

%% Distance to orientation relationships
for i=1:(n_clusters - 1)
    angle(center_mori(i), or_al_alfesi) / degree
end

%% Axis-angle representation of misorientation cluster centers
aa_axis = round(center_mori.axis);
aa_angle = center_mori.angle / degree;
for i=1:(n_clusters - 1)
    fprintf('Cluster %i:', i)
    aa_axis(i).uvw
    aa_angle(i)
end

%% Max. and standard deviation of misorientation clusters
for i=2:n_clusters
    mask_i = (c_id == c_id_unique(i));
    angles = angle(center_mori(i - 1), mori_al_alfesi(mask_i)) / degree;
    fprintf('Cluster %i has max. dev. %.5f and std. dev %.5f\n',...
        i - 1, max(angles), std(angles))
end

%% Plot clusters in axis-angle space
cluster_color_names = colornames('Alphabet');

figure
for i=1:n_clusters
    if i == 1
        name_i = {'Black'};
        rgb_i = 'k';
    else    
        [name_i, rgb_i] = colornames('Alphabet', cluster_color_names(i));
    end
    fprintf('Cluster %i is colored %s\n', c_id_unique(i), name_i{1})
    mask_i = (c_id == c_id_unique(i));
    mP = plot(mori_al_alfesi(mask_i), 'markeredgecolor', rgb_i);
    hold on
end

% Add orientation relationships from theory
or_color_names = colornames('Tableau');
for i=1:n_or
    [~, rgb_i] = colornames('Tableau', or_color_names(i));
    plot(or_al_alfesi(i), 'markersize', 20, 'markerfacecolor', rgb_i,...
        'markeredgecolor', 'k')
end

export_fig(fullfile(dir_mtex, 'or_al_alfesi_axis_angle.png'), res)
views = [0, 0; 0, 90; 90, 0; 90, 90];
for i=1:length(views)
    view(views(i, 1), views(i, 2))
    export_fig(fullfile(dir_mtex,...
        ['or_al_alfesi_axis_angle_az' num2str(views(i, 1)) '_el'...
        num2str(views(i, 2)) '.png']), res)
end

%% Plot clusters in map
figure
for i=2:length(phase_names)
    rgb = ori_colors{i - 1};
    plot(ebsd2(phase_names{i}), rgb, 'facealpha', 0.3)
    hold on
end

for i=1:n_clusters
    if i == 1
        rgb_i = 'k';
    else
        [~, rgb_i] = colornames('Alphabet', cluster_color_names(i));
    end
    mask_i = (c_id == c_id_unique(i));
    plot(gb_al_alfesi(mask_i), 'linecolor', rgb_i, 'linewidth', 3)
    hold on
end

export_fig(fullfile(dir_mtex, 'maps_al_alfesi_mori_clusters.png'), res)