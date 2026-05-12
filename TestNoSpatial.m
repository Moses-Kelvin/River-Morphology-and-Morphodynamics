%% RivMAP Analysis (Multiple Years - Fully Corrected & Robust)
% This script analyzes multi-year river channel masks using the RivMAP toolbox.
% It computes centerlines, banklines, widths, sinuosity, curvature, cutoffs,
% migration rates, and spatial-temporal patterns.
%
% Requirements:
%   - RivMAP toolbox in MATLAB path.
%   - Mask images named 'maskYYYY.tif' (e.g., mask2000.tif, mask2005.tif).
%   - Pixel size (meters per pixel) must be set correctly.
%   - Signal Processing Toolbox (for findpeaks).
%
% Author: [Your Name]
% Date: [Current Date]

close all; clear; clc;

% Add RivMAP path (adjust if necessary)
addpath('RivMAP');          % adds the RivMAP folder to the path

%% User Settings
pixel_size = 30;            % meters per pixel (adjust to your image resolution)
% Options
plot_individual = true;     % set to true to display each year's masks
plot_width_profile = false;  % set to true to save width profile figures
plot_migration_maps = true; % set to true to show migration maps

%% 1. Load Masks Automatically
files = dir('mask_KJ*.tif');
n = length(files);
if n == 0
    error('No mask files found');
end

% Sort files by year (assuming filename contains year)
years = zeros(n,1);
for i = 1:n
    tokens = regexp(files(i).name, '\d+', 'match');
    if ~isempty(tokens)
        years(i) = str2double(tokens{1});
    else
        error('Could not extract year from filename: %s', files(i).name);
    end
end
[~, idx] = sort(years);
files = files(idx);
years = years(idx);

es = 'NS';         % adjust if river flows north‑south (options: 'WE' or 'NS')
plotornot = 0;     % suppress internal RivMAP plots (set to 1 if needed)

% Preallocate structure array
riv = struct('meta', struct(), 'im', struct(), 'vec', struct(), 'mig', struct());

for i = 1:n
    fname = files(i).name;
    I = imread(fname) > 0;   % binary mask
    
    % Basic cleaning: fill holes, remove small isolated islands (adjust threshold)
    I = imfill(I, 'holes');
    I = bwareaopen(I, 20);   % remove objects <100 pixels

    I = imclose(I, strel('disk',3));
    I = imopen(I, strel('disk',2));
    
    % Store original full channel mask
    riv(i).im.hc = I;
    
    % Braiding index (number of separate channel threads)
    CC = bwconncomp(I);
    stats = regionprops(CC, 'Area');

    areas = [stats.Area];

    % Keep only significant channels (threshold)
    % threshold = 0.001 * sum(I(:));   % relative to river area
    threshold = 500;   % fixed pixel threshold
    braiding_index = sum(areas > threshold);
    riv(i).vec.braiding_index = braiding_index;
    
    % For centerline extraction, use the largest connected component (main channel)
    % This may lose braiding information but is necessary for a single centerline.
   I_main = bwareafilt(I, 1);   % largest component only
   I_main = imfill(I_main, 'holes');
    I_main = imclose(I_main, strel('disk',2));
    riv(i).im.st = I_main;
    
    riv(i).meta.year = years(i);
    riv(i).meta.exit_sides = es;
end

fprintf('Loaded %d mask files (years: %s).\n', n, num2str(years', ' %d'));

%% 2. Visualize Masks (Optional)
if plot_individual
    figure('Name', 'Loaded Masks');
    for i = 1:n
        imshow(riv(i).im.st);
        title(sprintf('Year %d', riv(i).meta.year));
        pause(0.2);
    end
end

%% 3. Per-Year Analysis (Centerline, Banks, Width, Curvature, etc.)
for i = 1:n
    fprintf('Processing year %d...\n', riv(i).meta.year);
    
    I = riv(i).im.st;          % logical mask (main channel)
    I_double = double(I);      % convert to double for RivMAP functions
    es = riv(i).meta.exit_sides;
    
    % --- Centerline extraction ---
    % Initial estimate of nominal width Wn
    D = bwdist(~I);
    width_map = 2 * D;
    channel_widths = width_map(I);
    channel_widths = channel_widths(channel_widths > 0);
    channel_widths = channel_widths(~isnan(channel_widths));
    Wmed = median(channel_widths);
    Wn_initial = round(0.4 * Wmed);
    Wn_initial = max(round(0.25 * Wmed), min(Wn_initial, 80));
    Wn_initial = double(Wn_initial);   % ensure double for bwareaopen
    
    % Extract centerline (using double mask and double threshold)
    [cl, Icl] = centerline_from_mask(I_double, es, Wn_initial, plotornot);
    if isempty(cl)
        error('Centerline extraction failed for year %d', riv(i).meta.year);
    end
    
    % --- Banklines ---
    banks = banklines_from_mask(I_double, es, plotornot);
    if length(banks) < 2
        error('Bankline extraction failed for year %d', riv(i).meta.year);
    end
    lb = banks{1};
    rb = banks{2};
    
    % --- Width from banklines ---
    W = width_from_banklines(cl, lb, rb, Wn_initial);
    W_clean = W(~isnan(W));
    W_clean = W_clean(W_clean > 0);
    Wavg = mean(W_clean);
   % ✅ Convert to meters
   W = W * pixel_size;
   Wavg = Wavg * pixel_size;
   riv(i).meta.Wavg = Wavg;
   riv(i).meta.Wn = Wn_initial;
    
    % --- Streamwise distance ---
    S = [0; cumsum(sqrt(diff(cl(:,1)).^2 + diff(cl(:,2)).^2))] * pixel_size;
    riv(i).vec.S = S;
    riv(i).vec.cl = cl;
    riv(i).vec.lb = lb;
    riv(i).vec.rb = rb;
    riv(i).vec.W = W;
    riv(i).vec.Wavg = Wavg;
    
    % --- Width from mask (optional) ---
    try
        [Wm, SWm] = width_from_mask(I_double, cl, round(Wavg / pixel_size));
        Wm = Wm * pixel_size;
        riv(i).vec.Wm = Wm;
        SWm = SWm * pixel_size;
        riv(i).vec.SWm = SWm;
    catch ME
        warning('width_from_mask failed for year %d: %s', riv(i).meta.year, ME.message);
        riv(i).vec.Wm = [];
        riv(i).vec.SWm = [];
    end
    
    % --- Smooth centerline for curvature ---
    % cls = savfilt(cl, round(1.5 * Wn_initial));
    cls = savfilt(cl, round(Wn_initial));
    % cls = savfilt(cl, Wn_initial);
    riv(i).vec.cls = cls;
    
    % --- Curvature and angles ---
    A = angles(cls);
    C = curvatures(cls);
    riv(i).vec.A = A;
    riv(i).vec.C = C;
    
    % --- Sinuosity ---
    valley_len = sqrt((cl(end,1)-cl(1,1))^2 + (cl(end,2)-cl(1,2))^2) * pixel_size;
    channel_len = S(end);
    sinuosity = channel_len / valley_len;
    riv(i).vec.sinuosity = sinuosity;
    riv(i).vec.cl_len = channel_len;
    
    % --- Wavelength from curvature peaks (using movmean instead of smooth) ---
    window = round(Wn_initial/2);
    if mod(window,2) == 0
        window = window + 1;          % make odd for symmetric smoothing
    end
    C_smooth = movmean(C, window, 'Endpoints', 'shrink');
    [pk, loc] = findpeaks(abs(C_smooth), 'MinPeakHeight', 0.2*max(abs(C_smooth)));
    if length(loc) > 1
        wavelength = mean(diff(S(loc)));
    else
        wavelength = NaN;
    end
    riv(i).vec.wavelength = wavelength;
    
    % --- Average absolute curvature ---
    riv(i).vec.avg_curvature = mean(abs(C(isfinite(C))));
    
    % --- Reach-averaged width ---
    riv(i).vec.Wra = (sum(I(:)) * pixel_size^2) / S(end);
    
    % --- (Optional) Plot width profile ---
    if plot_width_profile
        figure('Name', sprintf('Width Profile %d', riv(i).meta.year));
        plot(S, W, 'b-', 'LineWidth', 1.5); hold on;
        if ~isempty(riv(i).vec.Wm)
            plot(SWm, Wm, 'r-', 'LineWidth', 1.5);
            legend('Banklines', 'Mask', 'Location', 'best');
        else
            legend('Banklines', 'Location', 'best');
        end
        xlabel('Streamwise distance (m)');
        ylabel('Width (m)');
        title(sprintf('Year %d - Width Profile', riv(i).meta.year));
        grid on;
        saveas(gcf, sprintf('width_profile_%d.png', riv(i).meta.year));
        close(gcf);
    end
end
fprintf('Per‑year analysis completed.\n');

%% 4. Multi-Year Visualizations (Optional)
if plot_individual
    figure('Name', 'All Years Overlay');
    for i = 1:n
        imshow(riv(i).im.st); hold on;
        plot(riv(i).vec.cls(:,1), riv(i).vec.cls(:,2), 'b', 'LineWidth', 1.5);
        plot(riv(i).vec.lb(:,1), riv(i).vec.lb(:,2), 'm', 'LineWidth', 1.5);
        plot(riv(i).vec.rb(:,1), riv(i).vec.rb(:,2), 'm', 'LineWidth', 1.5);
        title(sprintf('Year %d', riv(i).meta.year));
        pause(0.3);
    end
end

%% 5. Migration and Cutoff Detection (Between Consecutive Years)
% Create a separate structure for migration intervals
migration = struct('year1', [], 'year2', [], 'Imig', [], 'Ie', [], 'Ia', [], ...
                   'cutarea', [], 'cutlen', [], 'chutelen', [], 'Icuts', []);

% For each consecutive pair of years
for i = 1:(n-1)
    fprintf('Computing migration from %d to %d...\n', riv(i).meta.year, riv(i+1).meta.year);
    
    % Get centerlines (smoothed)
    cl1 = riv(i).vec.cls;
    cl2 = riv(i+1).vec.cls;
    es1 = riv(i).meta.exit_sides;
    es2 = riv(i+1).meta.exit_sides;
    
    % Use the later year's Wn for migration functions (often recommended)
    Wn = riv(i+1).meta.Wn;
    sizeI = size(riv(i).im.st);
    
    % Centerline-based migration and cutoffs
    [Imig, Icutsall, cutidcs, cutarea, cutlen, chutelen] = ...
        migration_cl(cl1, cl2, es1, es2, Wn, sizeI);
    
    % Mask-based erosion/accretion
    I1 = riv(i).im.hc;   % full channel mask (original, not main component)
    I2 = riv(i+1).im.hc;
    [Ie, Ia, Inc, Icutsm] = migration_mask(I1, I2, Wn);
    
    % Store in migration structure
    migration(i).year1 = riv(i).meta.year;
    migration(i).year2 = riv(i+1).meta.year;
    migration(i).Imig = Imig;
    migration(i).Ie = Ie;
    migration(i).Ia = Ia;
    migration(i).Icuts = Icutsall;
    migration(i).cutarea = cutarea;
    migration(i).cutlen = cutlen;
    migration(i).chutelen = chutelen;
    
    % Compute migration areas
    migration(i).MAcl = sum(Imig(:));
    migration(i).MAe = sum(Ie(:));
    migration(i).MAa = sum(Ia(:));
    
    % For reach-averaged migration rates, we need the centerline length of the earlier year
    migration(i).cl_len = riv(i).vec.cl_len;
end

% Display cutoff summary
fprintf('Migration intervals computed: %d\n', length(migration));
for i = 1:length(migration)
    if ~isempty(migration(i).cutarea)
        fprintf('Interval %d–%d: %d cutoffs, total area %.2f pixels\n', ...
            migration(i).year1, migration(i).year2, length(migration(i).cutarea), sum(migration(i).cutarea));
    else
        fprintf('Interval %d–%d: no cutoffs\n', migration(i).year1, migration(i).year2);
    end
end

%% 6. Migration Maps (Accumulated over All Intervals)
Imig_total = false(size(riv(1).im.st));
Ie_total = false(size(riv(1).im.st));
Ia_total = false(size(riv(1).im.st));

for i = 1:length(migration)
    Imig_total = Imig_total | migration(i).Imig;
    Ie_total = Ie_total | migration(i).Ie;
    Ia_total = Ia_total | migration(i).Ia;
end

if plot_migration_maps
    figure('Name', 'Accumulated Migration Maps');
    subplot(1,3,1); imshow(Imig_total); title('Total Migration Area');
    subplot(1,3,2); imshow(Ie_total); title('Total Erosion');
    subplot(1,3,3); imshow(Ia_total); title('Total Accretion');
end

%% 7. Reach-Averaged Migration Rates (per Interval)
% Compute migration rates (pixels per year) using centerline length of earlier year
% Convert to meters per year using pixel_size
for i = 1:length(migration)
    dt = migration(i).year2 - migration(i).year1;
    if dt <= 0
        error('Non‑positive time interval at index %d', i);
    end
   % migration(i).Mrcl = (migration(i).MAcl / migration(i).cl_len) * pixel_size / dt;
    % Corrected migration rates (meters per year)
    migration(i).Mrcl = (migration(i).MAcl / migration(i).cl_len) * pixel_size / dt;
    migration(i).Mre  = (migration(i).MAe  / migration(i).cl_len) * pixel_size / dt;
    migration(i).Mra  = (migration(i).MAa  / migration(i).cl_len) * pixel_size / dt;
end

%% 7b. TOTAL (CUMULATIVE) MIGRATION

total_migration = 0;

for i = 1:length(migration)
    dt = migration(i).year2 - migration(i).year1;
    total_migration = total_migration + migration(i).Mrcl * dt;
end

fprintf('Total cumulative migration (m): %.2f\n', total_migration);
%% 7c. NET MIGRATION (START → END)

cl_start = riv(1).vec.cls;
cl_end   = riv(end).vec.cls;

Dmat = pdist2(cl_start, cl_end);
minDist = min(Dmat, [], 2);

net_migration = mean(minDist) * pixel_size;

fprintf('Net migration (1985–2025): %.2f m\n', net_migration);

% Plot reach-averaged migration rates
years_mig = arrayfun(@(x) x.year2, migration);
mig_cl = arrayfun(@(x) x.Mrcl, migration);
mig_erosion = arrayfun(@(x) x.Mre, migration);
figure('Name', 'Reach‑Averaged Migration Rates');
plot(years_mig, mig_cl, 'ro-', 'LineWidth', 2); hold on;
plot(years_mig, mig_erosion, 'bo-', 'LineWidth', 2);
legend('Centerline Migration', 'Erosion', 'Location', 'best');
xlabel('Year (end of interval)'); ylabel('Migration Rate (m/yr)');
title('Reach‑Averaged Migration Rate');
grid on;

%% 8. Cutoff Analysis
% Accumulate all cutoffs and create a colored map
Icutoffs_total = false(size(riv(1).im.st));
cut_years = [];
cut_areas = [];

for i = 1:length(migration)
    if ~isempty(migration(i).Icuts)
        Icutoffs_total = Icutoffs_total | migration(i).Icuts;
        ncuts = length(migration(i).cutarea);
        cut_years = [cut_years, repmat(migration(i).year2, 1, ncuts)];
        cut_areas = [cut_areas, migration(i).cutarea(:)'];
    end
end

% Figure: Binary cutoff map, color‑by‑year, and bar chart
figure('Name', 'Cutoff Analysis');
subplot(1,3,1); imshow(Icutoffs_total); title('Binary map of cutoffs','FontSize',14);

% Color‑by‑year
Icutc = zeros(size(riv(1).im.st));
for i = 1:length(migration)
    if ~isempty(migration(i).Icuts)
        Icutc(migration(i).Icuts) = i;   % assign interval index
    end
end
subplot(1,3,2);
cmap = parula(length(migration));
imshow(Icutc, cmap);
cb = colorbar;
% Set colorbar ticks to show years
ytick_pos = linspace(1, length(migration), 5);
ytick_labels = round(interp1(1:length(migration), years_mig, ytick_pos));
set(cb, 'Ticks', ytick_pos/length(migration), 'TickLabels', ytick_labels);
title('Cutoffs colored by interval','FontSize',14);

% Bar chart of cutoff areas over time
subplot(1,3,3);
if ~isempty(cut_areas)
    bar(cut_years, cut_areas * pixel_size^2 / 1e6, 'FaceColor','b');
else
    bar(0,0);
end
xlabel('Year'); ylabel('Cutoff area (km^2)');
title('Cutoff areas through time','FontSize',14);
% Extract years into a numeric array
all_years = arrayfun(@(x) x.meta.year, riv);
xlim([min(all_years)-1, max(all_years)+1]);
set(gcf, 'Position', get(0,'Screensize') .* [0 0 1 0.4]);


%% ================================
% TABLE 4.3: CHANNEL WIDTH STATISTICS
% ================================

T_width = table;

for i = 1:numel(riv)
    W = riv(i).vec.W;   % already in meters
    
    W_clean = W(~isnan(W) & W > 0);
    
    T_width.Year(i,1) = riv(i).meta.year;
    T_width.MeanWidth_m(i,1) = mean(W_clean);
    T_width.MinWidth_m(i,1)  = min(W_clean);
    T_width.MaxWidth_m(i,1)  = max(W_clean);
    T_width.StdDev_m(i,1)    = std(W_clean);
end

disp('Table 4.3: Channel Width Statistics');
disp(T_width);

writetable(T_width, 'Table_4_3_Channel_Width_Statistics.csv');

%% ================================
% TABLE 4.4: SINUOSITY INDEX
% ================================

T_sinuosity = table;

for i = 1:numel(riv)
    T_sinuosity.Year(i,1) = riv(i).meta.year;
    T_sinuosity.ChannelLength_m(i,1) = riv(i).vec.cl_len;
    
    % recompute valley length (in meters)
    cl = riv(i).vec.cl;
    valley_len = sqrt((cl(end,1)-cl(1,1))^2 + ...
                      (cl(end,2)-cl(1,2))^2) * pixel_size;
    
    T_sinuosity.ValleyLength_m(i,1) = valley_len;
    T_sinuosity.Sinuosity(i,1) = riv(i).vec.sinuosity;
end

disp('Table 4.4: Sinuosity Index');
disp(T_sinuosity);

writetable(T_sinuosity, 'Table_4_4_Sinuosity.csv');

%% ================================
% TABLE 4.5: RIVER NETWORK METRICS
% ================================

T_network = table;

for i = 1:numel(riv)

    I = riv(i).im.hc;

    % Skeletonize
    skel = bwmorph(I, 'skel', Inf);

    % Branch points (nodes)
    nodes = bwmorph(skel, 'branchpoints');
    num_nodes = sum(nodes(:));

    % Endpoints
    endpoints = bwmorph(skel, 'endpoints');
    num_endpoints = sum(endpoints(:));

    % Approximate links
    num_links = num_nodes + num_endpoints;

    % Area (km²)
    area_km2 = sum(I(:)) * pixel_size^2 / 1e6;

    % Density
    node_density = num_nodes / area_km2;

    % Branching intensity
    branching_intensity = num_links / max(num_nodes,1);

    % Store
    T_network.Year(i,1) = riv(i).meta.year;
    T_network.Nodes(i,1) = num_nodes;
    T_network.Links(i,1) = num_links;
    T_network.NodeDensity(i,1) = node_density;
    T_network.BranchingIntensity(i,1) = branching_intensity;

end

disp('Table 4.5: River Network Metrics');
disp(T_network);

writetable(T_network, 'Table_4_5_Network_Metrics.csv');


%% ================================
% TABLE 4.6: CHANNEL MIGRATION
% ================================

T_migration = table('Size',[numel(migration) 4], ...
    'VariableTypes',{'string','double','double','string'}, ...
    'VariableNames',{'YearPair','MeanMigration_m','MaxMigration_m','Direction'});

for i = 1:numel(migration)

    % Mean migration
    mean_mig = migration(i).Mrcl;

    % Max migration (estimate from distance map)
    cl1 = riv(i).vec.cls;
    cl2 = riv(i+1).vec.cls;

   Dmat = pdist2(cl1, cl2);
   minDist = min(Dmat, [], 2);
   minDist_m = minDist * pixel_size;

   max_mig = max(minDist_m);
   mean_mig_cl = mean(minDist_m);

    % Direction (simple interpretation)
    if migration(i).MAa > migration(i).MAe
        direction = "Expansion";
    else
        direction = "Erosion-dominated";
    end

    % Store
    T_migration.YearPair(i,1) = ...
        string(migration(i).year1) + "-" + string(migration(i).year2);

    T_migration.MeanMigration_m(i,1) = mean_mig;
    T_migration.MaxMigration_m(i,1)  = max_mig;
    T_migration.Direction(i,1) = direction;

end

disp('Table 4.6: Channel Migration');
disp(T_migration);

writetable(T_migration, 'Table_4_6_Migration.csv');

%% ================================
% EXTRA TABLE: MORPHODYNAMICS SUMMARY
% ================================

T_extra = table;

for i = 1:numel(riv)

    T_extra.Year(i,1) = riv(i).meta.year;
    T_extra.Width_m(i,1) = riv(i).vec.Wavg;
    T_extra.Sinuosity(i,1) = riv(i).vec.sinuosity;
    T_extra.Curvature(i,1) = riv(i).vec.avg_curvature;
    T_extra.Wavelength_m(i,1) = riv(i).vec.wavelength;
    T_extra.BraidingIndex(i,1) = riv(i).vec.braiding_index;

end

disp('Additional Morphodynamics Table');
disp(T_extra);

writetable(T_extra, 'Morphodynamics_Summary.csv');

%% 10. Time Series of Key Metrics
% Extract time series data
years_all = arrayfun(@(x) x.meta.year, riv);
sinuosity = arrayfun(@(x) x.vec.sinuosity, riv);
braiding = arrayfun(@(x) x.vec.braiding_index, riv);
width = arrayfun(@(x) x.vec.Wavg, riv);      % convert to meters
curvature = arrayfun(@(x) x.vec.avg_curvature, riv);
wavelength = arrayfun(@(x) x.vec.wavelength, riv);        

% Migration rates per interval (already computed)
mig_cl = arrayfun(@(x) x.Mrcl, migration);
mig_erosion = arrayfun(@(x) x.Mre, migration);
mig_accretion = arrayfun(@(x) x.Mra, migration);
years_mig = arrayfun(@(x) x.year2, migration);

net_migration = mig_accretion - mig_erosion;

figure;
plot(years_mig, net_migration, 'k-o','LineWidth',2);
xlabel('Year');
ylabel('Net Migration (m/yr)');
title('Net Channel Migration Balance');
grid on;

% Cumulative erosion and accretion areas (in km²)
cum_erosion = cumsum(arrayfun(@(x) x.MAe, migration)) * pixel_size^2 / 1e6;
cum_accretion = cumsum(arrayfun(@(x) x.MAa, migration)) * pixel_size^2 / 1e6;

% Plot all time series in one figure
figure('Name', 'Time Series of River Metrics', 'Position', [100, 100, 1400, 900]);

subplot(3,2,1);
plot(years_all, sinuosity, '-o', 'LineWidth', 2);
xlabel('Year'); ylabel('Sinuosity'); title('River Sinuosity'); grid on;

subplot(3,2,2);
plot(years_all, width, '-s', 'LineWidth', 2);
xlabel('Year'); ylabel('Average Width (m)'); title('Average Width'); grid on;

subplot(3,2,3);
plot(years_all, curvature, '-d', 'LineWidth', 2);
xlabel('Year'); ylabel('Curvature (pixel^{-1})'); title('Average Curvature'); grid on;

subplot(3,2,4);
plot(years_mig, mig_cl, 'ko-', 'LineWidth', 2); hold on;
plot(years_mig, mig_erosion, 'bs-', 'LineWidth', 2);
legend('Centerline', 'Erosion', 'Location', 'best');
xlabel('Year (end of interval)'); ylabel('Migration Rate (m/yr)'); title('Reach Migration'); grid on;

subplot(3,2,5);
plot(years_all, braiding, '-^', 'LineWidth', 2);
xlabel('Year'); ylabel('Braiding Index'); title('Braiding Index'); grid on;

subplot(3,2,6);
plot(years_all, wavelength, '-p', 'LineWidth', 2);
xlabel('Year'); ylabel('Wavelength (m)'); title('Meander Wavelength'); grid on;

% Additional figure for cumulative erosion/accretion
figure('Name', 'Cumulative Erosion vs Accretion');
plot(years_mig, cum_erosion, 'k-o', 'LineWidth', 2); hold on;
plot(years_mig, cum_accretion, 'r-o', 'LineWidth', 2);
legend('Erosion', 'Accretion', 'Location', 'best');
xlabel('Year (end of interval)'); ylabel('Cumulative Area (km^2)');
title('Erosion vs Accretion Balance'); grid on;

fprintf('Analysis complete.\n');