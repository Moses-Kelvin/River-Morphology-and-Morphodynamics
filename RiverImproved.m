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
files = dir('Segment1_activeChannel_mask*.tif');
n = length(files);
if n == 0
    error('No mask files found');
end

% Sort files by year (assuming filename contains year)
years = zeros(n,1);
for i = 1:n
    tokens = regexp(files(i).name, '\d+', 'match');
    if ~isempty(tokens)
        % years(i) = str2double(tokens{1});
        years(i) = str2double(tokens{end}); 
    else
        error('Could not extract year from filename: %s', files(i).name);
    end
end
[~, idx] = sort(years);
files = files(idx);
years = years(idx);

es = 'WE';         % adjust if river flows north‑south (options: 'WE' or 'NS')
plotornot = 0;     % suppress internal RivMAP plots (set to 1 if needed)

% Preallocate structure array
riv = struct('meta', struct(), 'im', struct(), 'vec', struct(), 'mig', struct());

for i = 1:n
    fname = files(i).name;
    I = imread(fname) > 0;   % binary mask
    
    % Basic cleaning: fill holes, remove small isolated islands (adjust threshold)
    I = imfill(I, 'holes');
    I = bwareaopen(I, 20);   % remove objects <100 pixels
    
    % Store original full channel mask
    riv(i).im.hc = I;
    
    % Braiding index (number of separate channel threads)
    CC = bwconncomp(I);
    braiding_index = CC.NumObjects;
    riv(i).vec.braiding_index = braiding_index;
    
    % For centerline extraction, use the largest connected component (main channel)
    % This may lose braiding information but is necessary for a single centerline.
    I_main = bwareafilt(I, 1);   % largest component only
    % Optional morphological smoothing to reduce noise
    % Enforce connectivity
      se = strel('disk', 2);
      I_main = imdilate(I_main, se);
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
    riv(i).meta.Wavg = Wavg;
    riv(i).meta.Wn = Wn_initial;
    
    % --- Streamwise distance ---
    S = [0; cumsum(sqrt(diff(cl(:,1)).^2 + diff(cl(:,2)).^2))];
    riv(i).vec.S = S;
    riv(i).vec.cl = cl;
    riv(i).vec.lb = lb;
    riv(i).vec.rb = rb;
    riv(i).vec.W = W;
    riv(i).vec.Wavg = Wavg;
    
    % --- Width from mask (optional) ---
    try
        [Wm, SWm] = width_from_mask(I_double, cl, round(0.5 * Wavg));
        riv(i).vec.Wm = Wm;
        riv(i).vec.SWm = SWm;
    catch ME
        warning('width_from_mask failed for year %d: %s', riv(i).meta.year, ME.message);
        riv(i).vec.Wm = [];
        riv(i).vec.SWm = [];
    end
    
    % --- Smooth centerline for curvature ---
    cls = savfilt(cl, Wn_initial);
    riv(i).vec.cls = cls;
    
    % --- Curvature and angles ---
    A = angles(cls);
    C = curvatures(cls);
    riv(i).vec.A = A;
    riv(i).vec.C = C;
    
    % --- Sinuosity ---
    valley_len = sqrt((cl(end,1)-cl(1,1))^2 + (cl(end,2)-cl(1,2))^2);
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
    [pk, loc] = findpeaks(abs(C_smooth), 'MinPeakHeight', 0.1*max(abs(C_smooth)));
    if length(loc) > 1
        wavelength = mean(diff(S(loc)));
    else
        wavelength = NaN;
    end
    riv(i).vec.wavelength = wavelength;
    
    % --- Average absolute curvature ---
    riv(i).vec.avg_curvature = mean(abs(C(isfinite(C))));
    
    % --- Reach-averaged width ---
    riv(i).vec.Wra = sum(sum(I)) / S(end);
    
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
        xlabel('Streamwise distance (pixels)');
        ylabel('Width (pixels)');
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
    migration(i).Mrcl = (migration(i).MAcl / migration(i).cl_len) * pixel_size / dt;
    migration(i).Mre  = (migration(i).MAe  / migration(i).cl_len) * pixel_size / dt;
    migration(i).Mra  = (migration(i).MAa  / migration(i).cl_len) * pixel_size / dt;
end

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

%% 9. Spatial Variation of Migration Rates (Segment-Based)
% Prepare data for spatial_changes function
% Icp: cell array of channel masks (full mask) for each year
% Icl: cell array of centerline images (binary) for each year
Icp = cell(1, n);
Icl = cell(1, n);
for i = 1:n
    Icp{i} = riv(i).im.hc;
    % Create binary centerline image from stored centerline coordinates
    Icl_img = false(size(riv(i).im.st));
    cl_x = round(riv(i).vec.cl(:,1));
    cl_y = round(riv(i).vec.cl(:,2));
    valid = cl_x >= 1 & cl_x <= size(Icl_img,2) & cl_y >= 1 & cl_y <= size(Icl_img,1);
    if any(valid)
        ind = sub2ind(size(Icl_img), cl_y(valid), cl_x(valid));
        Icl_img(ind) = true;
    end
    Icl{i} = Icl_img;
end

% Migration areas per interval
Imig_cell = cell(1, length(migration));
Ie_cell = cell(1, length(migration));
Ia_cell = cell(1, length(migration));
for i = 1:length(migration)
    Imig_cell{i} = migration(i).Imig;
    Ie_cell{i} = migration(i).Ie;
    Ia_cell{i} = migration(i).Ia;
end

% Pack the migration images into a cell array for spatial_changes
Ianalyze{1} = Imig_cell;
Ianalyze{2} = Ie_cell;
Ianalyze{3} = Ia_cell;

% Parameters for segmentation
all_Wn = arrayfun(@(x) x.meta.Wn, riv);
Wn = round(median(all_Wn));
% Cap Wn to avoid huge smoothing windows (causes savfilt errors)
if Wn > 20
    Wn = 20;
end
spacing = 2.1 * Wn;   % segment length (2 channel widths)
plotornot = 0;

% Run spatial_changes with error handling
try
    [Iout, Achan, cllen, belt] = spatial_changes(Icp, Icl, Ianalyze, spacing, Wn, es, plotornot);
    
    % Unpack results
    Acl = Iout{1};   % migration area per segment per interval (pixels)
    AE  = Iout{2};   % erosion area per segment per interval
    AA  = Iout{3};   % accretion area per segment per interval
    
    % Number of intervals
    nint = length(migration);
    
    % Create an image of all channel positions through time
    Icp_all = false(size(riv(1).im.st));
    for i = 1:n
        Icp_all = Icp_all | riv(i).im.st;
    end
    
    % Create an image of all migrated areas
    Imig_all = false(size(riv(1).im.st));
    for i = 1:nint
        Imig_all = Imig_all | migration(i).Imig;
    end
    
    % Plot results
    figure('Name', 'Spatial Migration Analysis');
    subplot(1,3,1);
    imshow(Icp_all); hold on;
    plot(belt.lb(:,1), belt.lb(:,2), 'm', 'LineWidth',1.5);
    plot(belt.rb(:,1), belt.rb(:,2), 'm', 'LineWidth',1.5);
    for i = 1:length(belt.cl)
        plot([belt.lb(i,1) belt.cl(i,1) belt.rb(i,1)], ...
             [belt.lb(i,2) belt.cl(i,2) belt.rb(i,2)], 'm');
    end
    title('Meander belt and segments','FontSize',14);
    
    subplot(1,3,2);
    imshow(Imig_all); hold on;
    plot(belt.cl(:,1), belt.cl(:,2), 'r.', 'MarkerSize',8);
    for i = 1:length(belt.cl)
        text(belt.cl(i,1)+5, belt.cl(i,2)+5, num2str(i), 'Color','r', 'FontSize',10);
    end
    title('Migrated areas with segments','FontSize',14);
    
    % Average annual migration rate per segment (m/yr)
    total_years = max(all_years) - min(all_years);
    M_spatial = (sum(Acl, 2) ./ cllen) * pixel_size / total_years;   % m/yr
    
    subplot(1,3,3);
    plot(M_spatial, 1:length(M_spatial), 'k-o', 'LineWidth',1.5);
    ylabel('Segment number');
    xlabel('Avg. annual migration rate (m/yr)');
    title('Spatial variation of migration','FontSize',14);
    grid on;
    set(gcf, 'Position', get(0,'Screensize') .* [0 0 1 0.4]);
    
    % Additional: spacetime plot of migration rates per segment over time
    % Compute per-interval migration rate per segment (m/yr)
    M_interval = zeros(size(Acl));
    for i = 1:nint
        dt = migration(i).year2 - migration(i).year1;
        M_interval(:,i) = (Acl(:,i) ./ cllen) * pixel_size / dt;
    end
    
    figure('Name', 'Spacetime Migration Rates');
    imagesc(years_mig, 1:size(M_interval,1), M_interval);
    xlabel('Year (end of interval)'); ylabel('Segment number');
    title('Migration rate (m/yr) per segment per interval');
    colorbar; colormap(jet);
    
catch ME
    warning('Spatial changes analysis failed: %s', ME.message);
    fprintf('Skipping spatial migration analysis (Section 9).\n');
end

%% 10. Time Series of Key Metrics
% Extract time series data
years_all = arrayfun(@(x) x.meta.year, riv);
sinuosity = arrayfun(@(x) x.vec.sinuosity, riv);
braiding = arrayfun(@(x) x.vec.braiding_index, riv);
width = arrayfun(@(x) x.vec.Wavg, riv) * pixel_size;      % convert to meters
curvature = arrayfun(@(x) x.vec.avg_curvature, riv);
wavelength = arrayfun(@(x) x.vec.wavelength, riv);        % in pixels (convert to meters if needed)

% Migration rates per interval (already computed)
mig_cl = arrayfun(@(x) x.Mrcl, migration);
mig_erosion = arrayfun(@(x) x.Mre, migration);
mig_accretion = arrayfun(@(x) x.Mra, migration);
years_mig = arrayfun(@(x) x.year2, migration);

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
xlabel('Year'); ylabel('Wavelength (pixels)'); title('Meander Wavelength'); grid on;

% Additional figure for cumulative erosion/accretion
figure('Name', 'Cumulative Erosion vs Accretion');
plot(years_mig, cum_erosion, 'k-o', 'LineWidth', 2); hold on;
plot(years_mig, cum_accretion, 'r-o', 'LineWidth', 2);
legend('Erosion', 'Accretion', 'Location', 'best');
xlabel('Year (end of interval)'); ylabel('Cumulative Area (km^2)');
title('Erosion vs Accretion Balance'); grid on;

fprintf('Analysis complete.\n');