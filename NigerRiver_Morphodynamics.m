%% RivMAP Analysis (Multiple Years)
% This script analyzes multi-year river channel masks using the RivMAP toolbox.
% It computes centerlines, widths, sinuosity, curvature, cutoffs,
% migration rates, and spatial-temporal patterns.
%
% Requirements:
%   - RivMAP toolbox in MATLAB path.
%   - Mask images in ".tif" format.
%   - Pixel size (meters per pixel) must be set correctly.
%   - Signal Processing Toolbox (for findpeaks).
%
% Author: [Kelvin Akhere Moses]

close all; clear; clc;

% Add RivMAP path (adjust if necessary)
addpath('RivMAP');          % adds the RivMAP folder to the path

%% User Settings
pixel_size = 30;            % meters per pixel (adjust to your image resolution)
% Options
plot_individual = true;     % set to true to display each year's masks
plot_width_profile = true;  % set to true to save width profile figures
plot_migration_maps = true; % set to true to show migration maps

%% 1. Load Masks Automatically
files = dir('Segment2_activeChannel_mask*.tif');
n = length(files);
if n == 0
    error('No mask files found');
end

% Sort files by year (assuming filename contains year)
years = zeros(n,1);
for i = 1:n
    tokens = regexp(files(i).name, '\d+', 'match');
    if ~isempty(tokens)
        years(i) = str2double(tokens{end});    
        % years(i) = str2double(tokens{1});
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

    I = imfill(I, 'holes');
    % Basic cleaning: remove small isolated islands, morphological smoothing
    I = bwareaopen(I, 20);           % remove objects <20 pixels
    I = imclose(I, strel('disk',3));
    I = imopen(I, strel('disk',2));
    
    % Store original full channel mask (with islands)
    riv(i).im.hc = I;
    
    % Braiding index (number of separate channel threads)
    CC = bwconncomp(I);
    stats = regionprops(CC, 'Area');
    areas = [stats.Area];
    threshold = 500;   % fixed pixel threshold for significant channels
    braiding_index = sum(areas > threshold);
    riv(i).vec.braiding_index = braiding_index;
    
    % For centerline extraction, use the largest connected component
    I_main = bwareafilt(I, 1);   % largest component only
    % I_main = I;   % keep full water mask
    % I_main = bwareaopen(I_main, 50);   % remove tiny noise only
    
    % Optional morphological smoothing (do NOT fill holes to preserve islands)
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

%% 3. Per-Year Analysis (Centerline, Width, Curvature, etc.)
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

    riv(i).im.cl = Icl;          % store the binary centerline image
    
    % --- Width from mask (primary method) ---
    spacing = round(Wn_initial);   % sample every ~Wn_initial pixels
    [W, SW] = width_from_mask(I_double, cl, spacing);

    % Force column vectors and consistent size
    W = W(:);
    SW = SW(:);

  % Remove invalid values BEFORE storing
    valid = ~isnan(W) & W > 0 & ~isnan(SW);
    W = W(valid);
    SW = SW(valid);

    W = W * pixel_size;            % convert to meters
    SW = SW * pixel_size;          % streamwise distance (m)
    
    % Clean NaN/zero values
    W_clean = W(~isnan(W) & W > 0);
    Wavg = mean(W_clean);
    
   
    riv(i).vec.W = W;
    riv(i).vec.SW = SW;            % streamwise coordinate for width profile
    riv(i).vec.Wavg = Wavg;
    riv(i).meta.Wavg = Wavg;
    riv(i).meta.Wn = Wn_initial;
    
    % --- Streamwise distance along centerline ---
    S = [0; cumsum(sqrt(diff(cl(:,1)).^2 + diff(cl(:,2)).^2))] * pixel_size;
    riv(i).vec.S = S;
    riv(i).vec.cl = cl;
    
    % No banklines stored
    riv(i).vec.lb = [];
    riv(i).vec.rb = [];
    
    % --- Smooth centerline for curvature ---
    cls = savfilt(cl, round(Wn_initial));
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
    
    % --- Wavelength from curvature peaks ---
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
    
    % --- Plot width profile ---
    if plot_width_profile
        figure('Name', sprintf('Width Profile %d', riv(i).meta.year));
        plot(SW, W, 'b-', 'LineWidth', 1.5);
        xlabel('Streamwise distance (m)');
        ylabel('Width (m)');
        title(sprintf('Year %d - Width Profile (mask)', riv(i).meta.year));
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
        % No banklines to plot
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
    
    % Use the later year's Wn for migration functions  
    % Wn = riv(i+1).meta.Wn;
    Wn = round(riv(i+1).meta.Wn);
    % sizeI = size(riv(i).im.st);
    % sizeI = round(size(riv(i).im.st));
    sizeI = size(riv(i).im.st);
     
    
    % Centerline-based migration and cutoffs
    [Imig, Icutsall, cutidcs, cutarea, cutlen, chutelen] = ...
        migration_cl(cl1, cl2, es1, es2, Wn, sizeI);
    
    % Mask-based erosion/accretion
    I1 = imfill(riv(i).im.hc, 'holes');
    I2 = imfill(riv(i+1).im.hc, 'holes');
    % I1 = riv(i).im.hc;   % full channel mask 
    % I2 = riv(i+1).im.hc;
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

%% Export total maps as GeoTIFF for GIS Software ( WITH CRS)

fprintf('Exporting total migration maps as GeoTIFF (UTM Zone 31N)...\n');

% Convert logical arrays to uint8
mig_uint8  = uint8(Imig_total);
eros_uint8 = uint8(Ie_total);
acc_uint8  = uint8(Ia_total);

try
    % Read spatial referencing from original mask
    [~, R] = readgeoraster(files(1).name);
    
    % Define CRS explicitly: WGS84 / UTM Zone 31N
    EPSG_CODE = 32631;

    % Write GeoTIFFs with CRS
    geotiffwrite('Total_Migration_Area.tif', mig_uint8, R, ...
        'CoordRefSysCode', EPSG_CODE);

    geotiffwrite('Total_Erosion.tif', eros_uint8, R, ...
        'CoordRefSysCode', EPSG_CODE);

    geotiffwrite('Total_Accretion.tif', acc_uint8, R, ...
        'CoordRefSysCode', EPSG_CODE);

    mig_intensity = double(Imig_total) * pixel_size;

    geotiffwrite('Migration_Intensity_m.tif', mig_intensity, R, ...
        'CoordRefSysCode', EPSG_CODE);

    fprintf('GeoTIFFs successfully exported with CRS (EPSG:32631).\n');

catch ME
    warning('GeoTIFF export failed: %s', ME.message);
    
    % Fallback (not recommended for GIS use)
    imwrite(mig_uint8,  'Total_Migration_Area.tif');
    imwrite(eros_uint8, 'Total_Erosion.tif');
    imwrite(acc_uint8,  'Total_Accretion.tif');
    
    fprintf('Fallback: Exported as plain TIFF (no georeferencing).\n');
end

if plot_migration_maps
    figure('Name', 'Accumulated Migration Maps');
    subplot(1,3,1); imshow(Imig_total); title('Total Migration Area');
    subplot(1,3,2); imshow(Ie_total); title('Total Erosion');
    subplot(1,3,3); imshow(Ia_total); title('Total Accretion');
end

%% 7. Reach-Averaged Migration Rates (per Interval)
% Compute migration rates (meters per year)
for i = 1:length(migration)
    dt = migration(i).year2 - migration(i).year1;
    if dt <= 0
        error('Non‑positive time interval at index %d', i);
    end
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


%% 9. Spatial Variation of Migration Rates (Meander Belt)
fprintf('Computing spatial migration patterns...\n');

% --- Collect necessary images for spatial_changes ---
Icp = cell(1, n);          % channel positions (single‑thread masks)
Icl = cell(1, n);          % centerline images
for i = 1:n
    Icp{i} = riv(i).im.st;          % single‑thread mask
    Icl{i} = riv(i).im.cl;          % centerline image (binary)
end

% Collect migration images over each interval
Imig = cell(1, n-1);
Ie   = cell(1, n-1);
Ia   = cell(1, n-1);
for i = 1:n-1
    Imig{i} = migration(i).Imig;
    Ie{i}   = migration(i).Ie;
    Ia{i}   = migration(i).Ia;
end

% Parameters for spatial analysis
Wn = riv(1).meta.Wn;               % nominal width (pixels)
es = riv(1).meta.exit_sides;       % exit sides (e.g., 'NS')
spacing = 2.1 * Wn;                % sampling interval along belt (pixels)
plotornot_spatial = 0;             % suppress intermediate plots

% Pack the three types of change into a cell for spatial_changes
Ianalyze{1} = Imig;
Ianalyze{2} = Ie;
Ianalyze{3} = Ia;

% Run spatial changes (this may take a moment)
[Iout, Achan, cllen, belt] = spatial_changes(Icp, Icl, Ianalyze, spacing, Wn, es, plotornot_spatial);

% Unpack results:
% Acl, AE, AA are matrices: rows = belt segments, columns = time intervals
Acl = Iout{1};   % centerline migration area per belt polygon per interval (pixels²)
AE  = Iout{2};   % erosion area per belt polygon per interval (pixels²)
AA  = Iout{3};   % accretion area per belt polygon per interval (pixels²)

% Sum across intervals to get total per belt segment
Acl_total = sum(Acl, 2);
AE_total  = sum(AE, 2);
AA_total  = sum(AA, 2);

% Compute average migration width (pixels) per segment
mig_width_pix = Acl_total ./ cllen;          % element‑wise division
% Convert to meters
mig_width_m   = mig_width_pix * pixel_size;
% Total elapsed years
total_years   = years(end) - years(1);
% Average annual migration rate (m/yr) per segment
mig_rate_m_per_yr = mig_width_m / total_years;

% Similarly for erosion and accretion
eros_width_pix = AE_total ./ cllen;
acc_width_pix  = AA_total ./ cllen;
eros_rate_m_per_yr = (eros_width_pix * pixel_size) / total_years;
acc_rate_m_per_yr  = (acc_width_pix * pixel_size) / total_years;

% --- Visualizations ---
% close all;

% Combine all channel masks
Icp_all = false(size(riv(1).im.st));
for i = 1:n
    Icp_all = Icp_all | Icp{i};
end

% Combine all migration areas
Imig_all = false(size(riv(1).im.st));
for i = 1:n-1
    Imig_all = Imig_all | Imig{i};
end

% 1. Meander belt and its centerline
 figure('Name', 'Meander Belt');
imshow(Icp_all); hold on;

plot(belt.lb(:,1), belt.lb(:,2), 'm', 'LineWidth', 1.5);
plot(belt.rb(:,1), belt.rb(:,2), 'm', 'LineWidth', 1.5);

for i = 1:length(belt.cl)
    plot([belt.lb(i,1) belt.cl(i,1) belt.rb(i,1)], ...
         [belt.lb(i,2) belt.cl(i,2) belt.rb(i,2)], 'm');
end

title('Meander Belt and Segments');
set(gca, 'FontSize', 12);
axis tight;

% Belt Centreline
figure('Name', 'Belt Centerline Nodes');

imshow(Imig_all); hold on;
plot(belt.cl(:,1), belt.cl(:,2), 'r.', 'MarkerSize', 10);

for i = 1:length(belt.cl)
    text(belt.cl(i,1)+5, belt.cl(i,2)+5, num2str(i), ...
        'Color', 'r', 'FontSize', 8);
end

title('Belt Centerline Nodes');
set(gca, 'FontSize', 12);
axis tight;

% spatial migration rate 
figure('Name', 'Spatial Migration Rate');

plot(mig_rate_m_per_yr, 1:length(mig_rate_m_per_yr), ...
    'k', 'LineWidth', 1.5);

ylim([0 length(mig_rate_m_per_yr)]);
ylabel('Along-stream distance');
xlabel('Migration rate (m/yr)');
title('Spatial Variation of Migration Rate');
grid on;
set(gca, 'FontSize', 12);
axis tight;
fprintf('Spatial analysis completed.\n');

%% Max Migration points
n_mig = length(migration);

max_points = zeros(n_mig, 6);   % [ID, X, Y, Lat, Lon, MaxMig]
year_pairs = strings(n_mig, 1); % string array (cleaner than cell)

for i = 1:n_mig

    cl1 = riv(i).vec.cls;
    cl2 = riv(i+1).vec.cls;

    [distances, ~] = bwdist(riv(i).im.cl);
    migration_values = distances(riv(i+1).im.cl) * pixel_size;

    [max_val, max_idx] = max(migration_values);
    max_coord = cl2(max_idx, :);

    [x_map, y_map] = intrinsicToWorld(R, max_coord(1), max_coord(2));
    [lat, lon] = projinv(R.ProjectedCRS, x_map, y_map);

    % Store directly (NO concatenation)
    max_points(i,:) = [ ...
        i, ...
        max_coord(1), ...
        max_coord(2), ...
        lat, ...
        lon, ...
        max_val];

    year_pairs(i) = string(migration(i).year1) + "-" + string(migration(i).year2);
end
T_max = array2table(max_points, ...
    'VariableNames', {'ID','Pixel_X','Pixel_Y','Latitude','Longitude','MaxMigration_m'});

T_max.YearPair = year_pairs;

writetable(T_max, 'Max_Migration_Points.csv');

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
    skel = bwmorph(I, 'skel', Inf);
    nodes = bwmorph(skel, 'branchpoints');
    num_nodes = sum(nodes(:));
    endpoints = bwmorph(skel, 'endpoints');
    num_endpoints = sum(endpoints(:));
    num_links = num_nodes + num_endpoints;
    area_km2 = sum(I(:)) * pixel_size^2 / 1e6;
    node_density = num_nodes / area_km2;
    branching_intensity = num_links / max(num_nodes,1);
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
    mean_mig = migration(i).Mrcl;
    cl1 = riv(i).vec.cls;
    cl2 = riv(i+1).vec.cls;
    Dmat = pdist2(cl1, cl2);
    minDist = min(Dmat, [], 2);
    minDist_m = minDist * pixel_size;
    max_mig = max(minDist_m);
    if migration(i).MAa > migration(i).MAe
        direction = "Expansion";
    else
        direction = "Erosion-dominated";
    end
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
%% 10. Time Series of Key Metrics (with corrected axes)
years_all = arrayfun(@(x) x.meta.year, riv);
sinuosity = arrayfun(@(x) x.vec.sinuosity, riv);
braiding = arrayfun(@(x) x.vec.braiding_index, riv);
width = arrayfun(@(x) x.vec.Wavg, riv);
curvature = arrayfun(@(x) x.vec.avg_curvature, riv);
wavelength = arrayfun(@(x) x.vec.wavelength, riv);        

mig_cl = arrayfun(@(x) x.Mrcl, migration);
mig_erosion = arrayfun(@(x) x.Mre, migration);
mig_accretion = arrayfun(@(x) x.Mra, migration);
years_mig = arrayfun(@(x) x.year2, migration);
net_migration = mig_accretion - mig_erosion;

% Net migration balance figure
figure;
plot(years_mig, net_migration, 'k-o','LineWidth',2);
xlabel('Year'); ylabel('Net Migration (m/yr)');
title('Net Channel Migration Balance');
grid on;
xlim([min(years_mig), max(years_mig)]);
xticks(min(years_mig):2:max(years_mig));

% Cumulative erosion/accretion
cum_erosion = cumsum(arrayfun(@(x) x.MAe, migration)) * pixel_size^2 / 1e6;
cum_accretion = cumsum(arrayfun(@(x) x.MAa, migration)) * pixel_size^2 / 1e6;

% Main time series figure
figure('Name', 'Time Series of River Metrics', 'Position', [100, 100, 1400, 900]);

subplot(3,2,1);
plot(years_all, sinuosity, '-o', 'LineWidth', 2);
xlabel('Year'); ylabel('Sinuosity'); title('River Sinuosity'); grid on;
xlim([min(years_all), max(years_all)]); xticks(min(years_all):2:max(years_all));

subplot(3,2,2);
plot(years_all, width, '-s', 'LineWidth', 2);
xlabel('Year'); ylabel('Average Width (m)'); title('Average Width'); grid on;
xlim([min(years_all), max(years_all)]); xticks(min(years_all):2:max(years_all));

subplot(3,2,3);
plot(years_all, curvature, '-d', 'LineWidth', 2);
xlabel('Year'); ylabel('Curvature (pixel^{-1})'); title('Average Curvature'); grid on;
xlim([min(years_all), max(years_all)]); xticks(min(years_all):2:max(years_all));

subplot(3,2,4);
plot(years_mig, mig_cl, 'ko-', 'LineWidth', 2); hold on;
plot(years_mig, mig_erosion, 'bs-', 'LineWidth', 2);
legend('Centerline', 'Erosion', 'Location', 'best');
xlabel('Year (end of interval)'); ylabel('Migration Rate (m/yr)'); title('Reach Migration'); grid on;
xlim([min(years_mig), max(years_mig)]); xticks(min(years_mig):2:max(years_mig));

subplot(3,2,5);
plot(years_all, braiding, '-^', 'LineWidth', 2);
xlabel('Year'); ylabel('Braiding Index'); title('Braiding Index'); grid on;
xlim([min(years_all), max(years_all)]); xticks(min(years_all):2:max(years_all));

subplot(3,2,6);
plot(years_all, wavelength, '-p', 'LineWidth', 2);
xlabel('Year'); ylabel('Wavelength (m)'); title('Meander Wavelength'); grid on;
xlim([min(years_all), max(years_all)]); xticks(min(years_all):2:max(years_all));

% Cumulative erosion/accretion figure
figure('Name', 'Cumulative Erosion vs Accretion');
plot(years_mig, cum_erosion, 'k-o', 'LineWidth', 2); hold on;
plot(years_mig, cum_accretion, 'r-o', 'LineWidth', 2);
legend('Erosion', 'Accretion', 'Location', 'best');
xlabel('Year (end of interval)'); ylabel('Cumulative Area (km^2)');
title('Erosion vs Accretion Balance'); grid on;
xlim([min(years_mig), max(years_mig)]); xticks(min(years_mig):2:max(years_mig));

fprintf('Analysis complete.\n');

% Find Max Migration Point for 2009-2011 (assuming i=5 is 2009 and i+1=6 is 2011)
idx1 = find(years == 2009);
idx2 = find(years == 2011);

cl2009 = riv(idx1).vec.cls;
cl2011 = riv(idx2).vec.cls;

% Simple distance check between centerlines
[distances, ~] = bwdist(riv(idx1).im.cl); % Distance from 2009 centerline
migration_values = distances(riv(idx2).im.cl) * pixel_size; 

[max_val, max_idx] = max(migration_values);
max_coord_pixel = cl2011(max_idx, :);

fprintf('Max Migration of %.2f m found at Pixel (X: %.1f, Y: %.1f)\n', max_val, max_coord_pixel(1), max_coord_pixel(2));
% Note: You will need to use 'pix2latlon' or your georeferencing info to convert these to Dec Degrees.
