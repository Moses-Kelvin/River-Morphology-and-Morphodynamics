%% RivMAP Analysis (Multiple Years - Demo Consistent)
close all; clear; clc;

addpath('RivMAP');          % adds the RivMAP folder to the path
pixel_size = 30;   % meters per pixel (adjust to your image resolution)

pwd                 % shows current folder
ls                  % lists files in current folder
%% 1 Load Masks Automatically

 %% 1 Load Masks Automatically
files = dir('mask*.tif');
n = length(files);
if n == 0
    error('No mask files found')
end

Wn = 10;           % nominal channel width (pixels)
es = 'WE';         % adjust if river flows north‑south

for i = 1:n
    fname = files(i).name;
    I = imread(fname) > 0;
    
    % Basic cleaning
    I = imfill(I, 'holes');
    I = bwareaopen(I, 100);
    
    % Braiding index on cleaned (multi‑thread) mask
    CC = bwconncomp(I);
    braiding_index = CC.NumObjects;
    
    % Extract main channel for analysis
    I_main = bwareafilt(I, 1);
    
    % Smoothing dilation
    se = strel('disk', 2);
    I_main = imdilate(I_main, se);
    
    year = str2double(fname(5:8));
    
    riv(i).meta.year = year;
    riv(i).meta.Wn = Wn;
    riv(i).meta.exit_sides = es;
    
    riv(i).im.st = I_main;
    riv(i).im.hc = I_main;% main channel for centerline/banks
    riv(i).vec.braiding_index = braiding_index;
end
%% 2 Plot Masks
 
close all
for i = 1:numel(riv)

    I = riv(i).im.st;

    imshow(I)
    title(riv(i).meta.year)
    pause(0.2)

end

%% 3 Analyze One Year (same as demo)

i = 1;
plotornot = 1;

Wn = riv(i).meta.Wn;
Ist = riv(i).im.st;
Ihc = riv(i).im.hc;
es = riv(i).meta.exit_sides;

%% Centerline

[cl,Icl] = centerline_from_mask(Ist,es,Wn,plotornot);

%% Banklines

banks = banklines_from_mask(Ist,es,plotornot);

%% Width

Wbl = width_from_banklines(cl,banks{1},banks{2},Wn);
Wavg = mean(Wbl(~isnan(Wbl)));

spacing = Wn/2;
[Wm,SWm] = width_from_mask(Ist,cl,spacing);

S = [0; cumsum(sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2))];

% Valley length (straight line)
valley_len = sqrt((cl(end,1)-cl(1,1))^2 + (cl(end,2)-cl(1,2))^2);

% Channel length
channel_len = S(end);

% Sinuosity
sinuosity = channel_len / valley_len;

riv(i).vec.sinuosity = sinuosity;

figure
plot(S,Wbl); hold on
plot(SWm,Wm,'r')
xlabel('streamwise distance')
ylabel('width')

%% Curvature

cls = savfilt(cl,Wn);

A = angles(cls);
C = curvatures(cls);

% Find curvature peaks
[pk,loc] = findpeaks(abs(C));

if length(loc) > 1

    % convert node spacing to distance
    wavelength = mean(diff(S(loc)));

else

    wavelength = NaN;

end

riv(i).vec.wavelength = wavelength;

figure
subplot(2,1,1)
plot(S,A)

subplot(2,1,2)
plot(S,C)

%% Reach Width

Wra = sum(sum(Ist))/S(end);

%% 4 Process All Years

%% 4 Process All Years
plotornot = 0;

for i = 1:numel(riv)
    fprintf('Processing year %d...\n', riv(i).meta.year);

    Wn = riv(i).meta.Wn;
    I = riv(i).im.st;
    es = riv(i).meta.exit_sides;

    % Centerline
    [cl, Icl] = centerline_from_mask(I, es, Wn, plotornot);
    cls = savfilt(cl, Wn);

    % Banklines
    banks = banklines_from_mask(I, es, plotornot);
    lb = banks{1};
    rb = banks{2};
    lbs = savfilt(lb, round(1.5*Wn));
    rbs = savfilt(rb, round(1.5*Wn));

    % Width from banklines
    W = width_from_banklines(cl, lb, rb, Wn);
    Wavg = mean(W(~isnan(W)));

    % Streamwise distance
    S = [0; cumsum(sqrt(diff(cl(:,1)).^2 + diff(cl(:,2)).^2))];

    % Width from mask (optional, for comparison)
    spacing = Wn/2;
    [Wm, SWm] = width_from_mask(I, cl, spacing);

    % Angles and curvature
    A = angles(cls);
    C = curvatures(cls);
    fprintf('Year %d: angles [%f %f], curvature [%f %f]\n', ...
        riv(i).meta.year, min(A), max(A), min(C), max(C));
    fprintf('Year %d: C has %d NaNs, mean abs = %f\n', ...
        riv(i).meta.year, sum(isnan(C)), mean(abs(C(isfinite(C)))));

    % Sinuosity
    valley_len = sqrt((cl(end,1)-cl(1,1))^2 + (cl(end,2)-cl(1,2))^2);
    channel_len = S(end);
    sinuosity = channel_len / valley_len;

    % Wavelength from curvature peaks
    [pk, loc] = findpeaks(abs(C));
    if length(loc) > 1
        wavelength = mean(diff(S(loc)));
    else
        wavelength = NaN;
    end

    % Store results
    riv(i).vec.cl = cl;
    riv(i).vec.cls = cls;
    riv(i).vec.lb = lb;
    riv(i).vec.lbs = lbs;
    riv(i).vec.rb = rb;
    riv(i).vec.rbs = rbs;
    riv(i).vec.W = W;
    riv(i).vec.Wavg = Wavg;
    riv(i).vec.Wra = sum(sum(I)) / S(end);
    riv(i).vec.S = S;
    riv(i).vec.A = A;
    riv(i).vec.C = C;
    riv(i).vec.cl_len = S(end);
    riv(i).vec.sinuosity = sinuosity;
    riv(i).vec.wavelength = wavelength;
    riv(i).im.cl = Icl;

    plotWidthProfile = 1;
    if plotWidthProfile > 0
         % --- Plot width profile for this year ---
        figure('Name', ['Width profile ' num2str(riv(i).meta.year)]);
        plot(S, W, 'b-', 'LineWidth', 1.5); hold on
        plot(SWm, Wm, 'r-', 'LineWidth', 1.5);
        xlabel('Streamwise distance (pixels)');
        ylabel('Width (pixels)');
        legend('Banklines', 'Mask', 'Location', 'best');
        title(['Year ' num2str(riv(i).meta.year) ' - Width profile']);
        grid on;
        saveas(gcf, ['width_profile_' num2str(riv(i).meta.year) '.png']);
    end

end

%% 5 Plot All Years

close all
for i = 1:numel(riv)

    imshow(riv(i).im.st); hold on

    plot(riv(i).vec.cls(:,1),riv(i).vec.cls(:,2),'b','linewidth',1.5)
    plot(riv(i).vec.lb(:,1),riv(i).vec.lb(:,2),'m','linewidth',1.5)
    plot(riv(i).vec.rb(:,1),riv(i).vec.rb(:,2),'m','linewidth',1.5)

    title(riv(i).meta.year)

    pause(0.3)

end

%% 6 Compute Migrations

for i = 1:numel(riv)-1

    for jj = (i+1):numel(riv)
        if isempty(riv(jj).vec) == 0
            i2 = jj;
            break
        end
    end

    cl1 = riv(i).vec.cls;
    cl2 = riv(i2).vec.cls;

    es1 = riv(i).meta.exit_sides;
    es2 = riv(i2).meta.exit_sides;

    sizeI = size(riv(i).im.cl);

    [Imig,Icutsall,cutidcs,cutarea,cutlen,chutelen] = migration_cl(cl1,cl2,es1,es2,Wn,sizeI);

    riv(i).mig.cl.cutyear = riv(i2).meta.year;   % year of the later image

    I1 = riv(i).im.st;
    I2 = riv(i2).im.st;

    [Ie,Ia,Inc,Icutsm] = migration_mask(I1,I2,Wn);

    riv(i).mig.cl.Imig = Imig;
    riv(i).mig.cl.Icuts = Icutsall;

    riv(i).mig.mask.Ie = Ie;
    riv(i).mig.mask.Ia = Ia;

    yrs(i) = riv(i).meta.year;

end

%% 7 Migration Maps

Imig = false(size(riv(1).im.st));
Ie = false(size(riv(1).im.st));
Ia = false(size(riv(1).im.st));

for i = 1:numel(riv)-1

    Imig(riv(i).mig.cl.Imig) = true;
    Ie(riv(i).mig.mask.Ie) = true;
    Ia(riv(i).mig.mask.Ia) = true;

end

figure

subplot(1,3,1)
subimage(Imig)
title('Migrated Area')

subplot(1,3,2)
subimage(Ie)
title('Erosion')

subplot(1,3,3)
subimage(Ia)
title('Accretion')

%% 8 Migration Rates

for i = 1:numel(riv)-1

    MAcl(i) = sum(sum(riv(i).mig.cl.Imig));
    MAe(i) = sum(sum(riv(i).mig.mask.Ie));
    MAa(i) = sum(sum(riv(i).mig.mask.Ia));

    len(i) = riv(i).vec.cl_len;

end

Mrcl = MAcl ./ len;
Mre = MAe ./ len;
Mra = MAa ./ len;
disp('=== Migration diagnostics ===');
disp(['Number of migration intervals: ', num2str(length(yrs))]);
disp(['Migration area (pixels): ', num2str(MAcl)]);
disp(['Erosion area (pixels): ', num2str(MAe)]);
disp(['Centerline length (pixels): ', num2str(len)]);

 
%% 9 Cutoff Analysis (Demo style)
% Combine all cutoff images and collect cutoff areas with their years
Icutoffs = false(size(riv(1).im.st));
cut_years = [];
cut_areas = [];

for i = 1:numel(riv)-1
    if isfield(riv(i).mig.cl, 'Icuts') && ~isempty(riv(i).mig.cl.Icuts)
        Icutoffs = Icutoffs | riv(i).mig.cl.Icuts;
    end
    if isfield(riv(i).mig.cl, 'cutarea') && ~isempty(riv(i).mig.cl.cutarea)
        ncuts = length(riv(i).mig.cl.cutarea);
        cut_years = [cut_years, repmat(riv(i).mig.cl.cutyear, 1, ncuts)];
        cut_areas = [cut_areas, riv(i).mig.cl.cutarea(:)'];
    end
end

% Figure 1: Binary cutoff map
figure;
subplot(1,3,1);
imshow(Icutoffs);
title('Binary map of cutoffs','FontSize',14);

% Figure 2: Cutoffs colored by year
Icutc = zeros(size(riv(1).im.st));   % initialize with zeros
% Assign interval index to cutoff pixels
for i = 1:numel(riv)-1
    if isfield(riv(i).mig.cl, 'Icuts') && ~isempty(riv(i).mig.cl.Icuts)
        Icutc(riv(i).mig.cl.Icuts) = i;
    end
end
subplot(1,3,2);
cmap = parula(numel(riv)-1);   % one color per interval
imshow(Icutc, cmap);
cb = colorbar;
% Set colorbar ticks to show years
ytick_pos = linspace(1, numel(riv)-1, 5);
ytick_labels = round(interp1(1:numel(riv)-1, yrs, ytick_pos));
set(cb, 'Ticks', ytick_pos/(numel(riv)-1), 'TickLabels', ytick_labels);
title('Cutoffs colored by year','FontSize',14);

% Figure 3: Cutoff areas over time (bar chart)
subplot(1,3,3);
if ~isempty(cut_areas)
    bar(cut_years, cut_areas * pixel_size^2 / 1e6, 'FaceColor','b');
else
    bar(0,0);   % empty plot if no cutoffs
end
xlabel('Year'); ylabel('Cutoff area (km^2)');
title('Cutoff areas through time','FontSize',14);
xlim([min(riv(1).meta.year)-1, max(riv(end).meta.year)+1]);

% Adjust figure size
set(gcf, 'Position', get(0,'Screensize') .* [0 0 1 0.4]);

%% 10 Spatial variation of migration rates (Demo style)
% Prepare cell arrays of images for all years
clear Icp Icl Imig Ie Ia
for i = 1:numel(riv)
    Icp{i} = riv(i).im.st;        % channel positions per year
    Icl{i} = riv(i).im.cl;         % centerline images per year
end
for i = 1:numel(riv)-1
    Imig{i} = riv(i).mig.cl.Imig;  % migration areas per interval
    Ie{i}   = riv(i).mig.mask.Ie;  % erosion areas per interval
    Ia{i}   = riv(i).mig.mask.Ia;  % accretion areas per interval
end

% Parameters
Wn = riv(1).meta.Wn;
es = riv(1).meta.exit_sides;
spacing = 2.1 * Wn;   % segment length (2 channel widths)
plotornot = 0;

% Pack the images to analyze into a cell array
Ianalyze{1} = Imig;
Ianalyze{2} = Ie;
Ianalyze{3} = Ia;

% Run spatial_changes
[Iout, Achan, cllen, belt] = spatial_changes(Icp, Icl, Ianalyze, spacing, Wn, es, plotornot);

% Unpack results
Acl = Iout{1};   % migration area per segment per interval
AE  = Iout{2};   % erosion area per segment per interval
AA  = Iout{3};   % accretion area per segment per interval

% Number of years (intervals)
nyrs = numel(riv)-1;

% Create an image of all channel positions through time
Icp_all = false(size(riv(1).im.st));
for i = 1:numel(riv)
    Icp_all = Icp_all | riv(i).im.st;
end

% Create an image of all migrated areas
Imig_all = false(size(riv(1).im.st));
for i = 1:nyrs
    Imig_all = Imig_all | riv(i).mig.cl.Imig;
end

% Plot results (three subplots as in demo 9b)
figure;

% Subplot 1: Meander belt and segments
subplot(1,3,1);
imshow(Icp_all); hold on;
plot(belt.lb(:,1), belt.lb(:,2), 'm', 'LineWidth',1.5);
plot(belt.rb(:,1), belt.rb(:,2), 'm', 'LineWidth',1.5);
for i = 1:length(belt.cl)
    plot([belt.lb(i,1) belt.cl(i,1) belt.rb(i,1)], ...
         [belt.lb(i,2) belt.cl(i,2) belt.rb(i,2)], 'm');
end
title('Meander belt and segments','FontSize',14);

% Subplot 2: Migrated areas with segment numbers
subplot(1,3,2);
imshow(Imig_all); hold on;
plot(belt.cl(:,1), belt.cl(:,2), 'r.', 'MarkerSize',8);
for i = 1:length(belt.cl)
    text(belt.cl(i,1)+5, belt.cl(i,2)+5, num2str(i), 'Color','r', 'FontSize',10);
end
title('Migrated areas with segments','FontSize',14);

% Subplot 3: Spatial migration rate (average over time)
% Average migration rate per segment = (sum of migration areas over intervals) / (cllen * nyrs) * pixel_size
M = sum(Acl, 1) ./ (cllen * nyrs) * pixel_size;   % in m/yr
subplot(1,3,3);
plot(M, 1:length(M), 'k-o', 'LineWidth',1.5);
ylabel('Segment number');
xlabel('Avg. annual migration rate (m/yr)');
title('Spatial variation of migration','FontSize',14);
grid on;

set(gcf, 'Position', get(0,'Screensize') .* [0 0 1 0.4]);

% Section 11
% After Section 10, you have:
% Acl - S x Q matrix (segments x intervals)
% yrs - years for each interval
% cllen - centerline length per segment (needed to convert area to rate)

% Convert migration area to migration rate (m/yr)
% Acl is in pixels, pixel_size = 30 m
migration_rate = (Acl ./ cllen) * pixel_size;   % each column: rates per interval

% Plot spacetime map
figure;
imagesc(yrs, 1:size(migration_rate,1), migration_rate);
xlabel('Year');
ylabel('Segment number');
title('Migration rate (m/yr) per segment');
colorbar;
colormap(jet);


figure
plot(yrs, Mrcl * pixel_size, 'ro-', 'LineWidth', 2); hold on
plot(yrs, Mre * pixel_size, 'bo-', 'LineWidth', 2);

legend('CL','Erosion')
ylabel('Migration Rate (m/yr)')
xlabel('Year')

%% Extract time series metrics
for i = 1:numel(riv)
    years(i) = riv(i).meta.year;
    sinuosity(i) = riv(i).vec.sinuosity;
    braiding(i) = riv(i).vec.braiding_index;
    wavelength(i) = riv(i).vec.wavelength;
    width(i) = riv(i).vec.Wavg * pixel_size;          % convert to meters
    curvature(i) = mean(abs(riv(i).vec.C(isfinite(riv(i).vec.C))));
end

%% Figure: All metrics in one (enlarged)
figure('Position', [100, 100, 1400, 900]);

subplot(3,2,1)
plot(years, sinuosity, '-o', 'LineWidth', 2)
xlabel('Year'); ylabel('Sinuosity'); title('River Sinuosity'); grid on;

subplot(3,2,2)
plot(years, width, '-s', 'LineWidth', 2)
xlabel('Year'); ylabel('Average Width (m)'); title('Average Width'); grid on;

subplot(3,2,3)
plot(years, curvature, '-d', 'LineWidth', 2)
xlabel('Year'); ylabel('Curvature (pixel^{-1})'); title('Average Curvature'); grid on;

subplot(3,2,4)
plot(yrs, Mrcl * pixel_size, 'ko-', 'LineWidth', 2)
xlabel('Year'); ylabel('Migration Rate (m/yr)'); title('Channel Migration'); grid on;

subplot(3,2,5)
plot(years, braiding, '-^', 'LineWidth', 2)
xlabel('Year'); ylabel('Braiding Index'); title('Braiding Index'); grid on;

subplot(3,2,6)
plot(years, wavelength, '-p', 'LineWidth', 2)
xlabel('Year'); ylabel('Wavelength (pixels)'); title('Meander Wavelength'); grid on;

%% Separate figure: Reach average migration rate (CL vs Erosion)
figure('Position', [100, 100, 800, 500]);
if ~isempty(Mrcl) && ~isempty(Mre)
    plot(yrs, Mrcl * pixel_size, 'ro-', 'LineWidth', 2); hold on
    plot(yrs, Mre * pixel_size, 'bo-', 'LineWidth', 2);
    legend('Centerline Migration', 'Erosion', 'Location', 'best');
    xlabel('Year'); ylabel('Migration Rate (m/yr)');
    title('Reach Average Migration Rate'); grid on;
end

%% Separate figure: Cumulative erosion vs accretion
figure('Position', [100, 100, 800, 500]);
if ~isempty(MAe) && ~isempty(MAa)
    plot(yrs, cumsum(MAe) * pixel_size^2 / 1e6, 'k-o', 'LineWidth', 2); hold on
    plot(yrs, cumsum(MAa) * pixel_size^2 / 1e6, 'r-o', 'LineWidth', 2);
    legend('Erosion', 'Accretion', 'Location', 'best');
    xlabel('Year'); ylabel('Cumulative Area (km^2)');
    title('Erosion vs Accretion Balance'); grid on;
end