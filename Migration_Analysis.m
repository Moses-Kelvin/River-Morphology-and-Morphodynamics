%% ==========================================================
% FULL STANDALONE: RIVMAP + MIGRATION + VECTORS + EXPORT
% ==========================================================

close all; clear; clc;
addpath('RivMAP');

%% SETTINGS
pixel_size = 30;
EPSG_CODE = 32631;
es = 'WE';

%% LOAD MASKS
files = dir('Segment2_activeChannel_mask*.tif');
n = length(files);

if n < 2
    error('Need at least 2 mask files');
end

% Extract years
years = zeros(n,1);
for i = 1:n
    tokens = regexp(files(i).name, '\d+', 'match');
    years(i) = str2double(tokens{end});
end

[years, idx] = sort(years);
files = files(idx);

fprintf('Loaded years: %s\n', num2str(years'));

%% GEOREFERENCE
[~, R] = readgeoraster(files(1).name);
img_size = size(imread(files(1).name));

%% STORAGE
riv = struct();

%% ==========================================================
% 1. PROCESS EACH YEAR (LIKE YOUR ORIGINAL SCRIPT)
%% ==========================================================
for i = 1:n

    I = imread(files(i).name);
    I = I > 0;
    I = logical(I);

    I = imfill(I,'holes');
    I = bwareaopen(I,20);
    I = imclose(I, strel('disk',3));

    I_main = bwareafilt(I,1);

    % Width estimation
    D = bwdist(~I_main);
    width_map = 2 * D;
    Wmed = median(width_map(I_main));

    Wn = round(0.4 * Wmed);
    Wn = max(round(0.25 * Wmed), min(Wn, 80));
    Wn = double(Wn);   % 🔥 critical fix

    % Centerline
    [cl, Icl] = centerline_from_mask(double(I_main), es, Wn, 0);

    if isempty(cl)
        error('Centerline failed for year %d', years(i));
    end

    riv(i).vec.cls = cl;
    riv(i).im.cl = Icl;
    riv(i).meta.year = years(i);

end

fprintf('Centrelines extracted successfully.\n');

%% ==========================================================
% 2. CENTRELINE MIGRATION MAP
%% ==========================================================

centreline_stack = zeros(img_size);

for i = 1:n
    cl = riv(i).vec.cls;

    x = round(cl(:,1));
    y = round(cl(:,2));

    valid = x > 0 & y > 0 & ...
            x <= img_size(2) & ...
            y <= img_size(1);

    idx = sub2ind(img_size, y(valid), x(valid));
    centreline_stack(idx) = i;
end

geotiffwrite('Centreline_Migration_Map.tif', ...
    uint16(centreline_stack), R, ...
    'CoordRefSysCode', EPSG_CODE);

fprintf('✔ Centreline migration map exported.\n');

%% ==========================================================
% 3. MIGRATION VECTORS
%% ==========================================================

X=[]; Y=[]; dX=[]; dY=[]; Mag=[]; YearPair=[];

for i = 1:n-1

    cl1 = riv(i).vec.cls;
    cl2 = riv(i+1).vec.cls;

    Dmat = pdist2(cl1, cl2);
    [~, idx] = min(Dmat, [], 2);

    matched = cl2(idx,:);

    dx = (matched(:,1) - cl1(:,1)) * pixel_size;
    dy = (matched(:,2) - cl1(:,2)) * pixel_size;

    mag = sqrt(dx.^2 + dy.^2);

    [x_map, y_map] = intrinsicToWorld(R, cl1(:,1), cl1(:,2));

    X = [X; x_map];
    Y = [Y; y_map];
    dX = [dX; dx];
    dY = [dY; dy];
    Mag = [Mag; mag];

    YearPair = [YearPair; ...
        repmat(string(years(i)) + "-" + string(years(i+1)), length(mag),1)];

end

T_vec = table(X,Y,dX,dY,Mag,YearPair,...
    'VariableNames',{'X','Y','dX_m','dY_m','Magnitude_m','YearPair'});

writetable(T_vec,'Migration_Vectors.csv');

fprintf('✔ Migration vectors exported.\n');

%% ==========================================================
% 4. VECTOR MAGNITUDE RASTER
%% ==========================================================

vector_raster = zeros(img_size);

for i = 1:length(X)

    [col,row] = worldToIntrinsic(R, X(i), Y(i));
    col = round(col); row = round(row);

    if row>0 && col>0 && row<=img_size(1) && col<=img_size(2)
        vector_raster(row,col) = Mag(i);
    end
end

geotiffwrite('Migration_Vector_Magnitude.tif', ...
    single(vector_raster), R, ...
    'CoordRefSysCode', EPSG_CODE);

fprintf('✔ Vector magnitude raster exported.\n');

%% ==========================================================
% 5. MAX MIGRATION POINTS
%% ==========================================================

max_map = zeros(img_size);
Xmax=[]; Ymax=[]; MaxVal=[]; YearPair_max=[];

for i = 1:n-1

    [distances, ~] = bwdist(riv(i).im.cl);
    migration_values = distances(riv(i+1).im.cl) * pixel_size;

    [max_val, idx] = max(migration_values);
    coord = riv(i+1).vec.cls(idx,:);

    x = round(coord(1));
    y = round(coord(2));

    if x>0 && y>0 && x<=img_size(2) && y<=img_size(1)
        max_map(y,x) = max_val;
    end

    [x_map, y_map] = intrinsicToWorld(R, coord(1), coord(2));

    Xmax = [Xmax; x_map];
    Ymax = [Ymax; y_map];
    MaxVal = [MaxVal; max_val];

    YearPair_max = [YearPair_max; ...
        string(years(i)) + "-" + string(years(i+1))];
end

geotiffwrite('Max_Migration_Points_Map.tif', ...
    single(max_map), R, ...
    'CoordRefSysCode', EPSG_CODE);

T_max = table(Xmax,Ymax,MaxVal,YearPair_max,...
    'VariableNames',{'X','Y','MaxMigration_m','YearPair'});

writetable(T_max,'Max_Migration_Points_Labelled.csv');

fprintf('✔ Max migration outputs exported.\n');

%% ==========================================================
% 6. VISUAL CHECK
%% ==========================================================

figure; imshow(centreline_stack, []); title('Centreline Migration');
figure; imshow(vector_raster, []); title('Vector Magnitude');
figure; imshow(max_map, []); title('Max Migration');

fprintf('ALL DONE SUCCESSFULLY.\n');gain 