%% 1. Load GeoTIFF Metadata
mask_file = 'Segment2_activeChannel_mask_2009.tif'; 
info = geotiffinfo(mask_file);
R = info.SpatialRef;

%% 2. Define Pixel (from RivMAP)
target_pixel_x = 1195.0; % column
target_pixel_y = 709.6;  % row

%% 3. Convert Pixel → Map (UTM)
[x_map, y_map] = intrinsicToWorld(R, target_pixel_x, target_pixel_y);


%% 4. Convert Map → Geographic (Lat/Lon)

% Get projection object (THIS is the correct one)
proj = R.ProjectedCRS;

% Convert UTM → Lat/Lon
[lat, lon] = projinv(proj, x_map, y_map);

%% 5. Print results
fprintf('--- Migration Analysis Results ---\n');
fprintf('Pixel: [%.2f, %.2f]\n', target_pixel_x, target_pixel_y);
fprintf('Latitude:  %.6f\n', lat);
fprintf('Longitude: %.6f\n', lon);
fprintf('Search String: %.6f, %.6f\n', lat, lon);