close all

flow = []; % Initialize array for X coordinates
loading = []; % Initialize array for Y coordinates
efficiency = []; % Initialize array for percentage values

percentages = [0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87, 0.86, 0.85];
for i = 1:10
    currentData = readmatrix("Efficiency_Data.xlsx", "Sheet", i);
    flow = [flow; currentData(:, 1)]; % append X coordinates
    loading = [loading; currentData(:, 2)]; % append Y coordinates
    efficiency = [efficiency; repmat(percentages(i), size(currentData, 1), 1)]; % append percentage values
end

% Define the grid for interpolation
step = 0.01;
[flow_query, loading_query] = meshgrid(min(flow):step:max(loading), min(loading):step:max(loading)); % step is the grid resolution

% Interpolate data
efficiency_interpolated = griddata(flow, loading, efficiency, flow_query, loading_query, 'cubic');

% % Apply Gaussian smoothing
% sigma = 1;
% efficiency_smooth = imgaussfilt(efficiency_interpolated, sigma); % sigma is the standard deviation of the Gaussian filter

% Vq now contains the interpolated 2D array over the XY plane
contour(flow_query, loading_query, efficiency_smooth)
shading interp % This makes the plot smoother
colorbar % Adds a color bar to indicate values
zlim([0,1])