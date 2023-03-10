function streamlineplot(airfoilcode, N, AoA, U_infinity)
% The streamlineplot() function produce any 4-digit NACA airfoil with
% velocity field vectors and streamlines.
% -------------------------------------------------------------------------
% inputs:
% airfoilcode - A four digit code that refers to the NACA airfoil e.g.
% 1410.
% N - The number of panels to be used e.g. 100.
% AoA - The angle of attack in degrees e.g. 10°.
% -------------------------------------------------------------------------
% outputs:
% C_L - The coefficient of lift.
% A streamlines plot.
% A velocity field vector plot.

%% generating and discretising the airfoil using the 'panelgen' function:
[x,z] = panelgen(airfoilcode, N, AoA);
AoA_rad = deg2rad(AoA); % convert the AoA from degrees to radians

%% calculating the panel angle, beta_i:
% initialise beta_i with an array full of zeros of dimension 1 by N+1:
beta_i = zeros(1, N+1);
for i = 1:length(x)-1
    beta_i(i) = atan2((z(i + 1) - z(i)), (x(i + 1) - x(i)));
end

%% Solving Ax = B, where A is (N+1)×(N+1) square matrix and B are column vectors of length (N+1) to find the doublet strength, panelstrengths:
% using the provided cdoublet.m script to find u_ij and v_ij of the midpoints of the panels to get A:
% initialise all the variables:
midpoints_coord = zeros(N+1, 2);
start_coord = zeros(N+1, 2);
end_coord = zeros(N+1, 2);
A = zeros(N+1, N+1);
B = zeros(N+1, 1);
u_ij = zeros(N, N+1);
v_ij = zeros(N, N+1);

% a for-loop to find the midpoints, midpoint_coord, of each panel and the start,
% start_coord, and the end points, end_coord:
for i = 1:N+1
    % generating midpoints:
    midpoints_coord(i, 1) = ((x(i)+x(i+1))/2);
    midpoints_coord(i, 2) = ((z(i)+z(i+1))/2);
    % generating start points:
    start_coord(i, 1) = x(i);
    start_coord(i, 2) = z(i);
    % generating end points:
    end_coord(i, 1) = x(i+1);
    end_coord(i, 2) = z(i+1);
end

% a for-loop to compute A and B using cdoublet():
for i = 1:N
    B(i, 1) = -U_infinity * sin(AoA_rad - beta_i(i));
    for j = 1:N+1
        [U_i, V_i] = cdoublet(midpoints_coord(i,:), start_coord(j,:), end_coord(j,:));
        u_ij(i, j) = U_i;
        v_ij(i, j) = V_i;
        A(i, j) = -U_i*sin(beta_i(i)) + V_i*cos(beta_i(i));
    end
end

% fulfilling the Kutta condition of, circulation at the TE is zero, by
% declaring the (N+1)th row:
A(N+1, 1) = 1;
A(N+1, N) = -1;
A(N+1, N+1) = 1;
B(N+1) = 0;

%% finally, solving for the doublet strength, panelstrengths i.e. mu:
panelstrengths = A\B;

%% generating velocity field vectors (arrows) and streamlines around the airfoil:
% producing panels for the plot using meshgrid() function:
resolution = 150;   % the higher the resolution, the finer the mesh, however the greater the compute time
[x_meshgrid, z_meshgrid] = meshgrid(linspace(-0.2, 1.2, resolution), linspace(-0.7, 0.7, resolution));

%% calculating the effect of the aerofoil panels on each grid:
% initialising U and V, each represent the horizontal and vertical
% componenets of the velocities:
U = zeros(resolution, resolution);
V = zeros(resolution, resolution);

for i = 1:resolution
    for j = 1:resolution
        U(i, j) = U_infinity * cos(AoA_rad);
        V(i, j) = U_infinity * sin(AoA_rad);
        for k = 1:N+1
            [U_grid, V_grid] = cdoublet([x_meshgrid(i,j), z_meshgrid(i,j)], [x(k), z(k)], [x(k+1), z(k+1)]);
            U(i, j) = U(i, j) + panelstrengths(k) * U_grid;
            V(i, j) = V(i, j) + panelstrengths(k) * V_grid;
        end
    end
end

%% identifying which U and V points are inside the airfoil and should thus not be plotted:
% checking for points inside the airfoil region:
points_inside = inpolygon(x_meshgrid, z_meshgrid, x, z);

% converting the values found to zero such that they are not plotted:
U(points_inside) = NaN;
V(points_inside) = NaN;

%% plotting the streamlines around the airfoil:
% defining the figure:
streamline_plot = figure;
hold on
% drawing the streamlines with well spacing and direction vectors:
streamslice(x_meshgrid, z_meshgrid, U, V);
axis equal;
% plotting the actual airfoil:
plot(x(1:length(x)-1), z(1:length(z)-1), 'r-', LineWidth=2) 
% specifying the suggested limit:
xlim([-0.2, 1.2])
ylim([-0.7, 0.7])
% labelling the axes and adding an appropriate title:
xlabel('Distance along x-axis')
ylabel('Distance along z-axis')
streamlineplot_title = ['Simulated streamlines around a NACA ', num2str(airfoilcode), ' with N = ', num2str(N), ', at \alpha = ', num2str(rad2deg(AoA)),'° and U_\infty = ', num2str(U_infinity),' m/s.'];
title(streamlineplot_title)
% changing the dimensions of the plot to fit the title:
set(0, 'defaultfigureposition', [1300 10 800 400])
% saving the plot programmatically: 
saveas(streamline_plot,['Streamlines_NACA_', num2str(airfoilcode), '_N_', num2str(N), '_AoA_', num2str(rad2deg(AoA)), '_U_infinity_', num2str(U_infinity), '.png'])
hold off

%% plotting the velocity field vectors around the airfoil:
% decrease the density of quiver field:
step_size = 1:10:length(x_meshgrid);
% defining the figure:
quiver_plot = figure;
hold on
% plotting the velocity field vectors: 
quiver(x_meshgrid(step_size, step_size), z_meshgrid(step_size, step_size), U(step_size, step_size), V(step_size, step_size))
axis equal;
% plotting the actual airfoil:
plot(x(1:length(x)-1), z(1:length(z)-1), 'r-', LineWidth=2)
% specifying the suggested limit:
xlim([-0.2, 1.2])
ylim([-0.7, 0.7])
% labelling the axes and adding an appropriate title:
quiverplot_title = ['Simulated velocity field vectors around a NACA ', num2str(airfoilcode), ' with N = ', num2str(N), ', at \alpha = ', num2str(rad2deg(AoA)),'° and U_\infty = ', num2str(U_infinity),' m/s.'];
title(quiverplot_title)
xlabel('Distance along x-axis')
ylabel('Distance along z-axis')
% changing the dimensions of the plot to fit the title:
set(0, 'defaultfigureposition', [1300 10 800 400])
% saving the plot programmatically:
saveas(quiver_plot, ['Velocity_arrows_NACA_', num2str(airfoilcode), '_N_', num2str(N), '_AoA_', num2str(rad2deg(AoA)), '_U_infinity_', num2str(U_infinity), '.png'])
hold off

%% calculating and displaying the lift coefficient, C_L:
C_L = (-2*panelstrengths(end))/(U_infinity);
fprintf('The coefficient of lift, C_L = %.6f.\n', C_L)
end
