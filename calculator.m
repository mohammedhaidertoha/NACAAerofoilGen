function [C_L, panelstrengths] = calculator(airfoilcode, N, AoA, U_infinity)
% The calculator() function calculates the coefficient of lift, C_L, and
% the strengths of the doublets, panelstrengths, i.e. mu.
% -------------------------------------------------------------------------
% inputs:
% airfoilcode - A four digit code that refers to the NACA airfoil e.g.
% 1410.
% N - The number of panels to be used e.g. 100.
% AoA - The angle of attack in degrees e.g. 10°.
% U_infinity - The freesream velocity in m/s e.g. 15 m/s.
% -------------------------------------------------------------------------
% outputs:
% C_L - The coefficient of lift.
% panelstrengths - The doublet strength.

%% generating and discretising the airfoil using the 'panelgen' function:
[x,z] = panelgen(airfoilcode, N, AoA);

%% calculating the panel angle, beta_i:
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
    B(i, 1) = -U_infinity * sin(AoA - beta_i(i));
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

%% finally, solving for the doublet strength, panelstrengths:
panelstrengths = A\B;

%% calculating the lift coefficient, C_L:
C_L = (-2*panelstrengths(end))/(U_infinity);
end
