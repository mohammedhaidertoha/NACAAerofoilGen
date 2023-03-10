% This script plots velocity field vectors (arrows) and streamlines around
% any airfoil using panelgen() function and Panel Code and saves the plots 
% programmatically. However, if the airfoil is NACA 2412, the script
% calculates C_L for range of AoA and N, reads the XFOIL data and plots C_L
% vs AoA.
% -------------------------------------------------------------------------
% inputs:
% airfoilcode - A four digit code that refers to the NACA airfoil e.g.
% 1410.
% U_infinity - The freesream velocity in m/s e.g. 15 m/s.
% AoA - The angle of attack in degrees e.g. 10°.
% N - The number of panels to be used e.g. 100.
% -------------------------------------------------------------------------
% outputs:
% C_L - The airfoil's lift coeffic2ient.
% Plots with the velocity field vectors and streamlines.
% C_L vs AoA plot for the special case when the airfoil is a NACA 2412.

%% housekeeping
clear
clc
close all
format default

%% prompting for user input for the airfoilcode in string format:
airfoilcode = input('Enter the 4-digit NACA airfoil code (e.g. 2415): ', 's');

%% raising errors to the user if the 'airfoilcode' input is invalid:
if size(airfoilcode) ~= 4    % if the input is not a 4-digit integer 
    error('The NACA airfoil code you have entered is invalid. Please enter a valid 4-digit NACA airfoil code.')
end

%% an if-statement that carries on the promp for inputs from the user, which generates and plots the streamlines around the inputted airfoil:
if str2double(airfoilcode) ~= 2412
    % input for Angle of Attack, AoA:
    AoA = input('Enter the angle of attack in degrees (e.g. 10°): ');
    AoA = deg2rad(AoA); % convert from degrees to radians because calculator() accepts radians as input

    % input for number of panels, N:
    N = input('Enter the number of panels to be used (e.g. 100): ');
    % raising errors if N is not a positive integer
    if N <= 0
        error('Please enter a valid N. The number of panels must be positive.')
    elseif mod(N,1) ~= 0
        error('Please enter a valid N. The number of panels must be integers.')
    elseif mod(N,2) ~= 0
        error('Please enter a valid N. The number of panels must be even.')
    end

    % input for freestream velocity, U_infinity:
    U_infinity = input('Enter the freestream velocity in m/s (e.g. 15 m/s): ');
    if U_infinity <= 0
        error('Please enter a valid freestream velocity. The value must be positive.')
    end

    %% generating and discretising the airfoil using the panelgen() and streamlineplot() functions:
    [x, z] = panelgen(airfoilcode, N, AoA);
    [C_L, panelstrengths] = calculator(airfoilcode,N, AoA, U_infinity);
    streamlineplot(airfoilcode, N, AoA, U_infinity);
    
%% The special case when the airfoil is the NACA 2412:
else
    %% read and extract data from the XFOIL file:
    XFOIL = readmatrix("xf-naca2412-il-1000000.txt");      % read the provided XFOIL file
    AoA_XFOIL = XFOIL(70:108, 1);                          % locate and extract the C_L values
    C_L_XFOIL = XFOIL(70:108, 2);                          % locate and extract the AoA values

    %% define necessary parameters such as N, AoA and U_
    N = [50 100 200];                           % the provided number of panels
    AoA = 0:10;                                 % the provided range of angles
    AoA_rad = deg2rad(AoA);                     % convert from degrees to radians because calculator() accepts radians as input
    U_infinity = 15;                            % the provided freestream velocity
    
    %% iterate through the values of AoA corresponding to N values thus the C_L values:
    for i = 1:length(N)
        for j = 1:length(AoA)
            [C_L, panelstrengths] = calculator(airfoilcode, N(i), AoA_rad(j), U_infinity);
            C_L_2412(i, j) = C_L;
        end
    end
    
    %% plotting the XFOIL data alongside the C_L vs AoA for each respective N:
    figure
    plot_XFOIL = plot(AoA_XFOIL, C_L_XFOIL, 'k-', LineWidth=1.5);   % plot the extracted XFOIL data
    hold on
    plot_N50 = plot(AoA, C_L_2412(1,:), '--', LineWidth=1.5);       % plot the curve for N = 50
    plot_N100 = plot(AoA, C_L_2412(2,:), '--x', LineWidth=1.5);     % plot the curve for N = 100
    plot_N200 = plot(AoA, C_L_2412(3,:), '-.', LineWidth=1.5);      % plot the curve for N = 200

    %% formatting and labelling the plots appropriately:
    figure_title = 'A plot of C_L versus \alpha for different cases of N against the XFOIL data';
    title(figure_title)
    figure_xlabel = 'Angle of Attack, \alpha (°)';
    xlabel(figure_xlabel)
    figure_ylabel = 'Lift Coefficient, C_L';
    ylabel(figure_ylabel)
    legend('XFOIL Data', 'N = 200', 'N = 50', 'N = 100', Location='bestoutside')
    set(0, 'defaultfigureposition', [1300 10 800 400])
    saveas(gcf, 'C_L versus AoA for different cases of N against the XFOIL data.png')
    hold off
    %% plotting the stream lines around the airfoil using the streamlineplot():
    streamlineplot(airfoilcode, N(end), deg2rad(AoA_XFOIL(end)), U_infinity);
end
