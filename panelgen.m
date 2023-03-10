function [X,Z] = panelgen(airfoilcode, N, angleofattack)
% The panelgen() function can generate and discretize any NACA 4-series 
% airfoil using Panel Code. This function will also calculate the panel 
% end points and return two one-dimensional arrays 'X' and 'Z'.
% -------------------------------------------------------------------------
% inputs:
% airfoilcode - A four digit code that refers to the NACA airfoil e.g.
% 1410.
% N - The number of panels to be used e.g. 100.
% angleofattack - The angle of attack in degrees e.g. 10°.
% -------------------------------------------------------------------------
% outputs:
% X - An one-dimensional array with the x-coordinates of the discretised
% airfoil
% Z - An one-dimensional array with the z-coordinates of the discretised
% airfoil

%% convert user input 'airfoilcode' into a string and 'angleofattack' to radians:
aoa = deg2rad(angleofattack);

%% extract the m, p and t values from 'airfoilcodestr':
m = str2double(airfoilcode(1))*0.01;     % max camber in % of chord
p = str2double(airfoilcode(2))*0.1;      % distance of maximum camber from the airfoil leading edge in tenths
t = str2double(airfoilcode([3,4]))*0.01; % max thickness

%% finding the panel end point, x, using a cosine distribution, mean camber line gradient, dyc_dx, and theta:
% initialise variables with arrays full of zeros of dimension N+1 by 1 
x = zeros(N+1, 1);
y_c = zeros(N+1, 1);
dyc_dx = zeros(N+1, 1);
theta = zeros(N+1, 1);
y_t = zeros(N+1,1);

for i = 1:N+1
    x(i) = 1 - 0.5*(1 - cos(2*pi*(i-1)/N));
    if (x(i) >= 0 && x(i) < p)
        dyc_dx(i)=(2*m/p^2)*(p-x(i));
        y_c(i) = (m/p^2)*(2*p*x(i) - x(i)^2);
    elseif (x(i) > p && x(i) <= 1)
        dyc_dx(i)=(2*m/(p-1)^2)*(p-x(i));
        y_c(i) = (m/(1-p)^2)*((1-2*p)+2*p*x(i)-x(i)^2);
    end
    theta(i)=atan(dyc_dx(i));
    y_t(i) = 5*t*(0.2969*x(i)^0.5-0.126*x(i)-0.3516*x(i)^2+0.2843*x(i)^3-0.1015*x(i)^4);
end

%% to close the trailing edge of the airfoil:
y_t(1)=0;
y_t(end)=0;

%% evaluating the coordinates of upper and lower surface panel endpoints:
% upper surface
x_U = x(:) - y_t(:).*sin(theta);
z_U = y_c(:) + y_t(:).*cos(theta);
% lower surface
x_L = x(:) + y_t(:).*sin(theta);
z_L = y_c(:) - y_t(:).*cos(theta);

%% finally, generating X and Z arrays from x_U, z_U, x_L and z_L:

zerothelement = find(x_L==0); % finding the 0th element in x_L:

% initialise the variables:
X = zeros(1, N+1);
Z = zeros(1, N+1);
for i=1:length(x_L)
    if i<=zerothelement
        X(i)=x_L(i);
        Z(i)=z_L(i);
    else
        X(i)=x_U(i);
        Z(i)=z_U(i);
    end
end

%% mimic the wake of an airfoil using the wake panel at the endpoint of X and Z arrays:
wake_x = 10;
X(length(X)+1) = wake_x * cos(-aoa);
Z(length(Z)+1) = wake_x * sin(-aoa);
end
