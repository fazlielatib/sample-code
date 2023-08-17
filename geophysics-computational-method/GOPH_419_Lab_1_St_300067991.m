%% *GOPH 419 Lab 1* 
%% *Fazlie Latib    30067991*
% 
%% Verifying the results from the launch_angle function at ve/v0 = 2, alpha = 0.25 and tol_alpha = 0.02

[expected_min_angle,expected_max_angle] = verify_launch_angle(2,0.25,0.02) % calculate the expected results using verify_launch_angle function

[min_angle,max_angle] = launch_angle(2,0.25,0.02) % compare the expected results with the real results using launch_angle function

%% Plotting graph of launch angles over a range of alpha with 4% alpha tolerance

range_alpha = 0:0.001:0.35; % set the range of alpha

min_angle = zeros(length(range_alpha),1);
max_angle = zeros(length(range_alpha),1);

% begin loop
for ii = 1:length(range_alpha)
    
    % estimate the minimum and maximum angles over the range
    [min_angle(ii),max_angle(ii)] = launch_angle(2,range_alpha(ii),0.04);
    
end

figure; hold on
plot(range_alpha,min_angle, 'r-') % plot minimum angles over range of alpha
plot(range_alpha,max_angle, 'b-') % plot maximum angles over range of alpha
xlim([0 0.45]) % set x-axis limits
ylim([0 pi/2]) % set y-axis limits
title('Launch angle over a range of alpha') % set the graph title
xlabel('alpha') % set the x-axis label
ylabel('launch angle (radian)') % set the y-axis label
legend('Minimum launch angle','Maximum launch angle','Location','east') % set the legends and its location
%% 
%% Plotting launch angles over a range of ve/v0 with 4% alpha tolerance

range_ve_v0 = 1.3:0.001:2.28; % set the range of alpha

min_angle = zeros(length(range_ve_v0),1);
max_angle = zeros(length(range_ve_v0),1);

% begin loop
for iii = 1:length(range_ve_v0)
    
     % estimate the minimum and maximum angles over the range
    [min_angle(iii),max_angle(iii)] = launch_angle(range_ve_v0(iii),0.25,0.04);
end

figure; hold on
plot(range_ve_v0,min_angle, 'r-') % plot minimum angles over range of alpha
plot(range_ve_v0,max_angle, 'b-') % plot maximum angles over range of alpha
xlim([0 4]) % set x-axis limits
ylim([0 pi/2]) % set y-axis limits
title('Launch angle over a range of ve/v0') % set the graph title
xlabel('ve/v0') % set the x-axis label
ylabel('launch angle (radian)') % set the y-axis label
legend('Minimum launch angle','Maximum launch angle','Location','east') % set the legends and its location
%% 
%% Calculating the error in sin phi

syms alpha ve_v0 % use symbolic functions

sin_phi = (1 + alpha)*sqrt( 1 - (alpha/(1 + alpha))*(ve_v0)^2 ); % state the sin phi equation

dfalpha = diff(sin_phi,alpha) % differentiate sin phi with respect to alpha
dfve_v0 = diff(sin_phi,ve_v0) % differentiate sin phi with respect to ve/v0
alpha = 0.25 % set value of alpha
ve_v0 = 2 % set value of ve/v0
subs(dfalpha) % calculate value of dfalpha
subs(dfve_v0) % calculate value of dfve_v0
error_sin_phi = abs(subs(dfalpha))*0.02 + abs(subs(dfve_v0))*0.05 % calculate error in sin phi
%% 
%% Calculating the condition number at ve/v0 = 2 and alpha = 0.25

cn_x = sqrt(alpha^2+ve_v0^2) % calculate values of ||x||
cn_j = sqrt(eval(dfalpha)^2+eval(dfve_v0)^2) % calculate values of ||J||
cn_f = eval(sin_phi) % calculate values of ||f||
cn = (cn_x*cn_j)/cn_f % calculate the condition number

%%
function [min_angle,max_angle] = launch_angle(ve_v0,alpha,tol_alpha)
%
% This function uses the Maclaurin series approximation for (sin-1 x)^2 
% which is the square of the rocket's launch angle and evaluates it at ve_v0, 
% alpha and tol_alpha to compute the allowable range of launch angles in radian.
%
% INPUTS:
% ve_v0 is the ratio of escaped velocity to terminal velocity
% alpha is the desired maximum altitude as a fraction of Earth's radius
% tol_alpha is the tolerance for the maximum altitude
%
% OUTPUT:
% min_angle is the approximate minimum allowable launch angle and max_angle
% is the approximate maximum launch angle at ve_v0, alpha and tol_alpha to
% 5 significant figures.

max_altitude = (1 + tol_alpha)*alpha; % calculate the maximum altitude for minimum launch angle
min_altitude = (1 - tol_alpha)*alpha; % calculate the minimum altitude for maximum launch angle

x_max_altitude = (1 + max_altitude)*sqrt( 1 - (max_altitude/(1 + max_altitude))*(ve_v0)^2 ); % x-value that we want to approximate minimum launch angle 
x_min_altitude = (1 + min_altitude)*sqrt( 1 - (min_altitude/(1 + min_altitude))*(ve_v0)^2 ); % x-value that we want to approximate maximum launch angle 

ea_min = 1; % define initial approximate relative error for minimum angle approximation
ea_max = 1; % define initial approximate relative error for maximum angle approximation
es = 0.5e-5; % define stopping criterion

min_angle_square = 0; % initiate Maclaurin series for minimum angle approximation
min_angle_prev = 0;
max_angle_square = 0; % initiate Maclaurin series for maximum angle approximation
max_angle_prev = 0;

n = 1; % define iteration number for minimum angle approximation
m = 1; % define iteration number for maximum angle approximation

% begin loop for minimum angle approximation
while ea_min > es
    
    % estimate square of the minimum angle
    min_angle_square = min_angle_square + 0.5*(2*x_max_altitude).^(2*n) ./ ( (n.^2)*((factorial(2*n)) / (factorial(n)).^2) );
    
    n = n + 1; % count the number of iterations
    
    % calculate approximate relative error
    ea_min = abs((min_angle_square - min_angle_prev) ./ min_angle_square);
    
    % store value of function for next iteration
    min_angle_prev = min_angle_square;
    
    % calculate minimum angle approximation
    min_angle = sqrt(min_angle_prev);
end

% begin loop for maximum angle approximation
while ea_max > es
    
    % estimate square of the maximum angle
    max_angle_square = max_angle_square + 0.5*(2*x_min_altitude).^(2*m) ./ ( (m.^2)*((factorial(2*m)) / (factorial(m)).^2) );
    
    m = m + 1; % count the number of iterations
    
    % calculate approximate relative error
    ea_max = abs((max_angle_square - max_angle_prev) ./ max_angle_square);
    
    % store value of function for next iteration
    max_angle_prev = max_angle_square;
    
    % calculate minimum angle approximation
    max_angle = sqrt(max_angle_prev);

end

end


function [expected_min_angle,expected_max_angle] = verify_launch_angle(ve_v0,alpha,tol_alpha)
%
% This function uses inverse sine function to compute the allowable 
% range of launch angles in radian at ve_v0, alpha and tol_alpha to
% verify the results of launch_angle function.
%
% INPUTS:
% ve_v0 is the ratio of escaped velocity to terminal velocity
% alpha is the desired maximum altitude as a fraction of Earth's radius
% tol_alpha is the tolerance for the maximum altitude
%
% OUTPUT:
% min_angle is the approximate minimum allowable launch angle and max_angle
% is the approximate maximum launch angle at ve_v0, alpha and tol_alpha

max_altitude = (1 + tol_alpha)*alpha; % calculate the maximum altitude for minimum launch angle
min_altitude = (1 - tol_alpha)*alpha; % calculate the minimum altitude for maximum launch angle

% calculate minimum angle approximation
expected_min_angle = asin((1 + max_altitude)*sqrt( 1 - (max_altitude/(1 + max_altitude))*(ve_v0)^2 ));

% calculate maximum angle approximation
expected_max_angle = asin((1 + min_altitude)*sqrt( 1 - (min_altitude/(1 + min_altitude))*(ve_v0)^2 ));

end