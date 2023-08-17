% GOPH 547 : Lab 3
% Fazlie Latib
% 30067991
%% 

% Load data given for Lab 3

load goph_547_lab3_data.mat

whos -file goph_547_lab3_data.mat
%% 

% Q1 - 2D contour plot of the raw Fz data and plots of the raw Fz data against x and y coordinates

figure();
subplot(2,2,1);
contourf(X,Y,Fz_raw);
axis equal;
grid on;
title('Raw total F_z field data');
xlabel('x [m]');
ylabel('y [m]');
c1 = colorbar;
colormap(parula);
ylabel(c1,'F_z [\gamma = nT]','fontweight','bold');

subplot(2,2,3);
plot(X,Fz_raw,'o',color = 'blue');
xlabel('x [m]');
ylabel('F_z [\gamma = nT]');

subplot(2,2,2);
plot(Fz_raw,Y,'o',color = 'blue');

% Q2 - Fit 1st order polynomial and add the line to Fz vs y data

a_y = polyfit(Y,Fz_raw,1);
linear_x = min(Y(:)):1:max(Y(:))+1;
linear_y = a_y(1) * linear_x + a_y(2);
hold on;
plot(linear_y,linear_x,'--',color = 'red');
xlim([60 90])
ylim([0 300])
xlabel('F_z [\gamma = nT]');
ylabel('y [m]');
caption1 = sprintf('F_z = %f + %f y', a_y(2), a_y(1));
text(62,250,caption1,"FontSize",10);
%% 

% Q3 - Remove the linear component of the regional variation in the y direction from the data

Fz = Fz_raw - a_y(1) * Y; % remove linear component of regional variation in y direction

figure();
subplot(2,2,1);
contourf(X,Y,Fz);
axis equal;
grid on;
title('Data with y-component of regional removed');
xlabel('x [m]');
ylabel('y [m]');
c1 = colorbar;
colormap(parula);
ylabel(c1,'F_z [\gamma = nT]','fontweight','bold');

subplot(2,2,2);
plot(Fz,Y,'o',color = 'blue');
title('');
xlim([60 70])
ylim([0 300])
xlabel('F_z [\gamma = nT]');
ylabel('y [m]');

subplot(2,2,3);
plot(X,Fz,'o',color = 'blue');

% Fit 1st order polynomial and add the line to Fz vs x data

a_x = polyfit(X,Fz,1);
linear_x = min(X(:)):1:max(X(:))+1;
linear_y = a_x(1) * linear_x + a_x(2);
hold on;
plot(linear_x,linear_y,'--',color = 'red');
title('');
xlabel('x [m]');
ylabel('F_z [\gamma = nT]');
xlim([0 450])
ylim([60 70])
caption1 = sprintf('F_z = %f + %f x', a_x(2), a_x(1));
text(40,62,caption1,"FontSize",10);
%% 

% Q4 - Remove the linear component of the regional variation in the x direction from the data 

Fz = Fz - a_x(1) * X; % remove linear component of regional variation in x direction

figure();
subplot(2,2,1);
contourf(X,Y,Fz);
axis equal;
grid on;
title('Data with x- and y-component of regional removed');
xlabel('x [m]');
ylabel('y [m]');
c1 = colorbar;
colormap(parula);
ylabel(c1,'F_z [\gamma = nT]','fontweight','bold');

subplot(2,2,2);
plot(Fz,Y,'o',color = 'blue');
title('');
xlim([62 70])
ylim([0 300])
xlabel('F_z [\gamma = nT]');
ylabel('y [m]');

subplot(2,2,3);
plot(X,Fz,'o',color = 'blue');
title('');
xlabel('x [m]');
ylabel('F_z [\gamma = nT]');
ylim([62 70])
%% 

% Q5 - Remove the constant component of the regional

Fz = Fz - min(Fz(:)); % subtract min Fz

figure();
contourf(X,Y,Fz); 
axis equal;
grid on;
c1 = colorbar;
title('Data with all regional components removed')
xlabel('X [m]');
ylabel('Y [m]');
ylabel(c1,'F_z [\gamma = nT]','fontweight','bold');
%% 

% Q6 - Upward continue by performing double integration

dx  = 30;
dy  = dx;
dA  = dx * dy;
h   = 30;

Fz_u = zeros(size(Fz));

for i = 1:length(Fz(:,1)) % use for loops for every column/row
    for j = 1:length(Fz(1,:))
        x_30 = [X(i,j),Y(i,j),-h];  % specify the x and y points in our data set
        for n = 1:length(Fz(:,1))
            for m = 1:length(Fz(1,:)) % specify to run loop for all Fz data values
                x_o = [X(n,m),Y(n,m),0]; % the integration itself
                r = sqrt((x_o(1) - x_30(1))^2 + (x_o(2) - x_30(2))^2 + h^2); % the integration itself
                up_Fz = (h/(2*pi)) * (Fz(n,m)/(r^3)) * dA; % the integration itself
                Fz_u(i,j) = up_Fz + Fz_u(i,j); % apply upward continuatiion    
            end
        end
    end
end

figure();
contourf(X,Y,Fz_u); 
axis equal;
grid on;
c1 = colorbar;
title('Reduced data upward continued by 30 m')
xlabel('X [m]');
ylabel('Y [m]');
ylabel(c1,'F_z [\gamma = nT]','fontweight','bold');
%% 

% Q7 - Downward continue by performing finite difference approximation
Fz_d = Fz;

for j3 = 2:size(Fz_d,1)-1 % here we have an algorithm to implement downward continuation
    for i3 = 2:size(Fz_d,2)-1 % specify the column and of what Fz values to 'hit'
        Fz_d(j3,i3) = 6 * Fz(j3,i3) - (Fz(j3-1,i3) + Fz(j3+1,i3) + Fz(j3,i3-1) + Fz(j3,i3+1) + Fz_u(j3,i3)); % applying given formula 
    end
end

figure();
contourf(X,Y,Fz_d); 
axis equal;
grid on;
c1 = colorbar;
title('Reduced data downward continued by 30 m')
xlabel('X [m]');
ylabel('Y [m]');
ylabel(c1,'F_z [\gamma = nT]','fontweight','bold');
