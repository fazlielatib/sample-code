% GOPH 547 : Lab 2 - Part 2
% Fazlie Latib
% 30067991
%% 
% Load data from Lab 2 Part 1
load g_raw.mat
load g_norm.mat
load g_tide.mat
load g_drft.mat
load g_corr.mat
load z_datum.mat
load Xg.mat
load Yg.mat
load Zg.mat

%% Q15 - Calculate free air correction

dz = Zg - z_datum;  
dg_FA = (-0.3086 * dz) * 10^3;  % convert mGal to uGal
dg_FA = reshape(dg_FA,[Nx,Ny]);
g_corr = g_corr - dg_FA;

%% Q16 - Create contour plot of g_FA

g_FA = g_corr;

figure();
contourf(Xg,Yg,g_FA);
axis equal;
title('Raw - Norm - Drift - Free Air');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%% Q17 - Calculate Bouguer plate correction

rho_B = 2.65 * 1000;
dg_BP=(0.04193 * rho_B * dz) * 10^3; %convert density from g/cc to kg/m^3
dg_BP((Zg - z_datum) > 0) = 0; 
dg_BP = reshape(dg_BP,[Nx,Ny]);
g_corr = g_corr - dg_BP;

%% Q18 - Create contour plot of g_elev

g_elev = g_corr;

figure();
contourf(Xg,Yg,g_elev);
axis equal;
title('Raw - Norm - Drift - Elevation');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%% Q19 - Calculate terrain corection

dg_terr = zeros(size(g_corr));
dX = 1000; dY = 1000;
dA = dX * dY;

G = 6.674e-11;

for i = 1:numel(dg_terr)
    xi = [Xg(i) * 1000,Yg(i) * 1000,z_datum];
    for j = 1:numel(Xt)
           xm = [Xt(j) * 1000,Yt(j) * 1000, 0.5 * (Zt(j) + z_datum)];
           dm = rho_B * (Zt(j) - z_datum) * dA;
           dg_terr(i) = dg_terr(i) + abs(grav_eff_point(xi,xm,dm,G));
    end
end

dg_terr = dg_terr * 10 ^ 8; % convert unit from m/s^2 to uGal
g_corr = g_corr + dg_terr;

%% Q20 - Create contour plot of g_terr

g_terr = g_corr;

figure();
contourf(Xg,Yg,g_terr);
axis equal;
title('Raw - Norm - Drift - Boueguer');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%% Q21 - Calculate regional corection and its plot

dg_rgnl = mean(mean(g_corr));
g_corr = g_corr - dg_rgnl;
g_anom = g_corr;

figure();
contourf(Xg,Yg,g_anom);
axis equal;
title('Raw - Norm - Drift - Boueguer - Regional');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

% 
figure();
subplot(3,3,1)
contourf(Xg,Yg,g_raw);
axis equal;
title('Raw');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,2)
contourf(Xg,Yg,g_drft);
axis equal;
title('Raw - Norm - Drift');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,3)
contourf(Xg,Yg,g_terr);
axis equal;
title('Raw - Norm - Drift - Elev - Terr');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,4)
contourf(Xg,Yg,g_norm);
axis equal;
title('Raw - Norm');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,5)
contourf(Xg,Yg,g_FA);
axis equal;
title('Raw - Norm - Drift - FA');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,6)
contourf(Xg,Yg,g_anom);
axis equal;
title('Raw - Norm - Drift - Elev - Terr - Rgnl');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,7)
contourf(Xg,Yg,g_tide);
axis equal;
title('Raw - Norm - Tide');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(3,3,8)
contourf(Xg,Yg,g_elev);
axis equal;
title('Raw - Norm - Drift - Elev');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%% Q22 - Create contour plots of dg/dx, dg/dy, d^2g/dz^2

dg_dx = zeros(size(g_anom));
dg_dy = zeros(size(g_anom));
d2g_dz2 = zeros(size(g_anom));
dg_dx(:,2:end) = (g_anom(:,2:end) - g_anom(:,1:end-1)) / dX;
dg_dy(2:end,:) = (g_anom(2:end,:)-g_anom(1:end-1,:)) / dY;
d2g_dz2(2:end-1,2:end-1) = (2 * (dX^2 + dY^2) * g_anom(2:end-1,2:end-1) - (dX^2*(g_anom(3:end,2:end-1) + g_anom(1:end-2,2:end-1)) + dY^2 * (g_anom(2:end-1,3:end) + g_anom(2:end-1,1:end-2)))) / (dX^2 * dY^2);

figure();
subplot(1,3,1)
contourf(Xg,Yg,dg_dx);
axis equal;
title('dg/dx');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
colorbar;
colormap(flipud(parula));

subplot(1,3,2)
contourf(Xg,Yg,dg_dy);
axis equal;
title('dg/dy');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
colorbar;
colormap(flipud(parula));

subplot(1,3,3)
contourf(Xg,Yg,d2g_dz2);
axis equal;
title('d^2g/dz^2');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
colorbar;
colormap(flipud(parula));

% estimate the total mass of the anomalies below the points of highest concentration of gravity effect
[xi_anom,yi_anom] = find(g_anom == max(max(g_anom)));
max_anom = g_anom(xi_anom,yi_anom);
x_anom = Xg(xi_anom,yi_anom);
y_anom = Yg(xi_anom,yi_anom);
masslocation = [x_anom,y_anom,z_datum-500];
slocation = [x_anom,y_anom,z_datum];
mass = max_anom * 10^(-8) * norm(masslocation - slocation)^3 / (G * 500);
disp('The total mass of anomalies is around:');
fprintf('%.4e kg\n',mass);
%% 
function [gz] = grav_eff_point(x,xm,m,G)
%
% This function computes the gravity effect due to a mass anomaly
%
% INPUTS:
% x is the investigation point coordinate
% xm is the mass anomaly coordinate
% m is the anomaly mass
% G is the univerdal gravitational constant
%
% OUTPUT:
% gz is the gravity effect at coordinate xm due to the mass anomaly
% m

    s = 0;
    
    for i = 1:3
        s = s + ( x(i) - xm(i) ) ^ 2;
    end
    
    r = sqrt(s); % distance from x to xm

    gz = G * m * (x(3) - xm(3)) / r ^ 3;
end

