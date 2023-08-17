% GOPH 547 : Lab 1
% Fazlie Latib
% 30067991

m = 1 * 10 ^ 7; % point mass
xm = [0,0,-10]; % point mass coordinate
G = 6.674e-11; % universal constant of gravitation

grid_space1 = 5; % 5 m grid spacings
grid_space2 = 10; % 10 m grid spacings 
grid_space3 = 25; % 25 m grid spacings

x_range1 = -100:grid_space1:100; % 5 m interval from -100 to 100 m
x_range2 = -100:grid_space2:100; % 10 m interval from -100 to 100 m
x_range3 = -100:grid_space3:100; % 25 m interval from -100 to 100 m

[X1,Y1] = meshgrid(x_range1); % grid with 5 m spacing
[X2,Y2] = meshgrid(x_range2); % grid with 10 m spacing
[X3,Y3] = meshgrid(x_range3); % grid with 25 m spacing

z1 = 0; % z = 0 m
z2 = 10; % z = 10 m
z3 = 100; % z = 100 m

% empty contour plot for U and gz with 5 m grid spacings
U1_grid1 = zeros(size(X1));
gz1_grid1 = zeros(size(X1));
U2_grid1 = zeros(size(X1));
gz2_grid1 = zeros(size(X1));
U3_grid1 = zeros(size(X1));
gz3_grid1 = zeros(size(X1));

% calculate U and gz at different z values with 5 m grid spacings
for i = 1:length(X1)
    for ii = 1:length(Y1)
        x = [X1(i,i),Y1(ii,ii),z1];
        [U] = grav_pot_point(x,xm,m,G);
        U1_grid1((length(x_range1)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz1_grid1((length(x_range1)+1)-ii,i) = gz;

        x = [X1(i,i),Y1(ii,ii),z2];
        [U] = grav_pot_point(x,xm,m,G);
        U2_grid1((length(x_range1)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz2_grid1((length(x_range1)+1)-ii,i) = gz;

        x = [X1(i,i),Y1(ii,ii),z3];
        [U] = grav_pot_point(x,xm,m,G);
        U3_grid1((length(x_range1)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz3_grid1((length(x_range1)+1)-ii,i) = gz;
    end
end

% empty contour plot for U and gz with 10 m grid spacings
U1_grid2 = zeros(size(X2));
gz1_grid2 = zeros(size(X2));
U2_grid2 = zeros(size(X2));
gz2_grid2 = zeros(size(X2));
U3_grid2 = zeros(size(X2));
gz3_grid2 = zeros(size(X2));

% calculate U and gz at different z values with 10 m grid spacings
for i = 1:length(X2)
    for ii = 1:length(Y2)
        x = [X2(i,i),Y2(ii,ii),z1];
        [U] = grav_pot_point(x,xm,m,G);
        U1_grid2((length(x_range2)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz1_grid2((length(x_range2)+1)-ii,i) = gz;

        x = [X2(i,i),Y2(ii,ii),z2];
        [U] = grav_pot_point(x,xm,m,G);
        U2_grid2((length(x_range2)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz2_grid2((length(x_range2)+1)-ii,i) = gz;

        x = [X2(i,i),Y2(ii,ii),z3];
        [U] = grav_pot_point(x,xm,m,G);
        U3_grid2((length(x_range2)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz3_grid2((length(x_range2)+1)-ii,i) = gz;
    end
end

% empty contour plot for U and gz with 25 m grid spacings
U1_grid3 = zeros(size(X3));
gz1_grid3 = zeros(size(X3));
U2_grid3 = zeros(size(X3));
gz2_grid3 = zeros(size(X3));
U3_grid3 = zeros(size(X3));
gz3_grid3 = zeros(size(X3));

% calculate U and gz at different z values with 25 m grid spacings
for i = 1:length(X3)
    for ii = 1:length(Y3)
        x = [X3(i,i),Y3(ii,ii),z1];
        [U] = grav_pot_point(x,xm,m,G);
        U1_grid3((length(x_range3)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz1_grid3((length(x_range3)+1)-ii,i) = gz;

        x = [X3(i,i),Y3(ii,ii),z2];
        [U] = grav_pot_point(x,xm,m,G);
        U2_grid3((length(x_range3)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz2_grid3((length(x_range3)+1)-ii,i) = gz;

        x = [X3(i,i),Y3(ii,ii),z3];
        [U] = grav_pot_point(x,xm,m,G);
        U3_grid3((length(x_range3)+1)-ii,i) = U;
        [gz] = grav_eff_point(x,xm,m,G);
        gz3_grid3((length(x_range3)+1)-ii,i) = gz;
    end
end

% plots of U and gz for 5 m grid spacings
figure();
subplot(3,2,1);
contourf(X1,Y1,U1_grid1);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 0 m')
c1 = colorbar;
c1.Limits = [0 6.5e-5];
c1.Label.String = 'U [J/kg = m^2/s^2]';
color = flipud(parula);
colormap(color);

subplot(3,2,2);
contourf(X1,Y1,gz1_grid1/1e-5);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 0 m')
c2 = colorbar;
c2.Limits = [0 0.65];
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';
color = flipud(parula);
colormap(color);

subplot(3,2,3);
contourf(X1,Y1,U2_grid1);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 10 m')
c1 = colorbar;
caxis([0 6.5e-5]);
c1.Label.String = 'U [J/kg = m^2/s^2]';

subplot(3,2,4);
contourf(X1,Y1,gz2_grid1/1e-5);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 10 m')
c2 = colorbar;
caxis([0 0.65]);
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';

subplot(3,2,5);
contourf(X1,Y1,U3_grid1);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 100 m')
c1 = colorbar;
caxis([0 6.5e-5]);
c1.Label.String = 'U [J/kg = m^2/s^2]';

subplot(3,2,6);
contourf(X1,Y1,gz3_grid1/1e-5);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 100 m')
c2 = colorbar;
caxis([0 0.65]);
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';

% plots of U and gz for 10 m grid spacings
figure();
subplot(3,2,1);
contourf(X2,Y2,U1_grid2);
hold on
plot(X2,Y2,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 0 m')
c1 = colorbar;
c1.Limits = [0 6.5e-5];
c1.Label.String = 'U [J/kg = m^2/s^2]';
color = flipud(parula);
colormap(color);

subplot(3,2,2);
contourf(X2,Y2,gz1_grid2/1e-5);
hold on
plot(X2,Y2,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 0 m')
c2 = colorbar;
c2.Limits = [0 0.65];
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';
color = flipud(parula);
colormap(color);

subplot(3,2,3);
contourf(X2,Y2,U2_grid2);
hold on
plot(X2,Y2,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 10 m')
c1 = colorbar;
caxis([0 6.5e-5]);
c1.Label.String = 'U [J/kg = m^2/s^2]';

subplot(3,2,4);
contourf(X2,Y2,gz2_grid2/1e-5);
hold on
plot(X2,Y2,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 10 m')
c2 = colorbar;
caxis([0 0.65]);
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';

subplot(3,2,5);
contourf(X2,Y2,U3_grid2);
hold on
plot(X2,Y2,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 100 m')
c1 = colorbar;
caxis([0 6.5e-5]);
c1.Label.String = 'U [J/kg = m^2/s^2]';

subplot(3,2,6);
contourf(X2,Y2,gz3_grid2/1e-5);
hold on
plot(X2,Y2,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 100 m')
c2 = colorbar;
caxis([0 0.65]);
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';

% plots of U and gz for 25 m grid spacings
figure();
subplot(3,2,1);
contourf(X3,Y3,U1_grid3);
hold on
plot(X3,Y3,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 0 m')
c1 = colorbar;
c1.Limits = [0 6.5e-5];
c1.Label.String = 'U [J/kg = m^2/s^2]';
color = flipud(parula);
colormap(color);

subplot(3,2,2);
contourf(X3,Y3,gz1_grid3/1e-5);
hold on
plot(X3,Y3,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 0 m')
c2 = colorbar;
c2.Limits = [0 0.65];
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';
color = flipud(parula);
colormap(color);

subplot(3,2,3);
contourf(X3,Y3,U2_grid3);
hold on
plot(X3,Y3,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 10 m')
c1 = colorbar;
caxis([0 6.5e-5]);
c1.Label.String = 'U [J/kg = m^2/s^2]';

subplot(3,2,4);
contourf(X3,Y3,gz2_grid3/1e-5);
hold on
plot(X3,Y3,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 10 m')
c2 = colorbar;
caxis([0 0.65]);
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';

subplot(3,2,5);
contourf(X3,Y3,U3_grid3);
hold on
plot(X3,Y3,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 100 m')
c1 = colorbar;
caxis([0 6.5e-5]);
c1.Label.String = 'U [J/kg = m^2/s^2]';

subplot(3,2,6);
contourf(X3,Y3,gz3_grid3/1e-5);
hold on
plot(X3,Y3,'xk','MarkerSize',2);
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 100 m')
c2 = colorbar;
caxis([0 0.65]);
c2.Label.String = 'g [mGal = 10^{-3} cm/s^2]';
%% 

function [U] = grav_pot_point(x,xm,m,G)
%
% This function computes the gravitational potential due to a mass anomaly
%
% INPUTS:
% x is the investigation point coordinate
% xm is the mass anomaly coordinate
% m is the anomaly mass
% G is the univerdal gravitational constant
%
% OUTPUT:
% U is the gravitational potential at coordinate xm due to the mass anomaly
% m

    s = 0;
    
    for i = 1:3
        s = s + ( x(i) - xm(i) ) ^ 2;
    end
    
    r = sqrt(s); % distance from x to xm
    
    U = G * m / r;
end

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