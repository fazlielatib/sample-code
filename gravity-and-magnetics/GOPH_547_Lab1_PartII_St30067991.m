% GOPH 547 : Lab 1 - Part 2
% Fazlie Latib
% 30067991

%% Question 4a and 4b

m = 1e7; % total mass of 5 point anomalies

% mean of masses and their locations
mean_m = m / 5; mean_x = 0; mean_y = 0; mean_z = -10;

% standard deviation of masses and their locations
sd_m = m / 100; sd_x = 20; sd_y = 20; sd_z = 2;

% generate 4 random masses and their locations (Set 1)
m_rand = sd_m .* randn(4,1) + mean_m;
x_rand = sd_x .* randn(4,1) + mean_x;
y_rand = sd_y .* randn(4,1) + mean_y;
z_rand = sd_z .* randn(4,1) + mean_z;

% calculate 5th mass and its location for Set 1
m5 = m - sum(m_rand);
x5 = - sum(m_rand.*x_rand) / m5;
y5 = - sum(m_rand.*y_rand) / m5;
z5 = (mean_z*m - sum(m_rand.*z_rand)) / m5;
m_rand(end+1) = m5; x_rand(end+1) = x5; y_rand(end+1) = y5; z_rand(end+1) = z5;

save mass_set_3 % save data for Set 1

% verify Set 1 satisfy the specified constraints

load mass_set_3.mat

mean_sd_m = [mean_m sd_m]
mean_sd_m_set1 = [mean(m_rand) std(m_rand)]
mean_sd_x = [mean_x sd_x]
mean_sd_x_set1 = [mean(x_rand) std(x_rand)]
mean_sd_y = [mean_y sd_y]
mean_sd_y_set1 = [mean(y_rand) std(y_rand)]
mean_sd_z = [mean_z sd_z]
mean_sd_z_set1 = [mean(z_rand) std(z_rand)]
z_rand <= -1

%% Question 4c

G = 6.674e-11; % universal constant of gravitation

% set position of each individual random mass point
xm1 = [x_rand(1),y_rand(1),z_rand(1)];
xm2 = [x_rand(2),y_rand(2),z_rand(2)];
xm3 = [x_rand(3),y_rand(3),z_rand(3)];
xm4 = [x_rand(4),y_rand(4),z_rand(4)];
xm5 = [x_rand(5),y_rand(5),z_rand(5)];

% set mass at each individual random point  
m1 = m_rand(1);
m2 = m_rand(2);
m3 = m_rand(3);
m4 = m_rand(4);
m5 = m_rand(5);

% toggle between these 3 grid spacings
%dx = 5; % set spacings for 5 m
%dx = 10; % set spacings for 10 m 
dx = 25; % set spacings for 25 m

x_range = -100:dx:100; % set chosen grid spacings from -100 to 100 m

[X1,Y1] = meshgrid(x_range); % set plot with chosen spacings

z1 = 0; % set z = 0 m
z2 = 10; % set z = 10 m
z3 = 100; % set z = 100 m

% build empty contour plot for U and gz
U1 = zeros(size(X1));
gz1 = zeros(size(X1));
U2 = zeros(size(X1));
gz2 = zeros(size(X1));
U3 = zeros(size(X1));
gz3 = zeros(size(X1));

% calculate U and gz at different z values with chosen grid spacings
for i = 1:length(X1)
    for ii = 1:length(Y1)

        % calculate U and gz at z = 0
        U_p1 = grav_pot_point(([X1(i,i),Y1(ii,ii),z1]),xm1,m1,G);
        U_p2 = grav_pot_point(([X1(i,i),Y1(ii,ii),z1]),xm2,m2,G);
        U_p3 = grav_pot_point(([X1(i,i),Y1(ii,ii),z1]),xm3,m3,G);
        U_p4 = grav_pot_point(([X1(i,i),Y1(ii,ii),z1]),xm4,m4,G);
        U_p5 = grav_pot_point(([X1(i,i),Y1(ii,ii),z1]),xm5,m5,G);
        % sum contributions of U from each of the individual mass 
        U1((length(x_range)+1)-ii,i) = U_p1 + U_p2 + U_p3 + U_p4 + U_p5; 

        gz_p1 = grav_eff_point(([X1(i,i),Y1(ii,ii),z1]),xm1,m1,G);
        gz_p2 = grav_eff_point(([X1(i,i),Y1(ii,ii),z1]),xm2,m2,G);
        gz_p3 = grav_eff_point(([X1(i,i),Y1(ii,ii),z1]),xm3,m3,G);
        gz_p4 = grav_eff_point(([X1(i,i),Y1(ii,ii),z1]),xm4,m4,G);
        gz_p5 = grav_eff_point(([X1(i,i),Y1(ii,ii),z1]),xm5,m5,G);
        % sum contributions of gz from each of the individual masses
        gz1((length(x_range)+1)-ii,i) = gz_p1 + gz_p2 + gz_p3 + gz_p4 + gz_p5; 

        % calculate U and gz at z = 10
        U_p1 = grav_pot_point(([X1(i,i),Y1(ii,ii),z2]),xm1,m1,G);
        U_p2 = grav_pot_point(([X1(i,i),Y1(ii,ii),z2]),xm2,m2,G);
        U_p3 = grav_pot_point(([X1(i,i),Y1(ii,ii),z2]),xm3,m3,G);
        U_p4 = grav_pot_point(([X1(i,i),Y1(ii,ii),z2]),xm4,m4,G);
        U_p5 = grav_pot_point(([X1(i,i),Y1(ii,ii),z2]),xm5,m5,G);
        % sum contributions of U from each of the individual mass 
        U2((length(x_range)+1)-ii,i) = U_p1 + U_p2 + U_p3 + U_p4 + U_p5; 
        
        gz_p1 = grav_eff_point(([X1(i,i),Y1(ii,ii),z2]),xm1,m1,G);
        gz_p2 = grav_eff_point(([X1(i,i),Y1(ii,ii),z2]),xm2,m2,G);
        gz_p3 = grav_eff_point(([X1(i,i),Y1(ii,ii),z2]),xm3,m3,G);
        gz_p4 = grav_eff_point(([X1(i,i),Y1(ii,ii),z2]),xm4,m4,G);
        gz_p5 = grav_eff_point(([X1(i,i),Y1(ii,ii),z2]),xm5,m5,G);
        % sum contributions of gz from each of the individual masses
        gz2((length(x_range)+1)-ii,i) = gz_p1 + gz_p2 + gz_p3 + gz_p4 + gz_p5; 

        % calculate U and gz at z = 100
        U_p1 = grav_pot_point(([X1(i,i),Y1(ii,ii),z3]),xm1,m1,G);
        U_p2 = grav_pot_point(([X1(i,i),Y1(ii,ii),z3]),xm2,m2,G);
        U_p3 = grav_pot_point(([X1(i,i),Y1(ii,ii),z3]),xm3,m3,G);
        U_p4 = grav_pot_point(([X1(i,i),Y1(ii,ii),z3]),xm4,m4,G);
        U_p5 = grav_pot_point(([X1(i,i),Y1(ii,ii),z3]),xm5,m5,G);
        % sum contributions of U from each of the individual mass 
        U3((length(x_range)+1)-ii,i) = U_p1 + U_p2 + U_p3 + U_p4 + U_p5; 

        gz_p1 = grav_eff_point(([X1(i,i),Y1(ii,ii),z3]),xm1,m1,G);
        gz_p2 = grav_eff_point(([X1(i,i),Y1(ii,ii),z3]),xm2,m2,G);
        gz_p3 = grav_eff_point(([X1(i,i),Y1(ii,ii),z3]),xm3,m3,G);
        gz_p4 = grav_eff_point(([X1(i,i),Y1(ii,ii),z3]),xm4,m4,G);
        gz_p5 = grav_eff_point(([X1(i,i),Y1(ii,ii),z3]),xm5,m5,G);
        % sum contributions of gz from each of the individual masses
        gz3((length(x_range)+1)-ii,i) = gz_p1 + gz_p2 + gz_p3 + gz_p4 + gz_p5; 

    end
end

figure();
subplot(3,2,1);
contourf(X1,Y1,U1);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
axis equal;
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 0 m')
c1 = colorbar;
caxis([0 max(U1,[],'all')]);
ylabel(c1, 'U [J/kg = m^2/s^2]','fontweight','bold');
color = flipud(parula);
colormap(color);

subplot(3,2,2);
contourf(X1,Y1,gz1/1e-5);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
axis equal;
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 0 m')
c2 = colorbar;
caxis([0 max(gz1,[],'all')/1e-5]);
ylabel(c2, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');
color = flipud(parula);
colormap(color);

subplot(3,2,3);
contourf(X1,Y1,U2);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
axis equal;
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 10 m')
c1 = colorbar;
caxis([0 max(U1,[],'all')]);
ylabel(c1, 'U [J/kg = m^2/s^2]','fontweight','bold');

subplot(3,2,4);
contourf(X1,Y1,gz2/1e-5);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
axis equal;
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 10 m')
c2 = colorbar;
caxis([0 max(gz1,[],'all')/1e-5]);
ylabel(c2, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

subplot(3,2,5);
contourf(X1,Y1,U3);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
axis equal;
xlabel('x [m]')
ylabel('y [m]')
title('Gravitational potential at z = 100 m')
c1 = colorbar;
caxis([0 max(U1,[],'all')]);
ylabel(c1, 'U [J/kg = m^2/s^2]','fontweight','bold');

subplot(3,2,6);
contourf(X1,Y1,gz3/1e-5);
hold on
plot(X1,Y1,'xk','MarkerSize',2);
axis equal;
xlabel('x [m]')
ylabel('y [m]')
title('Gravity effect at z = 100 m')
c2 = colorbar;
caxis([0 max(gz1,[],'all')/1e-5]);
ylabel(c2, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

%% Question 5a

load anomaly_data.mat

rho = rho * 1000;

% find j parameters
jz1 = find(z(1,1,:) == -8,1); % when z = -8 m
jz2 = find(z(1,1,:) == -10,1); % when z = -10 m
jz3 = find(z(1,1,:) == -12,1); % when z = -12 m

jy1 = find(y(:,1,1) == -8,1); % when y = -8 m
jy2 = find(y(:,1,1) == 0,1); % when y = 0 m
jy3 = find(y(:,1,1) == 8,1); % when y = 8 m

jx1 = find(x(1,:,1) == -16,1); % when x = -16 m
jx2 = find(x(1,:,1) == 0,1); % when x = 0 m
jx3 = find(x(1,:,1) == 16,1); % when x = 16 m

% load x,y and z grid to plot corresponding rho values in cross-section
ny = size(x,1); % number of y points is the number of ROWS
nx = size(x,2); % number of x points is the number of COLUMNS
nz = size(z,3); % number of z points is the number of 'PAGES' (3rd 'axis' dimension)

% For z = -8 m
x_cur_z1 = zeros(nx,ny);    x_cur_z1(:,:) = x(:,:,jz1);
y_cur_z1 = zeros(nx,ny);    y_cur_z1(:,:) = y(:,:,jz1);
rho_cur_z1 = zeros(nx,ny);  rho_cur_z1(:,:) = rho(:,:,jz1);

% For z = -10 m
x_cur_z2 = zeros(nx,ny);    x_cur_z2(:,:) = x(:,:,jz2);
y_cur_z2 = zeros(nx,ny);    y_cur_z2(:,:) = y(:,:,jz2);
rho_cur_z2 = zeros(nx,ny);  rho_cur_z2(:,:) = rho(:,:,jz2);

% For z = -12 m
x_cur_z3 = zeros(nx,ny);    x_cur_z3(:,:) = x(:,:,jz3);
y_cur_z3 = zeros(nx,ny);    y_cur_z3(:,:) = y(:,:,jz3);
rho_cur_z3 = zeros(nx,ny);  rho_cur_z3(:,:) = rho(:,:,jz3);

% For y = -8 m
x_cur_y1 = zeros(nx,nz);    x_cur_y1(:,:) = x(jy1,:,:);
z_cur_y1 = zeros(nx,nz);    z_cur_y1(:,:) = z(jy1,:,:);
rho_cur_y1 = zeros(nx,nz);  rho_cur_y1(:,:) = rho(jy1,:,:);

% For y = 0 m
x_cur_y2 = zeros(nx,nz);    x_cur_y2(:,:) = x(jy2,:,:);
z_cur_y2 = zeros(nx,nz);    z_cur_y2(:,:) = z(jy2,:,:);
rho_cur_y2 = zeros(nx,nz);  rho_cur_y2(:,:) = rho(jy2,:,:);

% For y = 8 m
x_cur_y3 = zeros(nx,nz);    x_cur_y3(:,:) = x(jy3,:,:);
z_cur_y3 = zeros(nx,nz);    z_cur_y3(:,:) = z(jy3,:,:);
rho_cur_y3 = zeros(nx,nz);  rho_cur_y3(:,:) = rho(jy3,:,:);

% For x = -16 m
y_cur_x1 = zeros(ny,nz);    y_cur_x1(:,:) = y(:,jx1,:);
z_cur_x1 = zeros(ny,nz);    z_cur_x1(:,:) = z(:,jx1,:);
rho_cur_x1 = zeros(ny,nz);  rho_cur_x1(:,:) = rho(:,jx1,:);

% For x = 0 m
y_cur_x2 = zeros(ny,nz);    y_cur_x2(:,:) = y(:,jx2,:);
z_cur_x2 = zeros(ny,nz);    z_cur_x2(:,:) = z(:,jx2,:);
rho_cur_x2 = zeros(ny,nz);  rho_cur_x2(:,:) = rho(:,jx2,:);

% For x = 16 m
y_cur_x3 = zeros(ny,nz);    y_cur_x3(:,:) = y(:,jx3,:);
z_cur_x3 = zeros(ny,nz);    z_cur_x3(:,:) = z(:,jx3,:);
rho_cur_x3 = zeros(ny,nz);  rho_cur_x3(:,:) = rho(:,jx3,:);

% plot the assigned sections
figure(); % For z = -8 m
subplot(3,3,1)
contourf(x_cur_z1,y_cur_z1,rho_cur_z1);
c3 = colorbar; 
color = flipud(parula);
colormap(color);
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-20,20]);
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('Density anomaly map at z = -8 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,4) % For z = -10 m
contourf(x_cur_z2,y_cur_z2,rho_cur_z2);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-20,20]);
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('Density anomaly map at z = -10 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,7) % For z = -12 m
contourf(x_cur_z3,y_cur_z3,rho_cur_z3);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-20,20]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('y [m]', 'fontweight','bold');
title('Density anomaly map at z = -12 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,2) % For y = -8 m
contourf(x_cur_y1,z_cur_y1,rho_cur_y1);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); 
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map at y = -8 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,5) % For y = 0 m
contourf(x_cur_y2,z_cur_y2,rho_cur_y2);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); 
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map at y = 0 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,8) % For y = 8 m
contourf(x_cur_y3,z_cur_y3,rho_cur_y3);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); 
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map at y = 8 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,3) % For x = -16 m
contourf(y_cur_x1,z_cur_x1,rho_cur_x1);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-15,-5]);
xlabel('y [m]', 'fontweight','bold'); 
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map at x = -16 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,6)% For x = 0 m
contourf(y_cur_x2,z_cur_x2,rho_cur_x2);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-15,-5]);
xlabel('y [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map at x = 0 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

subplot(3,3,9) % For x = 16 m
contourf(y_cur_x3,z_cur_x3,rho_cur_x3);
c3 = colorbar; 
caxis([0 max(rho,[],'all')]);
axis equal;
axis([-20,20,-15,-5]);
xlabel('y [m]', 'fontweight','bold'); 
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map at x = 16 m', 'fontweight','bold');
ylabel(c3, '\rho [kg/m^3]', 'fontweight','bold');

%% Question 5b and 5c

% calculate volume of a small cube (2m x 2m x 2m)
dV = 2 * 2 * 2;

ind_nz = find(rho); % find non-zero rho
num_non_zero = length(ind_nz); % number of non-zero locations

dx1 = 5; % set spacings for 5 m
dx2 = 10; % set spacings for 10 m 
dx3 = 25; % set spacings for 25 m

x_range1 = -100:dx1:100; % set 5m grid spacings from -100 to 100 m
x_range2 = -100:dx2:100; % set 10m grid spacings from -100 to 100 m
x_range3 = -100:dx3:100; % set 25m grid spacings from -100 to 100 m

[X1,Y1] = meshgrid(x_range1); % set plot with 5m spacings
[X2,Y2] = meshgrid(x_range2); % set plot with 10m spacings
[X3,Y3] = meshgrid(x_range3); % set plot with 25m spacings

% build empty contour plot gz
forward_gz1a = zeros(size(X1)); forward_gz1b = zeros(size(X1));
forward_gz2a = zeros(size(X2)); forward_gz2b = zeros(size(X2));
forward_gz3a = zeros(size(X3)); forward_gz3b = zeros(size(X3));

Nx1 = size(X1,2); Nx2 = size(X2,2); Nx3 = size(X3,2); 
Ny1 = size(X1,1); Ny2 = size(X2,1); Ny3 = size(X3,1);

m_tot = 0; % set initial total mass to 0

for j = 1:Nx1
    for i = 1:Ny1
        for k = 1:num_non_zero
            kk = ind_nz(k);
            xm = [x(kk), y(kk), z(kk)]; % coordinates of non-zero
            dm = rho(kk)*dV;
            m_tot = m_tot + dm;
            forward_gz1a(i,j) = forward_gz1a(i,j) + grav_eff_point(([X1(i,j),Y1(i,j),z1]),xm,dm,G);
            forward_gz1b(i,j) = forward_gz1b(i,j) + grav_eff_point(([X1(i,j),Y1(i,j),z3]),xm,dm,G);
        end
    end
end

m_tot % total mass of density anomalies

for j = 1:Nx2
    for i = 1:Ny2
        for k = 1:num_non_zero
            kk = ind_nz(k);
            xm = [x(kk), y(kk), z(kk)]; % coordinates of non-zero
            dm = rho(kk)*dV;
            forward_gz2a(i,j) = forward_gz2a(i,j) + grav_eff_point(([X2(i,j),Y2(i,j),z1]),xm,dm,G);
            forward_gz2b(i,j) = forward_gz2b(i,j) + grav_eff_point(([X2(i,j),Y2(i,j),z3]),xm,dm,G);
        end
    end
end

for j = 1:Nx3
    for i = 1:Ny3
        for k = 1:num_non_zero
            kk = ind_nz(k);
            xm = [x(kk), y(kk), z(kk)]; % coordinates of non-zero
            dm = rho(kk)*dV;
            forward_gz3a(i,j)= forward_gz3a(i,j) + grav_eff_point(([X3(i,j),Y3(i,j),z1]),xm,dm,G);
            forward_gz3b(i,j)= forward_gz3b(i,j) + grav_eff_point(([X3(i,j),Y3(i,j),z3]),xm,dm,G);
        end
    end
end

% For z = 0 m
figure(); 
sgtitle('Forward Modelling Gravity Effect (g_z) with the Anomaly')
subplot(2,3,1)
contourf(X1,Y1,forward_gz1a/1e-5);
hold on; 
plot(X1,Y1,'xk','MarkerSize',2);
c4 = colorbar; 
colormap(parula);
caxis([0 max(forward_gz1a,[],'all')/1e-5]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 5 m', 'fontweight','bold');
ylabel(c4, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

subplot(2,3,2)
contourf(X2,Y2,forward_gz2a/1e-5);
hold on; 
plot(X2,Y2,'xk','MarkerSize',2);
c4 = colorbar; 
caxis([0 max(forward_gz1a,[],'all')/1e-5]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 10 m', 'fontweight','bold');
ylabel(c4, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

subplot(2,3,3)
contourf(X3,Y3,forward_gz3a/1e-5);
hold on; 
plot(X3,Y3,'xk','MarkerSize',2);
c4 = colorbar; 
caxis([0 max(forward_gz1a,[],'all')/1e-5]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 25 m', 'fontweight','bold');
ylabel(c4, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

% For z = 100 m
subplot(2,3,4)
contourf(X1,Y1,forward_gz1b/1e-5);
hold on; 
plot(X1,Y1,'xk','MarkerSize',2);
c4 = colorbar; 
colormap(parula);
caxis([0 max(forward_gz1b,[],'all')/1e-5]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 5 m', 'fontweight','bold');
ylabel(c4, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

subplot(2,3,5)
contourf(X2,Y2,forward_gz2b/1e-5);
hold on; 
plot(X2,Y2,'xk','MarkerSize',2);
c4 = colorbar; 
caxis([0 max(forward_gz1b,[],'all')/1e-5]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 10 m', 'fontweight','bold');
ylabel(c4, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

subplot(2,3,6)
contourf(X3,Y3,forward_gz3b/1e-5);
hold on; 
plot(X3,Y3,'xk','MarkerSize',2);
c4 = colorbar; 
caxis([0 max(forward_gz1b,[],'all')/1e-5]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 25 m', 'fontweight','bold');
ylabel(c4, 'g [mGal = 10^{-3} cm/s^2]','fontweight','bold');

%% Question 5d

% build empty contour plot gz
forward_gz1c = zeros(size(X1)); forward_gz1d = zeros(size(X1));
forward_gz2c = zeros(size(X2)); forward_gz2d = zeros(size(X2));
forward_gz3c = zeros(size(X3)); forward_gz3d = zeros(size(X3));

z4 = 1; z5 = 110;

for j = 1:Nx1
    for i = 1:Ny1
        for k = 1:num_non_zero
            kk = ind_nz(k);
            xm = [x(kk), y(kk), z(kk)]; % coordinates of non-zero
            dm = rho(kk)*dV;
            forward_gz1c(i,j)= forward_gz1c(i,j) + grav_eff_point(([X1(i,j),Y1(i,j),z4]),xm,dm,G);
            forward_gz1d(i,j)= forward_gz1d(i,j) + grav_eff_point(([X1(i,j),Y1(i,j),z5]),xm,dm,G);
        end
    end
end

for j = 1:Nx2
    for i = 1:Ny2
        for k = 1:num_non_zero
            kk = ind_nz(k);
            xm = [x(kk), y(kk), z(kk)]; % coordinates of non-zero
            dm = rho(kk)*dV;
            forward_gz2c(i,j)= forward_gz2c(i,j) + grav_eff_point(([X2(i,j),Y2(i,j),z4]),xm,dm,G);
            forward_gz2d(i,j)= forward_gz2d(i,j) + grav_eff_point(([X2(i,j),Y2(i,j),z5]),xm,dm,G);
        end
    end
end

for j = 1:Nx3
    for i = 1:Ny3
        for k = 1:num_non_zero
            kk = ind_nz(k);
            xm = [x(kk), y(kk), z(kk)]; % coordinates of non-zero
            dm = rho(kk)*dV;
            forward_gz3c(i,j)= forward_gz3c(i,j) + grav_eff_point(([X3(i,j),Y3(i,j),z4]),xm,dm,G);
            forward_gz3d(i,j)= forward_gz3d(i,j) + grav_eff_point(([X3(i,j),Y3(i,j),z5]),xm,dm,G);
        end
    end
end

% calculate dg/dz at z = 0 m
dg_dz1a = (forward_gz1c - forward_gz1a) / (z4 - z1);
dg_dz2a = (forward_gz2c - forward_gz2a) / (z4 - z1);
dg_dz3a = (forward_gz3c - forward_gz3a) / (z4 - z1);

% calculate dg/dz at z = 100 m
dg_dz1b = ((forward_gz1d - forward_gz1b) / (z5 - z3)) / 1e-5;
dg_dz2b = ((forward_gz2d - forward_gz2b) / (z5 - z3)) / 1e-5;
dg_dz3b = ((forward_gz3d - forward_gz3b) / (z5 - z3)) / 1e-5;

% For z = 0 m
figure(); 
sgtitle('Approximation of dg/dz (mGal/m) at z = 0 and 100 m')
subplot(2,3,1)
contourf(X1,Y1,dg_dz1a);
hold on; 
plot(X1,Y1,'xk','MarkerSize',2);
c5 = colorbar; 
colormap(parula);
caxis([min(dg_dz1a,[],'all') 0]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 5 m', 'fontweight','bold');
ylabel(c5, 'dg/dz [mGal/m]', 'fontweight','bold');

subplot(2,3,2)
contourf(X2,Y2,dg_dz2a);
hold on; 
plot(X2,Y2,'xk','MarkerSize',2);
c5 = colorbar; 
caxis([min(dg_dz1a,[],'all') 0]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 10 m', 'fontweight','bold');
ylabel(c5, 'dg/dz [mGal/m]', 'fontweight','bold');

subplot(2,3,3)
contourf(X3,Y3,dg_dz3a);
hold on; 
plot(X3,Y3,'xk','MarkerSize',2);
c5 = colorbar; 
caxis([min(dg_dz1a,[],'all') 0]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 25 m', 'fontweight','bold');
ylabel(c5, 'dg/dz [mGal/m]', 'fontweight','bold');

% For z = 100 m
subplot(2,3,4)
contourf(X1,Y1,dg_dz1b);
hold on; 
plot(X1,Y1,'xk','MarkerSize',2);
c5 = colorbar; 
colormap(parula);
caxis([min(dg_dz1b,[],'all') 0]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 5 m', 'fontweight','bold');
ylabel(c5, 'dg/dz [mGal/m]', 'fontweight','bold');

subplot(2,3,5)
contourf(X2,Y2,dg_dz2b);
hold on; 
plot(X2,Y2,'xk','MarkerSize',2);
c5 = colorbar; 
caxis([min(dg_dz1b,[],'all') 0]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 10 m', 'fontweight','bold');
ylabel(c5, 'dg/dz [mGal/m]', 'fontweight','bold');

subplot(2,3,6)
contourf(X3,Y3,dg_dz3b);
hold on; 
plot(X3,Y3,'xk','MarkerSize',2);
c5 = colorbar; 
caxis([min(dg_dz1b,[],'all') 0]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 25 m', 'fontweight','bold');
ylabel(c5, 'dg/dz [mGal/m]', 'fontweight','bold');

%% Question 5e

% build empty contour plot d2g/dgz2
lapx1a = zeros(length(x_range1)-1,length(x_range1)-1); lapy1a = zeros(length(x_range1)-1,length(x_range1)-1);
lapx2a = zeros(length(x_range2)-1,length(x_range2)-1); lapy2a = zeros(length(x_range2)-1,length(x_range2)-1);
lapx3a = zeros(length(x_range3)-1,length(x_range3)-1); lapy3a = zeros(length(x_range3)-1,length(x_range3)-1);

lapx1b = zeros(length(x_range1)-1,length(x_range1)-1); lapy1b = zeros(length(x_range1)-1,length(x_range1)-1);
lapx2b = zeros(length(x_range2)-1,length(x_range2)-1); lapy2b = zeros(length(x_range2)-1,length(x_range2)-1);
lapx3b = zeros(length(x_range3)-1,length(x_range3)-1); lapy3b = zeros(length(x_range3)-1,length(x_range3)-1);

% calculate d2g/dx2 and d2g/dy2 with dx = 5 m
i = 2;
while i <= length(forward_gz1a) - 1
    j = 2;
    while j <= length(forward_gz1a) - 1
        lapx1a(i,j) = (forward_gz1a(i,j+1) - 2*forward_gz1a(i,j) + forward_gz1a(i,j-1)) / (1e-5 * dx1^2);
        lapy1a(i,j) = (forward_gz1a(i+1,j) - 2*forward_gz1a(i,j) + forward_gz1a(i-1,j)) / (1e-5 * dx1^2);

        lapx1b(i,j) = (forward_gz1b(i,j+1) - 2*forward_gz1b(i,j) + forward_gz1b(i,j-1)) / (1e-5 * dx1^2);
        lapy1b(i,j) = (forward_gz1b(i+1,j) - 2*forward_gz1b(i,j) + forward_gz1b(i-1,j)) / (1e-5 * dx1^2);
        j = j + 1;
    end
    i = i + 1;
end

% calculate d2g/dx2 and d2g/dy2 with dx = 10 m
i = 2;
while i <= length(forward_gz2a) - 1
    j = 2;
    while j <= length(forward_gz2a) - 1
        lapx2a(i,j) = (forward_gz2a(i,j+1) - 2*forward_gz2a(i,j) + forward_gz2a(i,j-1)) / (1e-5 * dx2^2);
        lapy2a(i,j) = (forward_gz2a(i+1,j) - 2*forward_gz2a(i,j) + forward_gz2a(i-1,j)) / (1e-5 * dx2^2);

        lapx2b(i,j) = (forward_gz2b(i,j+1) - 2*forward_gz2b(i,j) + forward_gz2b(i,j-1)) / (1e-5 * dx2^2);
        lapy2b(i,j) = (forward_gz2b(i+1,j) - 2*forward_gz2b(i,j) + forward_gz2b(i-1,j)) / (1e-5 * dx2^2);
        j = j + 1;
    end
    i = i + 1;
end

% calculate d2g/dx2 and d2g/dy2 with dx = 25 m
i = 2;
while i <= length(forward_gz3a) - 1
    j = 2;
    while j <= length(forward_gz3a) - 1
        lapx3a(i,j) = (forward_gz3a(i,j+1) - 2*forward_gz3a(i,j) + forward_gz3a(i,j-1)) / (1e-5 * dx3^2);
        lapy3a(i,j) = (forward_gz3a(i+1,j) - 2*forward_gz3a(i,j) + forward_gz3a(i-1,j)) / (1e-5 * dx3^2);

        lapx3b(i,j) = (forward_gz3b(i,j+1) - 2*forward_gz3b(i,j) + forward_gz3b(i,j-1)) / (1e-5 * dx3^2);
        lapy3b(i,j) = (forward_gz3b(i+1,j) - 2*forward_gz3b(i,j) + forward_gz3b(i-1,j)) / (1e-5 * dx3^2);
        j = j + 1;
    end
    i = i + 1;
end

% calculate d2g/dz2
lap1a = -(lapx1a + lapy1a); lap1b = -(lapx1b + lapy1b);
lap2a = -(lapx2a + lapy2a); lap2b = -(lapx2b + lapy2b);
lap3a = -(lapx3a + lapy3a); lap3b = -(lapx3b + lapy3b);

% For z = 0 m
figure(); 
sgtitle('Approximation of d^2g/dz^2 (mGal/m^2) at z = 0 and 100 m')
subplot(2,3,1)
contourf(X1(2:end-1,2:end-1),Y1(2:end-1,2:end-1),lap1a(2:end,2:end));
hold on; 
plot(X1(2:end-1,2:end-1),Y1(2:end-1,2:end-1),'xk','MarkerSize',2);
c5 = colorbar; 
caxis([0 max(lap1a,[],'all')]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 5 m', 'fontweight','bold');
ylabel(c5, 'd^2g/dz^2 [mGal/m^2]', 'fontweight','bold');

subplot(2,3,2)
contourf(X2(2:end-1,2:end-1),Y2(2:end-1,2:end-1),lap2a(2:end,2:end));
hold on; 
plot(X2(2:end-1,2:end-1),Y2(2:end-1,2:end-1),'xk','MarkerSize',2);
c5 = colorbar; 
caxis([0 max(lap1a,[],'all')]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 10 m', 'fontweight','bold');
ylabel(c5, 'd^2g/dz^2 [mGal/m^2]', 'fontweight','bold');

subplot(2,3,3)
contourf(X3(2:end-1,2:end-1),Y3(2:end-1,2:end-1),lap3a(2:end,2:end));
hold on; 
plot(X3(2:end-1,2:end-1),Y3(2:end-1,2:end-1),'xk','MarkerSize',2);
c5 = colorbar; 
caxis([0 max(lap1a,[],'all')]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 0 m and dx = 25 m', 'fontweight','bold');
ylabel(c5, 'd^2g/dz^2 [mGal/m^2]', 'fontweight','bold');

subplot(2,3,4)
contourf(X1(2:end-1,2:end-1),Y1(2:end-1,2:end-1),lap1b(2:end,2:end));
hold on; 
plot(X1(2:end-1,2:end-1),Y1(2:end-1,2:end-1),'xk','MarkerSize',2);
c5 = colorbar; 
caxis([0 max(lap1b,[],'all')]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 5 m', 'fontweight','bold');
ylabel(c5, 'd^2g/dz^2 [mGal/m^2]', 'fontweight','bold');

subplot(2,3,5)
contourf(X2(2:end-1,2:end-1),Y2(2:end-1,2:end-1),lap2b(2:end,2:end));
hold on; 
plot(X2(2:end-1,2:end-1),Y2(2:end-1,2:end-1),'xk','MarkerSize',2);
c5 = colorbar; 
caxis([0 max(lap1b,[],'all')]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 10 m', 'fontweight','bold');
ylabel(c5, 'd^2g/dz^2 [mGal/m^2]', 'fontweight','bold');

subplot(2,3,6)
contourf(X3(2:end-1,2:end-1),Y3(2:end-1,2:end-1),lap3b(2:end,2:end));
hold on; 
plot(X3(2:end-1,2:end-1),Y3(2:end-1,2:end-1),'xk','MarkerSize',2);
c5 = colorbar; 
caxis([0 max(lap1b,[],'all')]);
axis equal;
xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold');
title('z = 100 m and dx = 25 m', 'fontweight','bold');
ylabel(c5, 'd^2g/dz^2 [mGal/m^2]', 'fontweight','bold');
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