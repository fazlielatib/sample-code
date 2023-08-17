% GOPH 547 : Lab 2 - Part 1
% Fazlie Latib
% 30067991

load goph_547_lab_2_data.mat

% define variables in grav_survey_data
base_point  = grav_survey_data(:,1);
day         = grav_survey_data(:,2);
time        = grav_survey_data(:,3);
total_time  = grav_survey_data(:,4);
x           = grav_survey_data(:,5);
y           = grav_survey_data(:,6);
z           = grav_survey_data(:,7);
g           = grav_survey_data(:,8);
dg_tide     = grav_survey_data(:,9);

%%% Data preparation

%% Q1 - Plot time variation of relative gravity measurements and the tidal variation
figure();
g(1) = 0;
p = polyfit(total_time,g,1);
f = polyval(p,total_time);
plot(total_time,g,'xk','MarkerSize',0.5,'LineWidth',1.5);
hold on
plot(total_time,dg_tide,'o','MarkerSize',0.7);
hold on
plot(total_time,f,'--','MarkerSize',0.7)
xlabel('Time [hour]');
ylabel('Relative gravity reading [\muGal]');
legend('raw gravity measurements','tidal variation','mean drift','Location','northwest');
title('Relative gravity measurement and tidal variations over time');

%% Q2 - Create a contour plot of the terrain survey data
figure();
contourf(Xt,Yt,Zt);
axis equal;
grid on;
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Elevation [m]','fontweight','bold');
xlabel('Eastings [km]'); ylabel('Northings [km]');
title('Elevation data');

% display elevation information
disp(string(strcat('Mean elevation:',{' '},mat2str(mean(min(Zt))))));
disp(string(strcat('Range of elevations:',{' '},mat2str(min(min(Zt))),{' '},'to',{' '},mat2str(max(max(Zt))))));

z_datum = min(min(Zt)); % set the datum

%% Q3 - Create a new three-column matrix x_sort with the x, y and z data in the columns
x_sort = [x y z];
[x_sorted,ind_sort] = sortrows(x_sort,[1 2]); % sort by x-coordinate then y-coordinate

%% Q4 - Extract the unique entries from x_sort
[uniq_x,ind_uniq,~] = unique(x_sorted,'rows');

%% Q5 - Save columns in x_sort as three separate column vectors
Xg = uniq_x(:,1);
Yg = uniq_x(:,2);
Zg = uniq_x(:,3);

Nx = sqrt(length(Xg));
Ny = sqrt(length(Yg));

%% Q6

% convert each of Xg, Yg and Zg into a grid
Xg = reshape(Xg,[Nx,Ny]);
Yg = reshape(Yg,[Nx,Ny]);
Zg = reshape(Zg,[Nx,Ny]);

% create the grids with meshgrid function 
Xgrid = 0:50;
Ygrid = 0:50;
[Xg2,Yg2] = meshgrid(Xgrid,Ygrid);

% verify step
x_diff = Xg - Xg2;
y_diff = Yg - Yg2;
verify_x = find(x_diff);
verify_y = find(y_diff);

%% Q7

g(1) = grav_survey_data(1,8);
g_corr = g;
g_raw = g_corr;
g_raw(2:end) = g_raw(2:end) + g_raw(1);
g_raw = g_raw(ind_sort);
g_raw = g_raw(ind_uniq);
g_raw = reshape(g_raw,[Nx,Ny]);

figure();
contourf(Xg,Yg,g_raw);
title('Raw');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%%% Normal gravity correction

%% Q8

ge = 9.780327; A = 0.0053024; B = 0.0000058;
gt = ge * (1 + A * ((sind(49.1286))^2) - B * ((sind(2 * 49.1286))^2)); % use IGF to compute gt
g_corr(1) = g_corr(1) - gt * 10^8;

%% Q9

g_norm = g_corr;
g_norm(2:end) = g_norm(2:end) + g_norm(1);
g_norm = g_norm(ind_sort);
g_norm = g_norm(ind_uniq);
g_norm = reshape(g_norm,[Nx,Ny]);

figure();
contourf(Xg,Yg,g_norm);
title('Raw - Norm');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%%% Drift corrections

%% Q10
g_corr = g_corr - dg_tide;

%% Q11
g_tide = g_corr;
g_tide(2:end) = g_tide(2:end) + g_tide(1);
g_tide = g_tide(ind_sort);
g_tide = g_tide(ind_uniq);
g_tide = reshape(g_tide,[Nx,Ny]);

figure();
contourf(Xg,Yg,g_tide);
title('Raw - Norm - Tide');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

%% Q12

st = 1;
for i = 2:length(base_point)  
    if base_point(i) == 1
        st = i;
        drift = g_corr(i); 
        g_corr(st:end) = g_corr(st:end) - drift; % subtract overnight drift
    end 
    
    if base_point(i) == 2
        drift = g_corr(i);       
        drift_rate = drift / (total_time(i) - total_time(st));
        g_corr(st+1:i-1) = g_corr(st+1:i-1) - drift_rate * (total_time(st+1:i-1) - total_time(st)); % subtract the time drift 
        g_corr(i:end) = g_corr(i:end) - drift;  % subtract the drift from points that are not the start nor end of loop
    end
end

%% Q13

% verify g_corr = 0
for i=2:length(g_corr)
    if (base_point(i)==1||base_point(i)==2)
        if(g_corr(i)~=0)
            disp('error');
        end
    end
end

% add the value of g_corr(1) to all other points in g_corr
g_corr(2:end)=g_corr(2:end) + g_corr(1);

% plot g_corr against total_time
figure();
plot(total_time,g_corr,'xk','MarkerSize',0.5,'LineWidth',1.5);
xlabel('Time [hour]');
ylabel('Gravity effect [\muGal]');
title('Drift corrected gravity time series')

%% Q14

g_corr = g_corr(ind_sort);
g_corr = g_corr(ind_uniq);
g_corr = reshape(g_corr,[Nx,Ny]);
g_drft = g_corr;

figure();
contourf(Xg,Yg,g_drft);
title('Raw - Norm - Drift');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

% 
figure();
subplot(2,2,1)
contourf(Xg,Yg,g_raw);
axis equal;
title('Raw');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(2,2,2)
contourf(Xg,Yg,g_norm);
axis equal;
title('Raw - Norm');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(2,2,3)
contourf(Xg,Yg,g_tide);
axis equal;
title('Raw - Norm - Tide');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');

subplot(2,2,4)
contourf(Xg,Yg,g_drft);
axis equal;
title('Raw - Norm - Drift');
xlabel('Easting [km]');
ylabel('Northing [km]');
grid on
c1 = colorbar;
colormap(flipud(parula));
ylabel(c1,'Gravity effect [\muGal]','fontweight','bold');