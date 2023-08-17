%% *GOPH 517 Lab 1* 
%% *Fazlie Latib    30067991*
% 
%% 1) Loading the given data file

load('goph_517_lab1_data.mat')
%% 2) Calculating and plotting density, p-wave velocity, impedance and reflection coefficients against depth

z; % depth

I = rho .* vp; % calculate impedance at each interface

i = 2:length(I); % set numbers of layers

R = (I(i) - I(i-1)) ./ (I(i) + I(i-1)); % calculate reflection coefficient at each interface

% plot density against depth
figure();
subplot(1,4,1);
x = rho;
y = z;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('Density (kg/m^3)')
ylabel('Depth (m)')

% plot p-wave velocity against depth
subplot(1,4,2);
x = vp;
y = z;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('P-wave velocity (m/s)')

% plot impedance against depth
subplot(1,4,3);
x = I;
y = z;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('Impedance (kg/m^2 s)')

% plot reflection coefficient against depth
subplot(1,4,4);
x = R;
y = z(2:256);
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('Reflection coefficient')
xlabel('Amplitude')
%% 3) Defining source wavelet

% plot source wavelet, w
figure()
x = 0:0.001:1/32;
w = sin(32*x*2*pi);
plot(w,x)
ax = gca;
ax.YDir = 'reverse';
title('Source wavelet')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 1/32])
%% 4) Performing and comparing convolution using the self-built convolution function and the MATLAB built-in convolution function

conv_517(R,w) % perform convolution using written function between reflection coefficient and source wavelet as example
conv(R,w) % check the previous result with the built-in convolution fucntion
%% 5) Calculating the two-way traveltime to each interface

% use a for loop to calculate two-way travel time at each interface
t = 0;
twt = zeros(1,255);
for i = 2:length(z)
    twt(1,i-1) = t + 2 * (z(1,i) - z(1,i-1)) / vp(1,i);
    t = twt(1,i-1);
end

twt % vector of two-way travel time at each interface
%% 6) Defining r_R as the vector of reflection coefficients corresponding to each two-way travel time and its corresponding synthetic seismogram, s_R

% use a for loop to place each reflection coefficient according to its two-way travel time
n = round(twt(1,255)/0.001);
r_R = zeros(1,n);
for i = 1:length(twt)
    m = round(twt(1,i)/0.001);
    r_R(1,m) = R(i);
end

r_R; % vector for reflection coeffiecient corresponding to its two-way travel time
s_R = conv_517(r_R,w); % convolve r_R with source wavelet

% plot r_R against time
figure();
subplot(1,2,1);
y = 0.001:0.001:length(r_R)*0.001;
x = r_R;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('r_R')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])

% plot synthetic seismogram of r_R * w against time
subplot(1,2,2);
y = 0.001:0.001:length(s_R)*0.001;
x = s_R;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('s_R')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])
%% 7) Defining r_P as the vector of primary reflection amplitudes corresponding to each two-way travel time and its corresponding synthetic seismogram, s_P

% use a for loop to calculate the primary reflection amplitudes at each interface
RA = zeros(1,255);
RA(1,1) = R(1,1);
k = 1;
for i = 2:length(R)
    for ii = 1:i-1
        knew = k * (1 - (R(1,ii))^2); 
        k = knew;
    RA(1,i) = knew * R(1,i);
    end
    k = 1;
end

RA; % vector of primary reflection amplitudes

% use a for loop to place each primary reflection amplitude corresponding to its two-way travel time 
for i = 1:length(twt)
    m = round(twt(1,i)/0.001);
    r_P(1,m) = RA(i);
end

r_P; % vector for primary reflection amplitude corresponding to its two-way travel time
s_P = conv_517(r_P,w); % convolve r_P with source wavelet

% plot r_P against time
figure();
subplot(1,2,1);
y = 0.001:0.001:length(r_P)*0.001;
x = r_P;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('r_P')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])

% plot synthetic seismogram of r_P * w against time
subplot(1,2,2);
y = 0.001:0.001:length(s_P)*0.001;
x = s_P;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('s_P')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])
%% 7) Defining r_M as the vector of primary and first order multiples reflection amplitudes corresponding to each two-way travel time and its corresponding synthetic seismogram, s_M

r_M = zeros(1,length(r_R) .* 2); % set empty vector for r_M
N = length(r_R)-1;

% use nested for loop to calculate the first order multiple reflection amplitudes
for n = 1:N
    trans_coeff_n = 1;
    for i = 1:N
        trans_coeff_n = trans_coeff_n .* (1 - r_R(i));
    end

    for m = n-1:-1:1
        trans_coeff_m = 1;

        for ii = n:m
            trans_coeff_m = trans_coeff_m .* (1 + r_R(ii));
        end

        for p = m+1:N
            trans_coeff_p = 1;
            for iii = p:m
                trans_coeff_p = trans_coeff_p .* (1 - r_R(iii));
            end

          ampl_mult = (trans_coeff_n .* r_R(n)) .* (trans_coeff_m .* -r_R(m)) .* (trans_coeff_p .* r_R(p));
          sampl_ampl_n = ((n+(n-m)) + (p-m) + p) ./ 2;
          r_M(sampl_ampl_n) = r_M(sampl_ampl_n) + ampl_mult;

        end
    end
end

r_M = r_M(1:2030); % set r_M until the specified time
r_M = r_M + r_P; % add primary amplitudes to the first order multiple reflection amplitudes
s_M = conv_517(r_M,w); % convolve r_M with source wavelet

% plot r_M against time
figure();
subplot(1,2,1);
y = 0.001:0.001:length(r_M)*0.001;
x = r_M;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('r_M')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])

% plot synthetic seismogram of r_M * w against time
subplot(1,2,2);
y = 0.001:0.001:length(s_M)*0.001;
x = s_M;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('s_M')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])
%% 8) Calculating attenuation using quality factor, Q, and their plots

% perform interpolation to transform Q corresponding to its two way travel time
Q_new = Q(1:255);
timesample = 0:0.001:twt(1,255);
Qt = interp1(twt,Q_new,timesample);

% calculate attenuation
att = zeros(1,length(Qt));
for i = 1:length(Qt)
    att(1,i) = exp(-2*pi*32*timesample(i)/(2*Qt(i)));
end

% plot Q against depth
figure();
subplot(1,3,1);
y = z;
x = Q;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('Q')
xlabel('Quality factor')
ylabel('Depth (m)')

% plot Q against time
subplot(1,3,2);
y = timesample;
x = Qt;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('Q in time')
xlabel('Quality factor')
ylabel('Time (s)')
ylim([0 twt(1,255)])

% plot attenuation against time
subplot(1,3,3);
y = timesample;
x = att;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('Attenuation')
xlabel('Attenuation')
ylabel('Time (s)')
ylim([0 twt(1,255)])
%% 9) Defining r_Q as the vector of attenuated primary and first order multiples reflection amplitudes corresponding to each two-way travel time and its corresponding synthetic seismogram, s_Q

r_Q = r_M .* (att); % calculate new amplitude due to attenuation
r_Q = fillmissing(r_Q,"constant",0);

s_Q = conv_517(r_Q,w); % convolve r_M with source wavelet

% plot r_Q against time
figure();
subplot(1,2,1);
y = 0.001:0.001:length(r_Q)*0.001;
x = r_Q;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('r_Q')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])

% plot synthetic seismogram of r_Q * w against time
subplot(1,2,2);
y = 0.001:0.001:length(s_Q)*0.001;
x = s_Q;
plot(x,y)
ax = gca;
ax.YDir = 'reverse';
title('s_Q')
xlabel('Amplitude')
ylabel('Time (s)')
ylim([0 twt(1,255)])
%% Defining a convolution function

function [w]=conv_517(u,v)
    m = length(u); % define length of first input
    n = length(v); % define length of second input
    X = [u,zeros(1,n)]; 
    H = [v,zeros(1,m)]; 
    for i = 1:n+m-1
        w(i) = 0;
    for j = 1:m
        if(i-j+1>0)
            w(i) = w(i)+X(j)*H(i-j+1); % perform convolution
        end
    end
    end
end