%% *GOPH 517 Lab 2* 
%% *Fazlie Latib    30067991*

clear
%% Loading the given data file

load('goph_517_lab_2_data.mat')
%% Plotting the original data and the input wavelet

figure();
plot(t,data)
xlabel('Time (s)')
ylabel('Amplitude')
title('Original Data')
plot(tw,wavelet)
xlabel('Time (s)')
ylabel('Amplitude')
title('Source Wavelet')
%% Convolving provided wavelet with a delta function time series of length 500

timeseries = zeros(500,1);
timeseries(50,1) = 0.5;
timeseries(130,1) = -0.35;
timeseries(250,1) = 0.7;
timeseries(350,1) = 0.4;
timeseries(465,1) = -0.2;
conv_1 = conv_freq_mult(timeseries,wavelet);

figure();
subplot(2,1,1)
plot(timeseries)
xlabel('Length')
ylabel('Amplitude')
title('Delta function time series of length 500')

subplot(2,1,2)
plot(conv_1(65:565))
xlabel('Length')
ylabel('Amplitude')
title('Convolved delta function time series with source wavelet')
xlim([0 500])
%% Converting the original data into frequency domain and the corresponding plot

dt = t(2) - t(1);
fnyq = 0.5/dt;
df = fnyq / length(data);
f_data = -fnyq:2*df:fnyq-2*df;

figure();
plot(f_data,abs(fftshift(fft(data))))
title('Original data in frequency domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
%% Building Ormsby filter

% Ormsby
f_lc = 0;
f_hp = 10;
f_lp = 70;
f_hc = 90;

f_end = (length(data) - 1) * df;
f = 0:df:f_end;

A = zeros(length(f),1);

n = 1;
for f_1 = 0:df:f_lc
    ampl = 0;
    A(n,1) = ampl;
    n = n + 1;
end

n = n - 1;
for f_2 = f_lc:df:f_hp
    ampl = (f_2 - f_lc) / (f_hp - f_lc);
    A(n,1) = ampl;
    n = n + 1;
end

n = n - 1; 
for f_3 = f_hp:df:f_lp
    ampl = 1;
    A(n,1) = ampl;
    n = n + 1;
end

n = n - 1;
for f_4 = f_lp:df:f_hc
    ampl = (f_hc - f_4) / (f_hc - f_lp);
    A(n,1) = ampl;
    n = n + 1;
end

n = n - 1;
for f_5 = f_hc:df:f_end
    ampl = 0;
    A(n,1) = ampl;
    n = n + 1;
end

figure();
plot(f,A)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Ormsby filter')
%% Applying Ormsby filter to the original data and the corresponding plot

ormsby_f = fft(data) .* A;
ormsby_t = ifft(ormsby_f);

figure();
plot(t,data)
hold on
plot(t,real(ormsby_t))
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original data','Filtered data (Ormsby filter)','Location','northwest')
%% Building Butterworth filter

% Butterworth
f_h = 5;
n_h = 10;
f_l = 70;
n_l = 15;

B = zeros(length(data),1);

n = 1;
for f_b1 = 0:df:40
    w_h = f_b1 ./ f_h;
    B_h = sqrt(w_h .^ (2 * n_h) ./ (1 + w_h .^ (2 * n_h)));
    B(n,1) = B_h;
    n = n + 1;
end

n = n - 1;
for f_b2 = 41:df:1000
    w_l = f_b2 ./ f_l;
    B_l = sqrt(1 ./ (1 + w_l .^ (2 * n_l)));
    B(n,1) = B_l;
    n = n + 1;
end

figure();
plot(f,B)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Butterworth filter')
%% Applying Butterworth filter to the original data and the corresponding plot

butter_f = fft(data) .* B;
butter_t = ifft(butter_f);

figure();
plot(t,data)
hold on
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original data','Filtered data (Butterworth filter)','Location','northwest')
%% Comparing band-passed filtered results with original data

% Compare seismic trace with and without filter
figure();
subplot(3,1,1)
plot(t,data)
xlabel('Time (s)')
ylabel('Amplitude')
title('Original data')
xlim([0 3])

subplot(3,1,2)
plot(t,real(ormsby_t))
xlabel('Time (s)')
ylabel('Amplitude')
title('Filtered data using Orsmby filter')
xlim([0 3])

subplot(3,1,3)
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
title('Filtered data using Butterworth filter')
xlim([0 3])
%% Building short and long prediction filter

% Prediction filter
y = data(2:(length(data)/2));

length_pf1 = 30;
length_pf2 = 2000;

S = tril(toeplitz(data(1:(length(data)/2)-1)));

S_pf1 = S(:,(1:length_pf1));
S_pf2 = S(:,(1:length_pf2));

pf1 = (inv(S_pf1.' * S_pf1)) * (S_pf1.' * y);
pf2 = (inv(S_pf2.' * S_pf2)) * (S_pf2.' * y);
%% Applying prediction filter to original data and band-passed data, and its corresponding plot

data_f = tril(toeplitz(data));
bandpass_f = tril(toeplitz(real(butter_t)));

data_pf1 = data_f(:,(1:length(pf1)));
data_pf2 = data_f(:,(1:length(pf2)));

bandpass_pf1 = data_f(:,(1:length(pf1)));
bandpass_pf2 = data_f(:,(1:length(pf2)));

x_pf1 = data_pf1 * pf1;
x_pf2 = data_pf2 * pf2;

xbandpass_pf1 = bandpass_pf1 * pf1;
xbandpass_pf2 = bandpass_pf2 * pf2;

figure();
plot(t,data)
hold on
plot(t,x_pf1)
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original data','Filtered data (Short prediction filter)','Location','northwest')
figure();
plot(t,data)
hold on
plot(t,x_pf2)
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original data','Filtered data (Long prediction filter)','Location','northwest')
figure();
plot(t,xbandpass_pf1)
hold on
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
legend('Filtered data (Short prediction filter)','Band-passed data','Location','northwest')
figure();
plot(t,xbandpass_pf2)
hold on
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
legend('Filtered data (Long prediction filter)','Band-passed data','Location','northwest')
% Comparing short- and long-filtered original data results with original data

% Compare seismic trace with short and length filter
figure();
subplot(3,1,1)
plot(t,data)
xlabel('Time (s)')
ylabel('Amplitude')
title('Original data')
xlim([0 3])

subplot(3,1,2)
plot(t,x_pf1)
xlabel('Time (s)')
ylabel('Amplitude')
title('Short-filtered original data')
xlim([0 3])

subplot(3,1,3)
plot(t,x_pf2)
xlabel('Time (s)')
ylabel('Amplitude')
title('Long-filtered original data')
xlim([0 3])
%% Comparing short- and long-filtered band-passed data results with band-passed data

figure();
subplot(3,1,1)
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
title('Band-passed data')
xlim([0 3])

subplot(3,1,2)
plot(t,xbandpass_pf1)
xlabel('Time (s)')
ylabel('Amplitude')
title('Short-filtered band-passed data')
xlim([0 3])

subplot(3,1,3)
plot(t,xbandpass_pf2)
xlabel('Time (s)')
ylabel('Amplitude')
title('Long-filtered band-passed data')
xlim([0 3])
%% Building a matched filter and applying it to original data

mf = flip(wavelet);
x_mf = conv_freq_mult(data,mf);
%% Comparing original data, band-passed data and matched filtered original data

figure();
subplot(3,1,1)
plot(t,data)
xlabel('Time (s)')
ylabel('Amplitude')
title('Original data')
xlim([0 3])

subplot(3,1,2)
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
title('Band-passed data')
xlim([0 3])

subplot(3,1,3)
plot(t,x_mf(1:length(t)))
xlabel('Time (s)')
ylabel('Amplitude')
title('Matched filtered original data')
xlim([0 3])
%% Comparing band-passed data and matched filtered band-passed data

bandpass_mf = conv_freq_mult(real(butter_t),mf);

figure();
subplot(3,1,1)
plot(t,real(butter_t))
xlabel('Time (s)')
ylabel('Amplitude')
title('Band-passed data')
xlim([0 3])

subplot(3,1,2)
plot(t,bandpass_mf(1:length(t)))
xlabel('Time (s)')
ylabel('Amplitude')
title('Matched filtered band-passed data')
xlim([0 3])
%% Convolution function through multiplication of data in frequency domain

function [c] = conv_freq_mult(a,b)
    a_fft = fftshift(fft(a, length(a) + length(b) - 1)); 
    b_fft = fftshift(fft(b, length(a) + length(b) - 1)); 
    c = ifft(ifftshift(a_fft .* b_fft)); 
end