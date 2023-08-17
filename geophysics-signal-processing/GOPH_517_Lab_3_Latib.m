%% *GOPH 517 Lab 3* 
%% *Fazlie Latib    30067991*
% 
%% 1) Loading the given data file

load('goph_517_lab_3_data1.mat')
load('goph_517_lab_3_data2.mat')
load('goph_517_lab_3_wavelet.mat')
%% 2) Calculating minimum phase from the wavelet amplitude spectrum

dt = tw(2) - tw(1);
fnyq = 0.5/dt;
df = fnyq / length(w);
f_wavelet = -fnyq:2*df:fnyq-2*df;

ampl_f_wavelet = abs(fftshift(fft(w)));
[min_phase_t, time_min_phase] = hilbert_trans(ampl_f_wavelet,f_wavelet);

figure();
plot(tw(1:50),w(1:50))
xlabel('Time(s)')
ylabel('Amplitude')
title('Original Wavelet vs Minimum-Phase Wavelet')
hold on
plot(time_min_phase(1:50),min_phase_t(1:50))
legend('Original wavelet','Minimum-phase wavelet','Location','northeast')
%% 3) Finding wavelet autocorrelation for noise-free data

autocorr_w1 = xcorr(w,w);

t_autocorr_w1 = -length(autocorr_w1) * dt / 2 : dt : length(autocorr_w1) * dt / 2 - dt;

figure();
plot(t_autocorr_w1,autocorr_w1)
title('Noise-free data : Wavelet autocorrelation')
xlabel('Time lag (s)')
ylabel('Amplitude')

t_autocorr_w1_new = t_autocorr_w1(1971:2030);
autocorr_w1_new = autocorr_w1(1971:2030);

figure();
plot(t_autocorr_w1_new,autocorr_w1_new)
title('Noise-free data : Wavelet autocorrelation (certain time lags)')
xlabel('Time lag (s)')
ylabel('Amplitude')

[ampl_max1,x_max1] = max(autocorr_w1);

autocorr_w1 = autocorr_w1(x_max1:end);
%% 4) Creating noisy data autocorrelation and finding wavelet autocorrelation and wavelet

autocorr2 = xcorr(data2,data2);

[ampl_max2,x_max2] = max(autocorr2);

figure();
plot(t_autocorr_w1,autocorr2)
title('Noisy data autocorrelation')
xlabel('Time lag (s)')
ylabel('Amplitude')

index_left = x_max2 - 30;
index_right = x_max2 + 30;

window_length = index_right - index_left + 1;

Ns = length(autocorr2);

win_hann = zeros(1,Ns);
win_hann(index_left:index_right) = window(@hann,window_length);
win_hann = (win_hann).';
autocorr_hann2 = autocorr2 .* win_hann;

figure();
plot(t_autocorr_w1,win_hann)
title('Hann Windowing')
xlabel('Time lag (s)')
ylabel('Amplitude')

figure();
plot(t_autocorr_w1,autocorr_hann2)
title('Applying Hann windowing to the noisy data autocorrelation')
xlabel('Time lag (s)')
ylabel('Amplitude')

figure();
% t = -Ns/2 * dt : dt : (Ns/2)-1 * dt;
plot(t_autocorr_w1,autocorr_hann2)
xlim([-30*dt 30*dt])
title('Noisy data : Wavelet autocorrelation')
xlabel('Time lag (s)')
ylabel('Amplitude')

[ampl_max_w2,x_max_w2] = max(autocorr_hann2);

autocorr_w2 = autocorr_hann2(x_max_w2:end);

ampl_f_wavelet2 = abs(fftshift(fft(autocorr_hann2)));

wavelet2 = hilbert_trans(ampl_f_wavelet2,f_wavelet);
tw2 = 0 : dt : (length(wavelet2) - 1) * dt;

[ampl_max_wavelet2,x_max_wavelet2] = max(wavelet2);

figure();
plot(tw2,wavelet2 / ampl_max_wavelet2)
xlim([0 50*dt])
title('Noisy data : Wavelet')
xlabel('Time (s)')
ylabel('Amplitude')

wavelet2_new = wavelet2(1:2000)/ampl_max_wavelet2;
%% 5) Applying zero-lag Wiener deconvolution to data

stab1 = 0;
stab2 = 5e-1;

[wiener_deconv1] = zero_lag_wiener(autocorr_w1,stab1,data1);
[wiener_deconv2] = zero_lag_wiener(autocorr_w2,stab2,data2);

figure();
plot(t1,data1)
hold on
plot(t1,wiener_deconv1(1:2000)/100)
title('Zero-lag Wiener deconvolution on the noise-free data')
xlabel('Time(s)')
ylabel('Amplitude')
legend('Noise-free data','Deconvolved data','Location','northeast')

figure();
plot(t2,data2)
hold on
plot(t2,wiener_deconv2(1:2000)*2.5)
title('Zero-lag Wiener deconvolution on the noisy data')
xlabel('Time(s)')
ylabel('Amplitude')
legend('Noisy data','Deconvolved data','Location','northeast')
%% 6) Applying time-lag deconvolution to data

l1 = 0;
lag1 = 1;

[lag_wiener1] = lag_wiener(data1,w,l1,lag1);

figure();
plot(t1,data1)
hold on
plot(t1,lag_wiener1(1:2000))
title('Time-lag Wiener deconvolution on the noise-free data')
xlabel('Time(s)')
ylabel('Amplitude')
legend('Noise-free data','Deconvolved data','Location','northeast')

l2 = 5e-1;
lag2 = 8;

[lag_wiener2] = lag_wiener(data2,wavelet2_new,l2,lag2);

figure();
plot(t2,data2)
hold on
plot(t2,lag_wiener2(1:2000)*8)
title('Time-lag Wiener deconvolution on the noisy data')
xlabel('Time(s)')
ylabel('Amplitude')
legend('Noisy data','Deconvolved data','Location','northeast')

figure();
subplot(3,1,1)
plot(t1,data1)
xlabel('Time(s)')
ylabel('Amplitude')
title('Original noise-free data')

subplot(3,1,2)
plot(t1,wiener_deconv1(1:2000)/100)
xlabel('Time(s)')
ylabel('Amplitude')
title('Zero-lag Wiener deconvolution on noise-free data')

subplot(3,1,3)
plot(t1,lag_wiener1(1:2000))
xlabel('Time(s)')
ylabel('Amplitude')
title('Time-lag Wiener deconvolution on noise-free data')

figure();
subplot(3,1,1)
plot(t2,data2)
xlabel('Time(s)')
ylabel('Amplitude')
title('Original noisy data')

subplot(3,1,2)
plot(t2,wiener_deconv2(1:2000)*2.5)
xlabel('Time(s)')
ylabel('Amplitude')
title('Zero-lag Wiener deconvolution on noisy data')

subplot(3,1,3)
plot(t2,lag_wiener2(1:2000)*8)
xlabel('Time(s)')
ylabel('Amplitude')
title('Time-lag Wiener deconvolution on noisy data')
%% Defining a minimum phase time domain function

function [min_phase_t, time_min_phase] = hilbert_trans(ampl,freq)
    phase_f = imag(hilbert(log(ampl)));
    min_phase_f = ampl .* exp(phase_f * 1i);
    min_phase_t = flip(real(ifft(ifftshift(min_phase_f))));
    df = 0.5 * abs( freq(1) - freq(2) );
    fnyq = df * length(min_phase_t);
    dt = 0.5/fnyq;
    time_min_phase = 0:dt:(length(min_phase_t) - 1)*dt;
end
%% Defining a zero-lag Wiener deconvolution function

function [deconv] = zero_lag_wiener(autocorr,stability,data)
    autocorr(1) = autocorr(1) + stability;
    autocorr_matrix = toeplitz(autocorr);
    s = zeros(length(autocorr),1);
    s(1) = 1;
    filter = autocorr_matrix \ s;
    deconv = conv(data,filter);
end
%% Defining a time-lag Wiener deconvolution function

function [lag_wiener_deconv] = lag_wiener(dat,wavelet,lamda,lag)
    y = (zeros([length(dat) 1]));
    y(lag) = 1;
    W = tril(toeplitz(wavelet));
    phi = W' * W;
    I = eye(length(phi));
    filter = (phi + (lamda .* I)) \ (W' * y);
    lag_wiener_deconv = conv(filter,dat);
end
%% Defining a frequency domain deconvolution function

% function [R_time_dom] = freq_deconv(data,mu,)
%     R_freq_dom = (fft(data2)) .* conj(w) ./ ((abs(w).^2)+(mu * A_max)^2);
%     R_time_dom = ifft(R_freq_dom);
%     R_time_dom = real(R_time_dom);
%