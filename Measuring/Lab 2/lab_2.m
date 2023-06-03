clc
close all
clear

%% Question 2.3 - 2.5
% multisine

amount_of_samples = 4096;
K = 100;
f_constant = zeros(1, amount_of_samples);
f_constant(1:K) = 1;

figure
stem(f_constant);
xlabel('bins')
ylabel('amplitude')
title('frequency spectrum of multisine with constant phase')

% constant phase spectrum
x_constant = amount_of_samples*real(ifft(f_constant));

figure
plot(x_constant)
xlabel('samples')
ylabel('amplitude')
title('time domain of multisine with constant phase')

% random phase spectrum
f_random = f_constant.*exp(1i*2*pi*rand(size(f_constant)));
x_random = amount_of_samples*real(ifft(f_random));

figure
plot(x_random)
xlabel('samples')
ylabel('amplitude')
title('time domain of multisine with random phase')

%schroeder phase spectrum
phase = ((1:length(f_constant)).*(1 + (1:length(f_constant)))*pi)/K;
f_schroeder = f_constant.*exp(1i*phase);
x_schroeder = amount_of_samples*real(ifft(f_schroeder));

figure
plot(x_schroeder)
xlabel('samples')
ylabel('amplitude')
title('time domain of multisine with schroeder phase')

% crest factors
crest_constant = max(x_constant)/rms(x_constant);
crest_random = max(x_random)/rms(x_random);
crest_schroeder = max(x_schroeder)/rms(x_schroeder);

%crest for constant is high as most of the signal is low due to
%interference from the multisines. -> Low power in signal -> low snr

%crest for random is lower since the multisines do not interfere that much
%with eachother

% crest for schroeder is the lowest since the phase is specially chosen for
% a lower crest factor -> More power, highest snr


% constant plot shows peaks at the beginning and and of the sampling
% period. This is due to interference of the multisines inbetween those
% peaks

% random plot has less interference so the amplitude of the signal is more
% constant. but there are still some variations

% schroeder plot has the most constant amplitude in time domain.

%% lab signals

% schroeder
fs = 8e3; %hz
ts = 1/fs;
K = 500;
max_freq = 500;
freq_resolution = 1; % bin (and in this case Hz)
T = 1/freq_resolution; %measurement time
amount_of_samples = T/ts;
excited_freqs = 1:K;
f_schroeder = zeros(1, amount_of_samples);
f_schroeder(excited_freqs) = 1;
phase = ((1:length(f_schroeder)).*(1 + (1:length(f_schroeder)))*pi)/K;
f_schroeder = f_schroeder.*exp(1i*phase);

x_schroeder = amount_of_samples*real(ifft(f_schroeder));
x_schroeder = x_schroeder/rms(x_schroeder)*0.1; %get Vrms of 0.1 V

figure
stem((1:length(f_schroeder)), abs(f_schroeder))
xlabel('bins')
ylabel('amplitude')
title('[LAB] frequency domain of multisine with schroeder phase')

figure
plot((1:length(x_schroeder))*ts, x_schroeder)
xlabel('time [s]')
ylabel('amplitude')
title('[LAB] time domain of multisine with schroeder phase')


% constant phase
K = 500;
f_constant = zeros(1, amount_of_samples);
f_constant(1:K) = 1;

x_constant = amount_of_samples*real(ifft(f_constant));
x_constant = x_constant/rms(x_constant)*0.1; %get Vrms of 0.1 V

figure
stem((1:length(f_constant)), f_constant)
xlabel('bins')
ylabel('amplitude')
title('[LAB] frequency domain of multisine with constant phase')

figure
plot((1:length(x_constant))*ts, x_constant)
xlabel('time [s]')
ylabel('amplitude')
title('[LAB] time domain of multisine with constant phase')


% random phase
K = 500;
f_random = zeros(1, amount_of_samples);
f_random(1:500) = 1;
f_random = f_random.*exp(1i*2*pi*rand(size(f_random)));

x_random = amount_of_samples*real(ifft(f_random));
x_random = x_random/rms(x_random)*0.1; %get Vrms of 0.1 V

figure
stem((1:length(f_random)), abs(f_random))
xlabel('bins')
ylabel('amplitude')
title('[LAB] frequency domain of multisine with random phase')

figure
plot((1:length(x_random))*ts, x_random)
xlabel('time [s]')
ylabel('amplitude')
title('[LAB] time domain of multisine with random phase')


% periodic noise
x_per_noise = rand(1, amount_of_samples);
x_per_noise = x_per_noise/rms(x_per_noise)*0.1;

figure
plot((1:length(x_per_noise))*ts, x_per_noise)
xlabel('time [s]')
ylabel('amplitude')
title('[LAB] time domain of periodic noise')


% aperiodic 
periods = 10;
x_aper_noise = rand(1, amount_of_samples*periods);
x_aper_noise = x_aper_noise/rms(x_aper_noise)*0.1;

figure
plot((1:length(x_aper_noise))*ts, x_aper_noise)
xlabel('time [s]')
ylabel('amplitude')
title('[LAB] time domain of aperiodic noise')

hann = hanning(amount_of_samples*periods, 'periodic'); % creates window

%% Lab
close all; 

%% Schroeder phase

interesting_freqs = 1:500; %bins we are interested in
N = amount_of_samples; % amount of samples in 1 repetition
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'schroeder_data_nog_nekeer_nekeer.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

PlotData(umat, ymat, interesting_freqs, true, true, true, "Schroeder phase")
%% Constant phase

interesting_freqs = 1:500; %bins we are interested in
N = amount_of_samples; % amount of samples in 1 repetition
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'constant_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

PlotData(umat, ymat, interesting_freqs, false, false, true, "Constant phase")

%% Random phase

interesting_freqs = 1:500; %bins we are interested in
N = amount_of_samples; % amount of samples in 1 repetition
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'random_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

PlotData(umat, ymat, interesting_freqs, false, false, true, "Random phase")

%% Periodic noise

interesting_freqs = 1:500; %bins we are interested in
N = amount_of_samples; % amount of samples in 1 repetition
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'random_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

PlotData(umat, ymat, interesting_freqs, false, false, true, "Periodic noise")

%% Aperiodic noise

interesting_freqs = 1:500; %bins we are interested in
N = amount_of_samples*periods; % amount of samples in 1 repetition
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'random_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

PlotData(umat, ymat, interesting_freqs, false, false, true, "Aperiodic noise")

%% Aperiodic noise hann

interesting_freqs = 1:500; %bins we are interested in
N = amount_of_samples; % amount of samples in 1 repetition
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'random_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

PlotData(umat.*hann, ymat.*hann, interesting_freqs, false, false, true, "Aperiodic noise with hann")

