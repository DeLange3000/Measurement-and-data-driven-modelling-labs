 clc
 close all
 clear

 %% generating signals

 % random phase
 fs = 16e3; %Hz
 min_freq = 4; %Hz
 max_freq = 1000; %Hz
 Vrms = 0.5; %V
 N1 = 4000;
 T = N1/fs; %measurment period
 bin_size = 1/T;

 x_random_f = zeros(N1, 1);
 x_random_f((min_freq/bin_size):(max_freq/bin_size)) = 1*exp(1i*2*pi*rand(size((min_freq/4):(max_freq/4))));

 x_random_t = real(ifft(x_random_f));
 x_random_t = x_random_t/rms(x_random_t)*Vrms; %get vrms of 0.5V

 figure
 plot(1:N1, abs(x_random_f))
 xlabel('bins')
 ylabel('amplitude [V]')
 title('amplitude of multisine in frequency domain')

 figure
 plot((1:N1)/fs, x_random_t)
 xlabel('time [s]')
 ylabel('amplitude [V]')
 title('multisine in time domain')


 % aperiodic noise
 fs = 16e3; %Hz
 Vrms = 0.5; %V
 N2 = 40*N1;
 x_noise_t = randn(1, N2);
 x_noise_t = x_noise_t/rms(x_noise_t)*Vrms;

 figure
 plot((1:N2)/fs, x_noise_t)
 xlabel('time [s]')
 ylabel('amplitude [V]')
 title('noise in time domain')

 
%% data processing

interesting_bins = 1:250; %(min_freq/bin_size):(max_freq/bin_size); %bins we are interested in
N = 160000; % amount of samples in 1 repetition
Nrep = 1; % how many repetitions were measured
Drep = 1; % how many of the last repetitions you want to keep
FileName = 'matlab data/x_noise_data_1.mat'; %name of the mat file
Avgs = 16; %average over how many repetitions
input_noise = true; % set to true if input is one long measurement (noise input)
Hfunction = 'AveragingFRF'; %name of function to call
% 'AvgTimeDomain'
% 'AveragingDFT'
% 'AveragingFRF'
% 'AveragingAutoPowerInput'
% 'AveragingAutoPowerOutput'


[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);


% if noise is used as input
if (input_noise)
    tempu = zeros(N1, N2/N1);
    tempy = zeros(N1, N2/N1);
    j = 1;
    for i = 1:N1:N2
        tempu(:, j) = umat(i:i+N1-1);
        tempy(:, j) = ymat(i:i+N1-1);
        j = j + 1;
    end
    umat = tempu;
    ymat = tempy;
end


[H, stdH] = TransferFunc(umat, ymat, Avgs, Hfunction);

[a, b] = size(H);
for i = 1:b
    figure
    plot(interesting_bins, db(abs(H(interesting_bins,i))))
     xlabel('bins')
     ylabel('amplitude [V]')
     title('FRF with averaging '+ string(Hfunction) + ' using '+ string(Avgs) + ' repetitions')
     hold on
     plot(interesting_bins, db(stdH(interesting_bins)))
     legend('FRF', 'stdH')
end

