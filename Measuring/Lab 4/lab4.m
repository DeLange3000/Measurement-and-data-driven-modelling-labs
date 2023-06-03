clc
close all
clear

%% generating nonlinear response

T = 44;
t = 0:0.001:T - 0.001;
x = sin(2*pi*4*t) + sin(2*pi*11*t);

figure
subplot(2,1,1)
plot(x)

y = x- 1/2*x.^2 - 1/4*x.^4;

subplot(2,1,2)
plot(y)

figure
plot(0:1/T:length(t)/T - 1/T ,abs(fft(y)))

%% lab signals

%% 4.2
fs = 10e3; %Hz
N = 4096;
T = N*1/fs;
bin = 1/T; % freq resolution Hz

f = 102.5391; %Hz f0 = fs/N
t = 0:1/fs:(N - 1)/fs;

x_4_2 = sin(2*pi*f*t);

figure
plot(t, x_4_2)

Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'x_4_2_102_100k_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

figure
subplot(2,1,1)
plot(0:1/T:N/T - 1/T ,db(abs(fft(ymat(:,9)))))
title('output')
subplot(2,1,2)
plot(0:1/T:N/T - 1/T ,db(abs(fft(umat(:,9)))))
title('input')

figure
subplot(2,1,1)
plot(0:1/T:N/T - 1/T ,angle(fft(ymat(:,9))))
title('output')
subplot(2,1,2)
plot(0:1/T:N/T - 1/T ,angle(fft(umat(:,9))))
title('input')

%% 4.3

fs = 10e3; %Hz
f = 102.5391; %Hz
N = 4096;
T = N/fs;
T_period = 1/f; 
N_period = T/T_period; % how many samples is 1 period

t = 0:1/fs:(N-1)/fs;

x_4_3 = sin(2*pi*f*t);

figure
plot(t, x_4_3)

Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'x_4_2_102_100k_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

T = N_period/fs;
figure
subplot(2,1,1)
plot(0:1/T:N_period/T - 1/T ,db(abs(fft(ymat(1:N_period,Nrep)))))
subplot(2,1,2)
plot(0:1/T:N_period/T - 1/T ,db(abs(fft(umat(1:N_period,Nrep)))))

%% 4.4

fs = 10e3; %Hz
f =  102.5391; %Hz
N = 4096;
T = N/fs;
t = 0:1/fs:(N-1)/fs;
T_period = 1/f;
N_period = T/T_period; % how many samples is 1 period
amplitudes = log10(logspace(0.1, 1.1, 10));

x_4_4 = [];
for i = amplitudes
    x_4_4 = [ x_4_4 i*sin(2*pi*f*t)];
end
figure
plot(x_4_4)

N1 = length(x_4_4);
Nrep = 1; % how many repetitions were measured
Drep = 1; % how many of the last repetitions you want to keep
FileName = 'x_4_4_data.mat'; %name of the mat file
[umat , ymat] = ReadDataLab2(10*N, Nrep, Drep, FileName);
T = N_period/fs;
x = [];
y  = [];

for i = 1:length(amplitudes)
    figure
    subplot(2,1,1)
    temp = abs(fft(umat(N-N_period+1:N,1)));
    x = [x temp(2)];
    plot(0:1/T:(N_period-1)/T ,db(abs(fft(ymat(N-N_period+1:N,1)))))
    title('DFT of output with amplitude ' +string(amplitudes(i)))
    subplot(2,1,2)
    temp = abs(fft(ymat(N-N_period+1:N,1)));
    y = [y temp(2)];
    plot(0:1/T:(N_period-1)/T ,db(abs(fft(umat(N-N_period+1:N,1)))))
    title('DFT of input with amplitude '+string(amplitudes(i)))

    figure
    plot(0:1/T:(N_period-1)/T ,db(abs(fft(umat(N-N_period+1:N,1))/(fft(umat(N-N_period+1:N,1))))))
    title('FRF for amplitude '+ string(amplitudes(i)))

    ymat(1:N, :) = [];
    umat(1:N, :) = [];



end

%% 4.5

figure
plot((db(x)), db(y))
title('Compression of output power')
xlabel('input [dB]')
ylabel('output [dB]')
hold on
%slope = 7;
%x_inter = 25:70
%plot(x_inter, (x_inter - 27.412)*slope + 46.70)

% figure
% plot((db(x.^2)), db(y.^2))

%%

fs = 10e3; %Hz
f = 102.5391; %Hz
N = 4096;
T = N/fs;

t = 0:1/fs:(N-1)/fs;
amplitude = 0.5;

x_4_5 = amplitude*sin(2*pi*f*t);

figure
plot(x_4_5)
Nrep = 10; % how many repetitions were measured
Drep = 10; % how many of the last repetitions you want to keep
FileName = 'x_4_5_data.mat'; %name of the mat file

[umat , ymat] = ReadDataLab2(N, Nrep, Drep, FileName);

T = 100/fs;
N = 100;
figure
subplot(2,1,1)
plot(0:1/T:N/T - 1/T ,db(abs(fft(ymat(N-99:N,9)))))
title('output')
subplot(2,1,2)
plot(0:1/T:N/T - 1/T ,db(abs(fft(umat(N-99:N,9)))))
title('input')
