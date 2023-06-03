clc
close all
clear

%% 1.2.1 -> 1.2.5

N = 1000; %this is equal to the amount of bins (bins = #samples in fft)
Ts = 1e-2;
f = 3/(N*Ts);
x = 0:N-1;

y = cos(2*pi*Ts*x*f + 2*pi*rand); %3/N for 3 whole periods in window

figure
plot(x.*Ts, y);
xlabel("time [s]")
ylabel('amplitude')

figure
plot(linspace(0, 1/(Ts), N), db(fft(y)), 'o') % freq of cosine is correct!!
xlabel("Frequency (Hz)")
ylabel("Amplitude (dB)")

figure
plot(0:N-1, db(fft(y)), 'o') % freq of cosine is correct!!
xlabel("Frequency (bin)")
ylabel("Amplitude (dB)")

Y = fft(y);
tolerance = 1e-6;
Y(abs(Y)  < tolerance) = 0;

figure
plot(linspace(0, 1/(Ts), N), angle(Y), 'o')
xlabel('frequency [Hz]')
ylabel("angle [rad]")

%% 1.3.1 -> 1.3.2

N = 1000;
K = 10;
A = 1;
phases = 2*pi*rand(1, K);
Fs = 100; %Hz
Ts = 1/Fs;

f = (1:K)/(N*Ts);

x = 1:N;
y = 0;
for i = 1:length(f)
    y = y + A*sin(2*pi*Ts*f(i)*x + phases(i));
end
figure
plot(x.*Ts, y)
xlabel('time [s]')
ylabel("Amplitude")

Y = fft(y);
figure
plot(0:N-1, db(Y), 'o')
xlabel('frequency [bins]')
ylabel('Amplitude')

figure
plot(linspace(0, 1/(Ts), N), db(Y), 'o') % freq of cosine is correct!!
xlabel("Frequency (Hz)")
ylabel("Amplitude (dB)")

%% 1.3.3

N = 1000;
Fs = 200;
Ts = 1/Fs;
f = [4, 8, 12, 16, 20, 24]; %Hz
A = 1;
phases = 2*pi*rand(1, length(f));

x = 0:N-1;
y = 0;
for i = 1:length(f)
    y = y + sin(2*pi*f(i)*Ts*x + phases(i));
end

figure
plot(x.*Ts, y)
xlabel('time [s]')
ylabel("Amplitude")

Y = fft(y);
figure
plot(0:N-1, db(Y), 'o')
xlabel('frequency [bins]')
ylabel('Amplitude')

figure
plot(linspace(0, 1/(Ts), N), db(Y), 'o') % freq of cosine is correct!!
xlabel("Frequency (Hz)")
ylabel("Amplitude (dB)")

%% 1.4

N = 1000;
X = zeros(1, N);
A = 1;
phases = 2*pi*rand(1,30);
X(1:30) = 1*exp(1i*phases);

figure
plot(abs(X), 'o')

x = ifft(X);
figure
plot(N*real(x))

%% 1.4.3

f_excited = linspace(5, 15, 30);
fs = 100; %anything above 60 Hz is okay
Ts = 1/fs;
N = 1000;
A = 1;
f = zeros(1, N);
f_excited = 30;
%f(f_excited*N/fs+1) = A*exp(1i*2*pi*rand); % use this to select right freq
f(N/fs*5+1:N/(2*fs):N/fs*5 + 30*N/(2*fs)) = A*exp(1i*2*pi*rand(1, 30));

figure
plot(linspace(0, 1/(Ts), N), abs(f), 'o') % freq of cosine is correct!!
xlabel("Frequency (Hz)")
ylabel("Amplitude (dB)")

figure
plot((0:N-1).*Ts, N*real(ifft(f)))
xlabel('time [s]')
ylabel("Amplitude")

%% 1.5

N = 500;
K = 60;
f = (1:K)/N;
A = 1;

x = 0:N-1;
y1 = 0;
y2 = 0;
y3 = 0;

for i = 1:length(f)
    y1 = y1 + A*sin(2*pi*f(i)*x + 2*pi*rand);
    y2  = y2 + A*sin(2*pi*f(i)*x + i*(i+1*pi/K));
    y3 = y3 + A*sin(2*pi*f(i)*x + i*pi);
end

figure
subplot(3, 1, 1)
plot(x, y1);
xlabel("time [s]")
ylabel('amplitude')
subplot(3, 1, 2)
plot(x, y2);
xlabel("time [s]")
ylabel('amplitude')
subplot(3, 1, 3)
plot(x, y3);
xlabel("time [s]")
ylabel('amplitude')

figure
subplot(3, 1, 1)
plot(x, db(fft(y1)), 'o') % freq of cosine is correct!!
xlabel("Frequency (bin)")
ylabel("Amplitude (dB)")
subplot(3, 1, 2)
plot(x, db(fft(y2)), 'o') % freq of cosine is correct!!
xlabel("Frequency (bin)")
ylabel("Amplitude (dB)")
subplot(3, 1, 3)
plot(x, db(fft(y3)), 'o') % freq of cosine is correct!!
xlabel("Frequency (bin)")
ylabel("Amplitude (dB)")

Crest_random = max(abs(y1))/rms(y1)
Crest_schroeder = max(abs(y2))/rms(y2)
Crest_linear = max(abs(y3))/rms(y3)


