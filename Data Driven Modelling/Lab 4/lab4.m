clc
close all
clear

%% 3.1 Build the discrete time systme

%variables
length_impulse_response = 30; % length of the impulse response before it becomes neglegably small
nr_range = 100; % maximum order used for estimation
snr_dB = [600 60 6]; % for 3.5 make 6 -> 26
loop_length = 100; %set to 1 for everything before 3.6

[B, A] = cheby2(2, 3, [0.3 0.6], "bandpass"); % for 3.5 change to [0.4 0.6]?
F0 = tf(B, A, 1); % get discrete transfer function of the filter

%plot impulse response of transfer function
figure
impulse(F0)
title('Impulse response F0')

% get significant elements from the impulse response
G0 = impulse(F0);
G0 = G0(1:length_impulse_response);

% plot significant elements
figure
stem(G0)
title('Impulse response G0')
xlabel('Time (seconds)')
ylabel('Amplitude')

%% 3.2 Generate the exictation and output signals

% these are used for creating the histogram
V_ls_avg = zeros(length(snr_dB), nr_range, loop_length);
V_val_avg = zeros(length(snr_dB), nr_range, loop_length);
V_aic_avg = zeros(length(snr_dB), nr_range, loop_length);

for a = 1:loop_length % part 3.6
    disp(['loop iteration ' + string(a) + ' of ' + string(loop_length)])

    % generate random input signals
    u0_e = randn(50*length_impulse_response,1);
    u0_v = randn(50*length_impulse_response,1);
    
    % calculate output of input signals through filter
    y0_e = filter(G0, 1, u0_e);
    y0_v = filter(G0, 1, u0_v);
    
    % create arrays to add noise to the signal
    y0_e_noise = zeros(length(y0_e), length(snr_dB));
    y0_v_noise = zeros(length(y0_v), length(snr_dB));
    noise_e = zeros(length(y0_e), length(snr_dB));
    noise_v = zeros(length(y0_v), length(snr_dB));
    for i = 1:length(snr_dB)
        noise_e(:,i) = randn(size(y0_e))*std(y0_e)/db2mag(snr_dB(i)); %ensures SNR has wanted value
        y0_e_noise(:,i) = y0_e + noise_e(:,i);
        noise_v(:,i) = randn(size(y0_v))*std(y0_v)/db2mag(snr_dB(i));
        y0_v_noise(:,i) = y0_v + noise_v(:,i);
    end
    
%     % plot noisey signals
%     figure
%     hold on
%     for i = 1:length(snr_dB)
%         plot(y0_e_noise(:,i))
%     end
%     xlabel('Time (seconds)')
%     ylabel('Amplitude')
%     title('y0_e disturbed by noise')
%     legend(string(snr_dB)+' dB')
%     
%     figure
%     hold on
%     for i = 1:length(snr_dB)
%         plot(y0_v_noise(:,i))
%     end
%     xlabel('Time (seconds)')
%     ylabel('Amplitude')
%     title('y0_v disturbed by noise')
%     legend(string(snr_dB)+' dB')
    
    %% 3.3 Implement the least squares for varying n
    
    % calculates least squares
    V_ls = zeros(length(snr_dB), length_impulse_response);
    V_val = zeros(length(snr_dB), length_impulse_response);
    for nr = 1:nr_range
        % create Hn
        Hn_e = toeplitz(u0_e); % puts input elements on the diagonal
        Hn_v = toeplitz(u0_v);
        Hn_e = tril(Hn_e); % makes everything above diagonal zero
        Hn_v = tril(Hn_v);
        Hn_e = Hn_e(:,1:nr); % remove elements to get right dimension
        Hn_v = Hn_v(:,1:nr);
    
        % calculates least squares and validation for certain SNR levels
        for i = 1:length(snr_dB)
            theta_e =Hn_e\y0_e_noise(:,i);
            V_ls(i, nr, :) = 1/(length(y0_e)*var(noise_e(:,i)))*norm((y0_e_noise(:,i) - Hn_e*theta_e).^2);
            V_val(i, nr, :) = 1/(length(y0_v)*var(noise_v(:,i)))*norm((y0_v_noise(:,i) - Hn_v*theta_e).^2);
        end
    end
    
    % store data to plot in histogram
    V_ls_avg(:,:,a) = V_ls;
     V_val_avg(:,:,a) = V_val;

    
    %% 3.4 Select the optimal model order using AIC and Validation
    
    % calculate AIC
    V_aic = zeros(size(V_ls));
    for i = 1:length(snr_dB)
        V_aic(i,:, :) = V_ls(i,:, :).*(1 + 2*(1:nr)/length(y0_e));
    end

    % store data to plot in histogram
    V_aic_avg(:,:,a) = V_aic;
    


end

%% plots

    % plot least squares, validation and AIC
    figure
    hold on
    plot(V_ls_avg(3, :,1))
    plot(V_val_avg(3, :, 1))
    plot(V_aic_avg(3,:, 1))
    legend('V_{ls}', 'V_{val}', 'V_{aic}')

    % calculates the order for which the curves are minimal
    disp('Minimum for Least Squares: ')
    find((V_ls_avg(3,:, 1) == min(V_ls_avg(3,:, 1), [], 'all')))
    disp('Minimum for validation: ')
    find((V_val_avg(3,:, 1) == min(V_val_avg(3,:, 1), [], 'all')))
    disp('Minimum for AIC: ')
    find((V_aic_avg(3,:, 1) == min(V_aic_avg(3,:, 1), [], 'all')))

    % find minimal elements for each iteration
    ls_min = zeros(loop_length, 1);
    val_min = zeros(loop_length, 1);
    aic_min = zeros(loop_length, 1);
    for i = 1:loop_length
        ls_min(i) =  find((V_ls_avg(3,:, i) == min(V_ls_avg(3,:, i), [], 'all')));
        val_min(i) = find((V_val_avg(3,:, i) == min(V_val_avg(3,:, i), [], 'all')));
        aic_min(i) = find((V_aic_avg(3,:, i) == min(V_aic_avg(3,:, i), [], 'all')));
    end

    % plot histograms
    figure
    histogram(ls_min, 1:nr)
    title('Order where Least squares is minimal')
    xlabel('Order')

    figure
    histogram(val_min, 1:nr)
    title('Order where validation is minimal')
    xlabel('Order')

    figure
    histogram(aic_min, 1:nr)
    title('Order where AIC is minimal')
    xlabel('Order')

