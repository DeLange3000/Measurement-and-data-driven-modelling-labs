function [] = PlotData(umat, ymat, interesting_freqs, all_experiments, time_domain, frequency_domain, plot_title)

% umat: input data
% ymat: output data
% interesting_freqs: for which bins should data be plotted
% all_experiments: plots all experiments if set to true, otherwise it just
% plots the last one
% time domain: plots time domain signal if set to true
% frequency domain: plots frequency domain if set to true
% title: string that contains title for plots

[a,b] = size(umat);

if (all_experiments)
    Drep = 1:b;
else
    Drep = b;
end

for i = Drep
    freqU = fft(umat(:,i));
    freqY = fft(ymat(:,i));

    if (time_domain)
        figure
        plot(umat(:,i))
        hold on
        plot(ymat(:,i))
        legend('input','output')
        title('input and output of '+ plot_title)
    end

    if (frequency_domain)
        figure
        plot(interesting_freqs, db(abs(freqY(interesting_freqs)./freqU(interesting_freqs))))
        title('FRF of ', string(plot_title))
        ylabel('[dB]')
        xlabel('bins')
    end
end
