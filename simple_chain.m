%% Generate BPSK signal
clear;
clc;

N = 1000;
pbrs = PRBS15();
bits = pbrs.Generate(N);

M = 2;
symbols = (1:log2(M))*reshape(bits, log2(M), length(bits)/log2(M));

mod_symbols = exp(j*(2*pi*symbols/M));

%% Plot Constellation
plot_constellation(mod_symbols, 'BPSK Modulated Symbols');

%% Channel Model
Nupsample = 40;
mod_upsample = [mod_symbols; zeros(Nupsample-1, length(mod_symbols))];
mod_upsample = mod_upsample(:)';

num = 5;
tx_filter = sinc((-num*Nupsample:num*Nupsample)/Nupsample);

mod_filtered = conv(tx_filter, mod_upsample);

carrier = exp(j*2*pi*10*(0:length(mod_filtered)-1)/Nupsample);

channel_symbols = mod_filtered.*carrier + 10^(-10/10)*randn(1, length(mod_filtered));

%% Plot Channel Spectrum
plot_spectrum(channel_symbols, 'Channel Symbols');

%% Plot Channel IQ
plot_iq(channel_symbols, 'Channel Symbols');

%% IF Mixing Stage
% rx_lo = exp(-j*2*pi*((0:length(channel_symbols)-1)+1)*0.25);
rx_lo = zeros(1, length(channel_symbols));
rx_mixed = zeros(1, length(channel_symbols));
phase_detect = zeros(1, length(channel_symbols));
err = zeros(1, length(channel_symbols));
% [loop_b, loop_a] = butter(2, 0.6);
loop_b = 0.1*[1 1];
loop_a = [1 -1];
loop_z = [];
loop_z2 = [];
for n = 1:length(channel_symbols)
    t = (n - 1)/Nupsample;
    if (n > 1)
        rx_lo(n) = exp(-j*(2*pi*9.95*t + err(n-1) + 2));
    else
        rx_lo(n) = exp(-j*(2*pi*10*t));
    end
    rx_mixed(n) = rx_lo(n)*channel_symbols(n);
    
    phase_detect(n) = real(rx_mixed(n))*imag(rx_mixed(n));
    [err(n) loop_z] = filter(loop_b, loop_a, phase_detect(n), loop_z);
%     [err(n) loop_z2] = filter(loop_b2, loop_a2, err(n), loop_z2);
%     if (n > 1)
%         err(n) = phase_detect(n)*4 + err(n-1);
%     else
%         err(n) = phase_detect(n)*4;
%     end
end
% rx_mixed = rx_lo.*channel_symbols;
%% Plot Error
figure;
plot(err);
% plot_iq(rx_lo, 'RX LO');

%% Plot Spectrum
plot_spectrum(rx_mixed, 'RX after IF Mixing');

%% Matched filter response
rx_matched = conv(tx_filter, rx_mixed);

rx_symbols = rx_matched(Nupsample:Nupsample:end);

rx_bits = zeros(1, length(rx_symbols));
rx_bits(rx_symbols<0) = 1;

%% Plot Matched Filter
plot_iq(rx_matched, 'RX Matched Filter');

%% Plot Correlation between TX bits and RX bits
figure;
plot(conv(2*bits(end:-1:1)-1, 2*rx_bits-1));
title(['Correlation of TX bits vs RX bits; N = ' num2str(N)]);

%%
freqz(loop_b, loop_a)