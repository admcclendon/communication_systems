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

figure;
plot(real(mod_symbols), imag(mod_symbols), '.');
axis([-1, 1, -1, 1]);

%% Channel Model
Nupsample = 8;
mod_upsample = [mod_symbols; zeros(Nupsample-1, length(mod_symbols))];
mod_upsample = mod_upsample(:)';

num = 5;
tx_filter = sinc((-num*Nupsample:num*Nupsample)/Nupsample);

mod_filtered = conv(tx_filter, mod_upsample);

channel_symbols = mod_filtered + 10^(-10/10)*randn(1, length(mod_filtered));

%% Plot Channel Spectrum
L = length(channel_symbols);
NFFT = 2^nextpow2(L);
Y_channel = 10*log10(abs(fftshift(fft(channel_symbols, NFFT)/L)));
f_channel = (1/2)*linspace(-1, 1, NFFT);

figure;
plot(f_channel, Y_channel);

%% Plot Channel IQ

figure;
subplot(2, 1, 1);
plot(real(channel_symbols));
subplot(2, 1, 2);
plot(imag(channel_symbols));

%% Matched filter response

rx_matched = conv(tx_filter, channel_symbols);

rx_symbols = rx_matched(Nupsample:Nupsample:end);

rx_bits = zeros(1, length(rx_symbols));
rx_bits(rx_symbols<0) = 1;

%% Plot Matched Filter

figure;
subplot(2, 1, 1);
plot(real(rx_matched));
subplot(2, 1, 2);
plot(imag(rx_matched));