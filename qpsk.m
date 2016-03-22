%% Generate QPSK signal
clear;
clc;

N = 1000;

pbrs = PRBS15();

bits = pbrs.Generate(N);

M = 4;
symbols = (1:log2(M))*reshape(bits, log2(M), length(bits)/log2(M));

mod_symbols = exp(j*(2*pi*symbols/M + pi/4));

%% Plot constellation
figure; 
plot(real(mod_symbols), imag(mod_symbols), '.');
title('QPSK Modulated Symbols');
xlabel('Real');
ylabel('Imag');

%% Upsample and filter
Nupsample = 8;
mod_upsampled = [mod_symbols; zeros(Nupsample-1, length(mod_symbols))];
mod_upsampled = mod_upsampled(:)';

num_of_periods = 6;
tx_filter = sinc((-Nupsample*num_of_periods:Nupsample*num_of_periods)/Nupsample);
mod_filtered = conv(tx_filter, mod_upsampled);

%% Plot Spectrum of upsampled and filtered
L = length(mod_upsampled);
NFFT = 2^nextpow2(L);
Y_upsampled = fft(conv(ones(1, Nupsample), mod_upsampled), NFFT)/L;
f_upsampled = (1/2)*linspace(-1, 1, NFFT);

figure;
subplot(2, 1, 1);
plot(f_upsampled, 10*log10(abs(fftshift(Y_upsampled))));

L = length(mod_filtered);
NFFT = 2^nextpow2(L);
Y_filtered = fft(mod_filtered, NFFT)/L;
f_filtered = (1/2)*linspace(-1, 1, NFFT);

subplot(2, 1, 2);
plot(f_filtered, 10*log10(abs(fftshift(Y_filtered))));

%% Plot tx_filter
figure;
stem(tx_filter);
title('Impulse Response of TX Filter');

figure;
subplot(2, 1, 1);
plot(real(mod_filtered));
title('Upsampled Modulated Symbols');
subplot(2, 1, 2);
plot(imag(mod_filtered));
title('Filtered Modulated Symbols');

%% Channel Model

carrier = cos(2*pi*(1:length(mod_filtered))/4);

noisePower = -30;
noise = 10^(noisePower/10)*randn(1, length(mod_filtered));

channel_symbols = mod_filtered.*carrier + noise;

%% Plot channel constellation
figure;
plot(real(channel_symbols), imag(channel_symbols), '.');

%% Plot Spectrum
L = length(channel_symbols);
NFFT = 2^nextpow2(L);
Y_channel = fft(channel_symbols, NFFT)/L;
f_channel = (1/2)*linspace(-1, 1, NFFT);

figure;
plot(f_channel, 10*log10(abs(fftshift(Y_channel))));

%% Costas Loop

% rx_lo = exp(j*2*pi*(0:length(channel_symbols)-1)/4);
[b, a] = butter(10, 0.25);
% rx_samples = filter(b, a, rx_lo.*channel_symbols);

% rx_filtered = conv(tx_filter, channel_symbols);
rx_filtered = channel_symbols;
[loop_b, loop_a] = butter(3, 0.1);
vco_in = 0;
z_loop = [];
z_lpf = [];
out = zeros(1, length(rx_filtered));
err = zeros(1, length(rx_filtered));
vco_actual = zeros(1, length(rx_filtered));
for i = 1:length(rx_filtered)
    t = (i-1)/Nupsample;
    vco_actual(i) = exp(j*2*pi*t*(vco_in + 0.25));
    
    [out(i), z_lpf] = filter(b, a, rx_filtered(i)*vco_actual(i), z_lpf);
    err(i) = real(out(i))*imag(out(i));
    [loop_out, z_loop] = filter(loop_b, loop_a, err(i), z_loop);
    vco_in = loop_out+vco_in;
end

%%
figure;
subplot(2, 1, 1);
plot(real(out));
subplot(2, 1, 2);
plot(imag(out));

L = length(out);
NFFT = 2^nextpow2(L);
Y_out = fft(out, NFFT)/L;
f_out = (1/2)*linspace(-1, 1, NFFT);

figure;
plot(f_out, 10*log10(abs(fftshift(Y_out))));

%%
figure;
plot(err)

figure;
plot(1:length(vco_actual), real(vco_actual), 'b', 1:length(vco_actual), imag(vco_actual),'g');
%%
L = length(rx_filtered);
NFFT = 2^nextpow2(L);
Y_rx_filtered = fft(rx_filtered, NFFT)/L;
f_rx_filtered = (1/2)*linspace(-1, 1, NFFT);

figure;
plot(f_rx_filtered, 10*log10(abs(fftshift(Y_rx_filtered))));

%%
figure;
subplot(2, 1, 1);
plot(real(rx_filtered));
subplot(2, 1, 2);
plot(imag(rx_filtered));

%% plot rx_lo

figure;
subplot(2, 1, 1);
plot(real(rx_lo));
subplot(2, 1, 2);
plot(imag(rx_lo));

%% plot vco error signal

figure;
plot(filter(b, a, real(rx_filtered).*imag(rx_filtered)));