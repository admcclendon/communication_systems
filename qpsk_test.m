%% Modulate at baseband
close all;
clear;
clc;

N = 1000; % number of bits

bits = randi([0 1], 1, N);

symbols = (reshape(bits, 2, length(bits)/2)'*[1;2])';

mod_symbols = exp(j*((2*symbols+1)*pi/4));

%% Constellation of modulated symbols
figure; plot(real(mod_symbols),imag(mod_symbols), '.')

%% Upsample and filter
Nupsample = 3;
mod_symbols_upsample = [mod_symbols; zeros(Nupsample-1, length(mod_symbols))];
mod_symbols_upsample = mod_symbols_upsample(:)';

mod_symbols_filtered = conv(ones(1, Nupsample), mod_symbols_upsample);

%% Spectrum
figure; pwelch(mod_symbols_filtered);

%% Channel model (10dB SNR)
% channel_impulse = rand(1, 20);
% channel_impulse = [channel_impulse(1:length(channel_impulse)/2) 1.5 channel_impulse((length(channel_impulse)/2+1):end)];
% channel_impulse = channel_impulse/max(channel_impulse); % normalize
channel_impulse = [0.01, -0.02, 0.05, -0.1, 0.2, 1, 0.15, -0.15, 0.05, -0.02, 0.005];

figure; stem(channel_impulse);
channel_symbols = conv(channel_impulse, mod_symbols);

channel_symbols_filtered = 10^(-40/10)*conv(channel_impulse, mod_symbols_filtered);
channel_symbols_filtered = channel_symbols_filtered + 10^(-50/10)*exp(j*1)*randn(1, length(channel_symbols_filtered));

%% Constellation of channel symbols
figure; plot(real(channel_symbols), imag(channel_symbols), '.')

axis([-1.5, 1.5, -1.5, 1.5]);
axis square;

%% Constellation of channel filtered symbols
figure; plot(real(channel_symbols_filtered), imag(channel_symbols_filtered), '.');


%% Spectrum
figure; pwelch(channel_symbols)

%% Equalize
Neq = 2;
Pc = zeros(2*Neq+1,2*Neq+1);
for k = 1:2*Neq+1
    Pc(k, :) = channel_impulse((0:-1:-2*Neq)+2*Neq+k+(length(channel_impulse)-(4*Neq+1))/2);
end
Peq = [zeros(Neq,1);1;zeros(Neq,1)];

c = (inv(Pc)*Peq)';

y = conv(c, channel_symbols);
figure; plot(real(y), imag(y), '.')

%% AGC

agc_symbols = channel_symbols_filtered;
G = 0; % dB
Navg = 10;
alpha = 1000;
Pd = -20;
mu = 0.5;
errors = [];
Gs = [];
for i = 0:floor(length(agc_symbols)/Navg)-1
    index = (i*Navg+1):(i+1)*Navg;
    agc_symbols(index) = 10^(G/20)*channel_symbols_filtered(index);
    
    error = mu*(Pd - 10*log10(sum(abs(agc_symbols(index)).^2)/Navg));
    G = G + error;
    
    errors = [errors error];
    Gs = [Gs G];
end

figure;
subplot(3, 1, 1);
plot((1:length(agc_symbols))/Nupsample, 10*log10(abs(agc_symbols)));
title(['Num of samples ' num2str(Navg) ', mu ' num2str(mu) ', Magnitude of received signal']);
subplot(3, 1, 2);
stem(Gs);
title('Receiver Gain (dB)');
subplot(3, 1, 3);
stem(errors);
title('Error Signal');
