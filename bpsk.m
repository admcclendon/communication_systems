close all
clear

%% Basic Parameters
Fs = 1e6;
Ts = 1/Fs;

Fb = 10e3;
Tb = 1/Fb;

N = 20;

fc = 100e3;
fo = 0;
t = 0:Ts:N*Tb-Ts;

%% Modulation @ baseband
bits = randi([0 1], 1, N);

z = exp(j*pi*bits);
Z_upsample = digital_upsample(z, Fs, Fb);

% Pulse shaping
yx = ones(1,Fs/Fb);

%tx = -3*Tb:Ts:3*Tb;
%yx = sinc(tx/Tb);

Z_shaped = filter(yx, 1, Z_upsample);

%% Mix up to fc
y = real(Z_shaped).*cos(2*pi*(fc+fo)*t) + imag(Z_shaped).*sin(2*pi*(fc+fo)*t);

%% Plot Modulated Signal with AWGN
figure; plot(t*1e3, y);
title('Modulated BPSK Signal')
xlabel('Time (ms)')
ylabel('Signal Level')

%% Calculate Spectrum
L = length(y);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

%% Plot single-sided amplitude spectrum.
figure; plot(f,10*log10(abs(Y(1:NFFT/2+1)))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)| (dB)')
axis([fc-4*Fb, fc+4*Fb, -50, 0])

%% Demodulate signal and filter
r = y.*cos(2*pi*fc*t);

[b, a] = butter(9, Fb/Fs);
rf = filter(b, a, r);

%% Calculate I/Q Samples
r_i = y.*cos(2*pi*fc*t);
r_q = y.*sin(2*pi*fc*t);

r_i_filtered = filter(b, a, r_i);
r_q_filtered = filter(b, a, r_q);

%% Plot I/Q Samples
figure;
subplot(2, 1, 1);
plot(t, r_i_filtered);
title('Inphase Samples');

subplot(2,1,2);
plot(t, r_q_filtered);
title('Quadrature Samples');

%% Plot demodulated signal
figure; plot(t*1e3, rf);
title('Demodulated Signal');
xlabel('Time (ms)');
ylabel('Signal Level');

%% Matched Filter
b_avg = [0.5, ones(1, Fs/Fb - 2), 0.5];
r_m = filter(b_avg, 1, rf);

figure; plot(t*1000, r_m);
title('Matched Filter Output');
xlabel('Time (ms)');
ylabel('Signal Level');

%% Resample matched filter output
bits_recv = r_m(floor((1:N)*Tb*Fs));
bits_recv(bits_recv > 0) = 0;
bits_recv(bits_recv < 0) = 1;
%for i = 1:length(bits_recv)
%    matched_sample = r_m(floor(i*Tb*Fs));
%    if matched_sample > 0
%        bits_recv(i) = 0;
%    else
%        bits_recv(i) = 1;
%    end
%end
sum(abs(bits-bits_recv))


r_z = r_i_filtered + j*r_q_filtered;
bits_upsample = repmat(bits, Fs/Fb, 1);
figure; plot(diff(unwrap(angle(r_z)))); hold on; plot(bits_upsample(:)); hold off;