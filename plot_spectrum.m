function f = plot_spectrum(x, ti)

    if nargin == 1
        ti = 'Spectrum';
    else
        ti = ['Spectrum of ' ti];
    end
    L = length(x);
    NFFT = 2^nextpow2(L);
    Y = 10*log10(abs(fftshift(fft(x, NFFT)/L)));
    freq = (1/2)*linspace(-1, 1, NFFT);

    f = figure;
    plot(freq, Y);
    title(ti);
end