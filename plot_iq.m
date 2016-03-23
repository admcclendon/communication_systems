function f = plot_iq(x, ti)
    if nargin == 1
        ti1 = 'I Channel';
        ti2 = 'Q Channel';
    else
        ti1 = ['I Channel of ' ti];
        ti2 = ['Q Channel of ' ti];
    end
    
    f = figure;
    subplot(2, 1, 1);
    plot(real(x));
    title(ti1);
    subplot(2, 1, 2);
    plot(imag(x));
    title(ti2);
end