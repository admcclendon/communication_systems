function f = plot_constellation(x, ti)
    
    if nargin == 1
        ti = 'Constellation';
    else
        ti = ['Constellation of ' ti];
    end
    f = figure;
    plot(real(x), imag(x), '.');
    title(ti);
    axis([-1 1 -1 1]);
    
end