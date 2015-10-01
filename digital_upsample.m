function y = digital_upsample(x, p, q)

    if mod(p, q) == 0
        y = repmat(x, p/q, 1);
        y = y(:)';
    else
        N = floor(length(x)*p/q)+1;
        y = zeros(1, N);
        for i = 1:length(x)
            n = (ceil((i-1)*p/q):floor(i*p/q))+1;
            y(n) = x(i);
        end
    end

end