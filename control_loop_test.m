alpha = 0.1;

b = alpha*[1 1];
a = [(alpha + 1), (alpha - 1)];

x = 5*ones(1, 200) + cumsum(0.01*ones(1, 200));

y = filter(b, a, x);

figure;
plot(x-y)