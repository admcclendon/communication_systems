alpha = 0.1;
tau1 = .65;
tau2 = .3;
K = 2;
b = 2*[(1+2*tau1), 2, (1-2*tau1)];
a = [(tau2*4+K+tau1*K*2), (2*K-tau2*8), (tau2*4+K-tau1*K*2)];

N = 50;
x = ones(1, N)+cumsum(ones(1, N));

y = filter(b, a, x);

figure;
plot(x-y)