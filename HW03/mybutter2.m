clear all;
clc;

%% Part 2 (a)
Ts = 0.002;
fs = 1/Ts;
f1 = 100;
M = 1000;
n = linspace(1, M, M);
x = zeros(1, M);

for m = 1:M
    x(m) = cos(2 * pi * (m - 1) * Ts) + 2 * cos(2 * pi * f1 * (m - 1) * Ts);
end

figure
subplot(3, 1, 1)
plot(n, x)
title('x[n] = cos(2π(m - 1)Ts) + 2cos(2πf1(m - 1)Ts)')
xlabel('n')
ylabel('x[n]')

%% Part 2 (b)
L = 16; % order
fc = 0.5 * f1/(0.5 * fs);
[b1, a1] = butter(L, fc);
[H1, w1] = freqz(b1, a1, M);

y1 = filter(b1, a1, x);

subplot(3, 1, 2)
%plot(n, y1, 'b', n, cos(2 * pi * (n - 1) * Ts), 'r')
plot(n, y1, 'b')
title('y[n] = cos(2π(m - 1)Ts)')
xlabel('n')
ylabel('y[n]')

%% Part 2 (c)
fc = f1 / (0.5 * fs);
[b2, a2] = butter(L, [fc * 0.707  fc * 1.414], 'bandpass');
[H2, w2] = freqz(b2, a2, M);

y2 = filter(b2, a2, x);

subplot(3, 1, 3)
%plot(n, y2, 'b', n, 2* cos(2 * pi * f1 * (n - 1) * Ts), 'r')
plot(n, y2, 'b')
title('y[n] = 2cos(2πf1(m - 1)Ts)')
xlabel('n')
ylabel('y[n]')