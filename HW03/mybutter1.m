clear all;
clc;

%% Part 1 (a)
N = 100;
fs = 20;
Ts = 1/fs;
n = linspace(1, 100,N);
x = zeros(1, N);
for m = 1:N
    x(m) = cos(2 * pi * (m - 1) * Ts);
end

figure
plot(n, x)
title('x[n] = cos(2π(n-1)Ts)')
xlabel('n')
ylabel('x[n]')

%% Part 1 (b)
L = 3; % order
fc = 0.05; % Normalized Cutoff Frequency
[b1, a1] = butter(L, fc);
[H, w] = freqz(b1, a1, N); % Frequency response H(e^jw)

figure
subplot(3, 3, 1)
plot(w, 20*log10(abs(H)))
title('Magnitude of H(e^{jω})')
xlabel('Normalized frequency ω (rad/sample)')
ylabel('20log|H(e^{jω})| (dB)')

subplot(3, 3, 4)
plot(w, unwrap(angle(H))*180/pi)
title('Phase of H(e^{jω})')
xlabel('Normalized frequency ω (rad/sample)')
ylabel('Phase (deg)')

y = filter(b1, a1, x);

subplot(3, 3, 7)
plot(n, y);
title('y[n] = x[n] * h[n]')
xlabel('n')
ylabel('y[n]')

%% Part 1 (c)
L = 7;
fc = 0.05;
fs = 20;

[b2, a2] = butter(L, fc);
[H, w] = freqz(b2, a2, N); % Frequency response H(e^jw)

subplot(3, 3, 2)
plot(w, 20*log10(abs(H)))
title('Magnitude of H(e^{jω})')
xlabel('Normalized frequency ω (rad/sample)')
ylabel('20log|H(e^{jω})| (dB)')

subplot(3, 3, 5)
plot(w, unwrap(angle(H))*180/pi)
title('Phase of H(e^{jω})')
xlabel('Normalized frequency ω (rad/sample)')
ylabel('Phase (deg)')

y = filter(b2, a2, x);

subplot(3, 3, 8)
plot(n, y);
title('y[n] = x[n] * h[n]')
xlabel('n')
ylabel('y[n]')

%% Part 1 (d)
L = 3;
fc = 0.5;
fs = 20;

[b3, a3] = butter(L, fc);
[H, w] = freqz(b3, a3, N); % Frequency response H(e^jw)

subplot(3, 3, 3)
plot(w, 20*log10(abs(H)))
title('Magnitude of H(e^{jω})')
xlabel('Normalized frequency ω (rad/sample)')
ylabel('20log|H(e^{jω})| (dB)')

subplot(3, 3, 6)
plot(w, unwrap(angle(H))*180/pi)
title('Phase of H(e^{jω})')
xlabel('Normalized frequency ω (rad/sample)')
ylabel('Phase (deg)')

y = filter(b3, a3, x);

subplot(3, 3, 9)
plot(n, y);
title('y[n] = x[n] * h[n]')
xlabel('n')
ylabel('y[n]')

%% Part 1 (e)
% See the report for explanation!
