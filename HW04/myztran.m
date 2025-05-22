clear all;
clc;

%% Question (a)
k0 = 0.09; % Gain
z0 = [1; 1; -1; -1]; % Zeros
p0 = [0.3 - 0.4 * 1i; 0.3 + 0.4 * 1i; 0.1 - 0.1 * 1i; 0.1 + 0.1 * 1i]; % Poles
subplot(2, 4, 1)
zplane(z0, p0) % Generate Pole-zero Plot


%% Question (b)
[b0, a0] = zp2tf(z0, p0, k0);
[r, p, k] = residuez(b0, a0);
N0 = 20;
n = linspace(0, N0, N0 + 1);
h = zeros(1, N0 + 1);
subplot(2, 4, 2)
for x = 1 : N0 + 1
    if x <= length(k)
        h(x) = h(x) + k(1, x);
    end
    for m = 1 : length(p)
        h(x) = h(x) + r(m, 1) * (p(m, 1)^(x - 1));
    end
end
stem(n, h(n + 1), 'filled');
title('h[n], inverse z-Transform of H(z)')
xlabel('n')
ylabel('h[n]')

%% Question (c)
[H, w] = freqz(b0, a0, N0 + 1);

subplot(2, 4, 3)
plot(w, 20*log10(abs(H)));
title('Magnitude Response of H(z)')
xlabel('Normalized Frequency ω (rad/sample)')
ylabel('20log|H(z)| (dB)')

subplot(2, 4, 4)
plot(w, unwrap(angle(H))*180/pi);
title('Phase Response of H(z)')
xlabel('Normalized Frequency ω (rad/sample)')
ylabel('Phase (deg)')

%% Question (d)
sos = zp2sos(z0, p0, k0);

%% Question (e)
subplot(2, 4, 5)
[H1, w] = freqz(sos(1, 1:3), sos(1, 4:6), N0 + 1);
plot(w, 20*log10(abs(H1)));
title('Magnitude Response of H_1(z)')
xlabel('Normalized Frequency ω (rad/sample)')
ylabel('20log|H_1(z)| (dB)')

subplot(2, 4, 6)
[H2, w] = freqz(sos(2, 1:3), sos(2, 4:6), N0 + 1);
plot(w, 20*log10(abs(H2)));
title('Magnitude Response of H_2(z)')
xlabel('Normalized Frequency ω (rad/sample)')
ylabel('20log|H_1(z)| (dB)')

subplot(2, 4, 7)
plot(w, 20*(log10(abs(H1.*H2))));
title('Magnitude Response of H_1(z)H_2(z)')
xlabel('Normalized Frequency ω (rad/sample)')
ylabel('20log|H_1(z)H_2(z)| (dB)')

%% Question (f)
x_impulse = zeros(1, N0 + 1);
x_impulse(1) = 1;
y = filter(b0, a0, x_impulse);

subplot(2, 4, 8)
stem(n, y(n + 1), 'filled');
title('y[n], Inverse z-Transform of Y(z) = H(z)X(z)');
xlabel('n');
ylabel('y[n]');

%%figure
%%plot(w, 20*log10(abs(H)), 'b', w, 20*log10(abs(H1 .* H2)), 'r')
%%title('Comparison of H(z) and H_1(z)H_2(z)')
%%xlabel('Normalized Frequency ω (rad/sample)')
%%ylabel('Magnitude (dB)')

%%figure
%%stem(n, h(n+1), 'filled', 'b')
%%hold on;
%%stem(n, y(n+1), 'filled', 'r')
%%title('Comparison of h[n] and y[n]')
%%xlabel('n')
%%ylabel('h[n] (blue) and y[n] (red)')
