clear all;
clc;

%% Question (a)
T = 100; 
N1 = 500;
T_s = T/N1;
N = 2 * N1 + 1;
n = linspace(-N1, N1, N);
x = sinc(2*T_s*n);
omega = 2 * pi / (2 * N1);

figure % Figure 1 - x[n] = sin(2πTs*n)/2πTs*n
plot(n, x) % Generate Continuous Signal
% stem(n, x, 'filled') % Generate Discrete Signal
title('x[n] = sin(2πTs*n)/2πTs*n')
xlabel('n')
ylabel('x[n]')

%% Question (b)
X_jw = fft(x);
f = linspace(-omega, omega, N);

figure % Figure 2 - Fourier transform of x[n] using fft
plot(f, abs(fftshift(X_jw)))
title('Fourier transform of x[n] using fft')
xlabel('ω (rad/s)')
ylabel('|X(e^{jω})|')

%% Question (c)
X_s = zeros(1, N); % Implementation of DTFT

for k = -N1 : N1
    for m = -N1 : N1
        X_s(1, k + N1 + 1) = X_s(1, k + N1 + 1) + x(m + N1 + 1) * exp(-1j * k * omega * m);
    end
end

% For comparing fft and Fourier transform from definition
% figure
% plot(f, abs(fftshift(X_jw)), 'b', f, abs(X_s), 'r')
% xlabel('ω (rad/s)')
% ylabel('|X(e^{jω})|')

figure % Figure 3 - Fourier transform of x[n], no fft
plot(f, abs(X_s))
title('Fourier transform of x[n], no fft')
xlabel('ω (rad/s)')
ylabel('|X(e^{jω})|')


%% Question (d)
T_w = T/2;
w = 0.5*(1 + cos(2*pi*T_s*n/T_w));
for m = -N1:N1
    if abs(m * T_s) > T_w / 2
        w(m + N1 + 1) = 0;
    end    
end

figure % Figure 4 - w[n]
plot(n, w)
title('w[n] = 0.5 * (1 + cos(2π|t|/Tw))')
xlabel('n')
ylabel('w[n]')

%% Question (e)
y = x.*w; % element-wise multiplication

% For comparing x[n] and y[n]
% figure
% plot(n, x, 'b', n, y, 'r')
% xlabel('n')

figure % Figure 5 - y[n] = x[n]w[n]
plot(n, y);
title('y[n] = x[n]w[n]')
xlabel('n')
ylabel('y[n]')

%% Question (f)
y_fft = fft(y);

figure % Figure 6 - Fourier transform of y[n] = x[n]w[n] using fft
plot(f, abs(fftshift(y_fft)))
title('Fourier transform of y[n] = x[n]w[n] using fft')
xlabel('ω (rad/s)')
ylabel('|Y(e^{jω})|')
