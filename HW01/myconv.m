clear all; 
clc;

%% Question (a)
n1 = 50; % number of elements in x1
n2 = 50; % number of elements in x2

x1 = zeros(n1, 1);
x2 = zeros(n2, 1);

for n = 1:n1
    if n <= 20
        x1(n, 1) = n;
    elseif n <= 39
        x1(n, 1) = 40 - n;
    else
        x1(n, 1) = 0;
    end
end


for n = 1:n2
    if n < 11
        x2(n, 1) = 1;
    else
        x2(n, 1) = 0;
    end
end

% Plot x1[n]
stem(x1, 'filled')
title('x1[n] = n, 1 <= n <= 20 \\ 40 - n, 21 <= n <= 39')
figure;

% Plot x2[n]
stem(x2, 'filled')
title('x2[n] = u[n - 1] - u[n - 11]')
figure;

%% Question (b)

y = conv(x1, x2);

% Plot conv(x1, x2)
stem(y, 'filled')
title('y[n] = conv(x1, x2)')
figure;

% Supp 1. Generate conv(x1, x2) and myconv(x1, x2) on the same graph
% stem(y, 'filled')
% title('conv(x1, x2) and myconv(x1, x2)')
% hold on;

% Supp 2. Shift x2 to obtain correct convolution
% x2_shifted = zeros(n2, 1);
% for n = 2:n2
%     x2_shifted(n) = x2(n - 1);
% end
% y_modified = conv(x1, x2_shifted);
% tiledlayout(2, 1)
% nexttile
% stem(y_modified, 'filled','r')
% title('y[n] = conv(x1, x2__modified)')

%% Question (c)

y_matrix = zeros(n1 + n2, 1);
x1_matrix = zeros(n1 + n2 - 1, n1);
x2_matrix = zeros(n2, 1);

for n = 1:n1
    for m = n:(n + n1 - 1)
        x1_matrix(m, n) = x1(m - n + 1);
    end
end

for m = 1:n2
    x2_matrix(m, 1) = x2(m);
end

y_matrix(2:n1 + n2, 1) = x1_matrix * x2_matrix;

% stem(y_matrix, 'filled')
% legend
% hold off;
% figure;
% Supp 2. end

% nexttile
% stem(y_matrix, 'filled')
% title('y[n] = myconv(x1, x2)')
% figure;
% Supp 1. end

% Plot myconv(x1, x2) = x1[n] * x2[n]
stem(y_matrix, 'filled')
title('y[n] = myconv(x1, x2) = x1[n] * x2[n]')
figure;

%% Question (d)

for n = 1:n1
    if n <= 3
        x1(n) = 3^(n);
    else
        x1(n) = 0;
    end
end

for n = 1:n2
    if n <= 5
        x2(n) = 2^(n);
    else
        x2(n) = 0;
    end 
end

% Plot x1[n] = 3^(n), 1 <= n <= 3
stem(x1, 'filled')
title('x1[n] = 3^n, 1 <= n <= 3')
figure;

% Plot x2[n] = 2^(n), 1 <= n <= 5
stem(x2, 'filled')
title('x2[n] = 2^n, 1 <= n <= 5')
figure;

y = conv(x1, x2);

% Plot y[n] = conv(x1[n], x2[n])
% stem(y, 'filled')
% title('y[n] = conv(x1, x2)')
% figure;

% Supp 3. Shift x2 to obtain correct convolution
% x2_shifted = zeros(n2, 1);
% for n = 2:n2
%     x2_shifted(n) = x2(n - 1);
% end
% y_modified = conv(x1, x2_shifted);
% tiledlayout(2, 1)
% nexttile
% stem(y_modified, 'filled', 'r')
% title('y[n] = conv(x1, x2__modified)')


% Supp 4. Generate conv(x1, x2) and myconv(x1, x2) on the same graph
stem(y, 'filled')
title('conv(x1, x2) and myconv(x1, x2)')
hold on;


for m = 1:(n1 + n2 - 1)
    for n = 1:n1
        x1_matrix(m, n) = 0;
    end
end

for m = 1:n2
    x2_matrix(m, 1) = 0;
end

for n = 1:n1
    for m = n:(n + n1 - 1)
        x1_matrix(m, n) = x1(m - n + 1);
    end
end

for m = 1:n2
    x2_matrix(m, 1) = x2(m);
end

y_matrix(2:n1 + n2, 1) = x1_matrix * x2_matrix;

% stem(y_matrix, 'filled')
% legend
% hold off;
% figure;
% Supp 4. end


% nexttile
% stem(y_matrix, 'filled')
% title('y[n] = myconv(x1, x2)')
% figure;
% Supp 3. end

% Plot myconv(x1, x2)
stem(y_matrix, 'filled')
title('y[n] = myconv(x1, x2) = x1[n] * x2[n]')
