
clear; clc; close all;

% Parameters
B = 12;
K = 8;
N = 2^16;
A = 0.9;

theta = 2*pi*(0:N-1)'/N;
phi = 2*pi*(0:K-1)/K;
step = 2/2^B;

% Generate quantization errors
q = zeros(N,K);

for k = 1:K
    x = A*sin(theta + phi(k));
    xq = step*round(x/step);
    q(:,k) = xq - x;
end

% Fourier coefficients of the quantization error
Q = fft(q)/N;

% Phase compensation and averaging
Qavg = mean(Q .* exp(-1j*phi), 2);

% Plot harmonic amplitudes
m = 1:2:19;

before = 20*log10(2*abs(Q(m+1,1))/A);
after = 20*log10(2*abs(Qavg(m+1))/A);

plot(m,before,'o-',m,after,'x-');
grid on;
xlabel('Harmonic order');
ylabel('Quantization error (dBc)');
legend('Single phase','Phase diversity');
xticks(m);

% Reconstruct the averaged time-domain error
Qreal = zeros(N,1);
Qreal(2:N/2) = Qavg(2:N/2);
Qreal(N/2+2:end) = conj(flipud(Qreal(2:N/2)));

q4 = real(ifft(Qreal*N));

% Compare DAC outputs
x1 = A*sin(theta) + q(:,1);
x4 = A*sin(theta) + q4;


% Plot quantization error
figure('Color','w','Position',[100 100 900 400]);

idx = 1:80;

plot(theta(idx),q(idx,1)/step,'-o', ...
    'LineWidth',1.5,'MarkerSize',3);
hold on;

plot(theta(idx),q4(idx)/step,'-s', ...
    'LineWidth',1.5,'MarkerSize',3);

grid on;
xlabel('Phase (rad)');
ylabel('Quantization error (LSB)');
title('Quantization error');
legend('One phase','8 phases');
ylim([-0.6 0.6]);
