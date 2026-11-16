clear; clc; close all;

% Parameters
B = 12;
K = 8;
A = 0.9;
Fs = 100e6;
Tc = 100e-6;
f0 = 1e6;
f1 = 2e6;

t = (0:round(Fs*Tc)-1)'/Fs;
N = length(t);
S = (f1-f0)/Tc;
phi = 2*pi*(0:K-1)/K;
step = 2/2^B;

% FMCW chirp phase
psi = 2*pi*(f0*t + 0.5*S*t.^2);

% Generate DAC quantization errors
q = zeros(N,K);

for k = 1:K
    x = A*sin(psi + phi(k));
    xq = step*round(x/step);
    q(:,k) = xq-x;
end

% Phase compensation and coherent averaging
qAvg = real(mean(hilbert(q).*exp(-1j*phi),2));

% Metrics
rms1 = sqrt(mean(q(:,1).^2));
rmsK = sqrt(mean(qAvg.^2));

suppression = 20*log10(rms1/rmsK);

sqnr1 = 20*log10((A/sqrt(2))/rms1);
sqnrK = 20*log10((A/sqrt(2))/rmsK);

fprintf('RMS error (1 phase): %.4f LSB\n',rms1/step);
fprintf('RMS error (8 phases): %.4f LSB\n',rmsK/step);
fprintf('Error reduction: %.2f dB\n',suppression);
fprintf('SQNR: %.2f -> %.2f dB\n',sqnr1,sqnrK);

% Quantization error FFT
f = (0:N/2)'*Fs/N;

E1 = fft(q(:,1));
EK = fft(qAvg);

dB1 = 20*log10(max(2*abs(E1(1:N/2+1))/(N*A),1e-15));
dBK = 20*log10(max(2*abs(EK(1:N/2+1))/(N*A),1e-15));

% Figure
figure('Color','w','Position',[100 100 1100 420]);
title('Quantization Error before and After Phase-Diversity')

idx = t >= 25e-6 & t <= 25.4e-6;
plot(t(idx)*1e6,q(idx,1)/step,'o-','LineWidth',1.4);
hold on;
plot(t(idx)*1e6,qAvg(idx)/step,'s-','LineWidth',1.4);
grid on;
xlabel('Time (us)');
ylabel('Quantization error (LSB)');
legend('Single phase','8 phases');

figure;

plot(f/1e6,dB1,'LineWidth',1);
title('FFT of Quantization Error Before and After Phase-Diversity')
hold on;
plot(f/1e6,dBK,'LineWidth',1);
grid on;
xlim([0 12]);
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB relative to A/bin)');
legend('Single phase','8 phases');

