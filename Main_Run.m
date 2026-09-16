clear; 
clc;

c = 299792458;

f_start = 58.1e9;
f_stop = 63.1e9;
B = f_stop - f_start;

t_chirp = 100e-6;
S = B / t_chirp;

dt = 1e-12;

antenna_gain = 3;
if_gain = 700;

adc_fs = 2e6;

target_range = 0.6;
target_rcs = 1;

lambda = c / ((f_start + f_stop) / 2);

t = 0:dt:(t_chirp - dt);

% Keep LO phase fixed
lo_phase = 0;

% Four transmitted start phases
phase_set = [0, pi/2, pi, 3*pi/2];

% Nonlinearity parameters
a1 = 1;
a2 = 0.2;

% Generate fixed LO
lo = Chirp_Gen(t, f_start, S, lo_phase);

% ADC downsampling factor
downsample_factor = round(1 / (dt * adc_fs));

% Number of ADC samples
num_adc_samples = ceil(length(t) / downsample_factor);

% Store one ADC capture for each phase
ADC_all = zeros(length(phase_set), num_adc_samples);

for k = 1:length(phase_set)

    tx_phase = phase_set(k);

    % Transmitted chirp
    tx = Chirp_Gen(t, f_start, S, tx_phase);

    % Propagation
    rx = channel_propagation_model( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, tx_phase, target_range, target_rcs);

    % Mixer
    mixed = mixer(lo, rx);

    % IF amplifier + low-pass filter
    if_linear = If_Amp_LowPass_Filter(dt, mixed, if_gain);

    % Nonlinearity
    if_signal = Nonlinearity(if_linear, a1, a2);

    % ADC
    ADC = downsample(if_signal, downsample_factor);
    ADC = ADC - mean(ADC);

    ADC_all(k, :) = ADC;

end

% Plot spectrum before and after phase cycling
FFT_ADC(ADC_all, phase_set, adc_fs, c, S);