
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

target_A_range = 0.6;
target_A_rcs = 1;

target_B_range = 2 * target_A_range;
target_B_rcs = 0.001;

lambda = c / ((f_start + f_stop) / 2);

t = 0:dt:(t_chirp - dt);

lo_phase = 0;

% phase states
phase_set = [0, pi/2, pi, 3*pi/2];

a1 = 1;
a2 = 0.2;

lo = Chirp_Gen(t, f_start, S, lo_phase);

downsample_factor = round(1 / (dt * adc_fs));
num_adc_samples = ceil(length(t) / downsample_factor);

ADC_all = zeros(length(phase_set), num_adc_samples);
ADC_B_only = zeros(length(phase_set), num_adc_samples);
ADC_linear_all = zeros(length(phase_set), num_adc_samples);
ADC_linear_B = zeros(length(phase_set), num_adc_samples);

for k = 1:length(phase_set)

    tx_phase = phase_set(k);
    
    % a
    rx_A = channel_propagation_model( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, tx_phase, ...
        target_A_range, target_A_rcs);
    
    % b
    rx_B = channel_propagation_model( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, tx_phase, ...
        target_B_range, target_B_rcs);

    % a + b
    rx = rx_A + rx_B;

    mixed = mixer(lo, rx);

    if_linear = If_Amp_LowPass_Filter( ...
        dt, mixed, if_gain);

    if_signal = Nonlinearity( ...
        if_linear, a1, a2);

    ADC = downsample(if_signal, downsample_factor);
    ADC = ADC - mean(ADC);

    ADC_all(k, :) = ADC;

    % Linear (a2 = 0) reference for A + B
    ADC_ref = downsample(a1 * if_linear, downsample_factor);
    ADC_linear_all(k, :) = ADC_ref - mean(ADC_ref);

    % b only
    mixed_B = mixer(lo, rx_B);

    if_linear_B = If_Amp_LowPass_Filter( ...
        dt, mixed_B, if_gain);

    if_signal_B = Nonlinearity( ...
        if_linear_B, a1, a2);

    ADC_B = downsample( ...
        if_signal_B, downsample_factor);

    ADC_B = ADC_B - mean(ADC_B);

    ADC_B_only(k, :) = ADC_B;

    % Linear (a2 = 0) reference for target B
    ADC_B_ref = downsample(a1 * if_linear_B, downsample_factor);
    ADC_linear_B(k, :) = ADC_B_ref - mean(ADC_B_ref);

end

FFT_ADC( ...
    ADC_all, ADC_B_only, phase_set, ...
    adc_fs, c, S, ADC_linear_all, ADC_linear_B, ...
    target_A_range, target_B_range);
