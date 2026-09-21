clear;
clc;
rng(1);

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

target_B_range = 1.2;
target_B_rcs = 0.01;

lambda = c / ((f_start + f_stop) / 2);

t = 0:dt:(t_chirp - dt);

phase_set = 2*pi*(0:15)/16;

phase_error_std_deg = 2;
phase_error = deg2rad(phase_error_std_deg) * randn(size(phase_set));

a1 = 1;
a2 = 0.2;

dac_step = 0.00578;

lo = Chirp_Gen(t, f_start, S, 0);

downsample_factor = round(1 / (dt * adc_fs));
num_adc_samples = ceil(length(t) / downsample_factor);

ADC_all = zeros(length(phase_set), num_adc_samples);
ADC_A = zeros(length(phase_set), num_adc_samples);
ADC_B = zeros(length(phase_set), num_adc_samples);
ADC_B_reference = zeros(length(phase_set), num_adc_samples);

for k = 1:length(phase_set)

    phi_actual = phase_set(k) + phase_error(k);

    rx_A = channel_propagation_model_DAC( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, phi_actual, ...
        target_A_range, target_A_rcs, dac_step);

    rx_B = channel_propagation_model_DAC( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, phi_actual, ...
        target_B_range, target_B_rcs, dac_step);

    rx_B_reference = channel_propagation_model_DAC( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, phase_set(k), ...
        target_B_range, target_B_rcs, dac_step);

    mixed = mixer(lo, rx_A + rx_B);
    if_linear = If_Amp_LowPass_Filter(dt, mixed, if_gain);
    if_signal = Nonlinearity(if_linear, a1, a2);

    ADC = downsample(if_signal, downsample_factor);
    ADC_all(k, :) = ADC - mean(ADC);

    mixed = mixer(lo, rx_A);
    if_linear = If_Amp_LowPass_Filter(dt, mixed, if_gain);
    if_signal = Nonlinearity(if_linear, a1, a2);

    ADC = downsample(if_signal, downsample_factor);
    ADC_A(k, :) = ADC - mean(ADC);

    mixed = mixer(lo, rx_B);
    if_linear = If_Amp_LowPass_Filter(dt, mixed, if_gain);
    if_signal = Nonlinearity(if_linear, a1, a2);

    ADC = downsample(if_signal, downsample_factor);
    ADC_B(k, :) = ADC - mean(ADC);

    mixed = mixer(lo, rx_B_reference);
    if_linear = If_Amp_LowPass_Filter(dt, mixed, if_gain);
    if_signal = Nonlinearity(if_linear, a1, a2);

    ADC = downsample(if_signal, downsample_factor);
    ADC_B_reference(k, :) = ADC - mean(ADC);

end

FFT_ADC( ...
    ADC_all, ADC_A, ADC_B, ADC_B_reference, ...
    phase_set, adc_fs, c, S, ...
    target_B_range, phase_error_std_deg);