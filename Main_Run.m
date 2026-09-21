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

target_range = 1.2;
target_rcs = 0.01;

lambda = c / ((f_start + f_stop) / 2);

t = 0:dt:(t_chirp - dt);

lo_phase = 0;

phase_set = 2*pi*(0:15)/16;

a1 = 1;
a2 = 0.2;

dac_step = 0.00578;

k_B = 1.380649e-23;
temperature = 290;
resistance = 50;
noise_bw = 800e3;

thermal_noise_rms = ...
    if_gain * sqrt(k_B * temperature * resistance * noise_bw);

lo = Chirp_Gen(t, f_start, S, lo_phase);

downsample_factor = round(1 / (dt * adc_fs));
num_adc_samples = ceil(length(t) / downsample_factor);

ADC_clean = zeros(length(phase_set), num_adc_samples);
ADC_noisy = zeros(length(phase_set), num_adc_samples);

for k = 1:length(phase_set)

    tx_phase = phase_set(k);

    rx = channel_propagation_model_DAC( ...
        t, c, lambda, antenna_gain, ...
        f_start, S, tx_phase, ...
        target_range, target_rcs, dac_step);

    mixed = mixer(lo, rx);

    if_linear = If_Amp_LowPass_Filter( ...
        dt, mixed, if_gain);

    if_signal = Nonlinearity( ...
        if_linear, a1, a2);

    ADC = downsample( ...
        if_signal, downsample_factor);

    ADC = ADC - mean(ADC);

    ADC_clean(k, :) = ADC;

    noise = Thermal_Noise( ...
        length(ADC), adc_fs, ...
        noise_bw, thermal_noise_rms);

    ADC_noisy(k, :) = ADC + noise;
    ADC_noisy(k, :) = ...
        ADC_noisy(k, :) - mean(ADC_noisy(k, :));

end

FFT_ADC( ...
    ADC_clean, ADC_noisy, ...
    phase_set, adc_fs, c, S);