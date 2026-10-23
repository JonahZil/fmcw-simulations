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
ADC_ref = zeros(length(phase_set), num_adc_samples);

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

    % Nonlinear signal
    if_nonlinear = Nonlinearity(if_linear, a1, a2);
    
    % Linear reference (a2 = 0)
    if_reference = Nonlinearity(if_linear, a1, 0);
    
    % ADC sampling
    ADC_nonlinear = downsample(if_nonlinear, downsample_factor);
    ADC_reference = downsample(if_reference, downsample_factor);
    
    % DC removal
    ADC_nonlinear = ADC_nonlinear - mean(ADC_nonlinear);
    ADC_reference = ADC_reference - mean(ADC_reference);
    
    % Store both signals
    ADC_all(k, :) = ADC_nonlinear;
    ADC_ref(k, :) = ADC_reference;

end

% Plot spectrum before and after phase cycling
FFT_ADC(ADC_all, ADC_ref, phase_set, adc_fs, c, S, target_range);