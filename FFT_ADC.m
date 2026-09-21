function FFT_ADC( ...
    ADC_all, ADC_B_only, ADC_B_ideal, ...
    phase_set, adc_fs, c, S)

num_phases = size(ADC_all, 1);
N = size(ADC_all, 2);

w = hann(N).';
w = w / mean(w);

Y_all = zeros(num_phases, N);
Y_B = zeros(num_phases, N);
Y_B_ideal = zeros(num_phases, N);

for k = 1:num_phases

    Y_all(k, :) = ...
        fft(ADC_all(k, :) .* w);

    Y_B(k, :) = ...
        fft(ADC_B_only(k, :) .* w);

    Y_B_ideal(k, :) = ...
        fft(ADC_B_ideal(k, :) .* w);

end


% =========================================================
% ONE-PHASE RESULTS
% =========================================================

Y_before = Y_all(1, :);

Y_B_before = Y_B(1, :);

Y_B_ideal_before = Y_B_ideal(1, :);


% =========================================================
% MULTI-PHASE RESULTS
% =========================================================

Y_after = zeros(1, N);

Y_B_after = zeros(1, N);

Y_B_ideal_after = zeros(1, N);

for k = 1:num_phases

    phi = phase_set(k);

    correction = exp(1j * phi);

    Y_after = Y_after ...
        + Y_all(k, :) * correction;

    Y_B_after = Y_B_after ...
        + Y_B(k, :) * correction;

    Y_B_ideal_after = Y_B_ideal_after ...
        + Y_B_ideal(k, :) * correction;

end

Y_after = Y_after / num_phases;

Y_B_after = Y_B_after / num_phases;

Y_B_ideal_after = ...
    Y_B_ideal_after / num_phases;


% =========================================================
% EXTRACT DAC-INDUCED ERROR
% =========================================================

E_B_before = ...
    Y_B_before - Y_B_ideal_before;

E_B_after = ...
    Y_B_after - Y_B_ideal_after;


% =========================================================
% POSITIVE FREQUENCIES
% =========================================================

num_positive = floor(N/2) + 1;

Y_before_mag = ...
    abs(Y_before(1:num_positive));

Y_after_mag = ...
    abs(Y_after(1:num_positive));

Y_B_before_mag = ...
    abs(Y_B_before(1:num_positive));

Y_B_after_mag = ...
    abs(Y_B_after(1:num_positive));

Y_B_ideal_before_mag = ...
    abs(Y_B_ideal_before(1:num_positive));

Y_B_ideal_after_mag = ...
    abs(Y_B_ideal_after(1:num_positive));

E_B_before_positive = ...
    E_B_before(1:num_positive);

E_B_after_positive = ...
    E_B_after(1:num_positive);

E_B_before_mag = ...
    abs(E_B_before_positive);

E_B_after_mag = ...
    abs(E_B_after_positive);


% =========================================================
% RANGE AXIS
% =========================================================

f = (0:num_positive-1) * adc_fs / N;

range_axis = c * f / (2*S);


% =========================================================
% ORIGINAL PLOT 1
% A + B HARMONIC CANCELLATION
% =========================================================

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(Y_before_mag(2:end) + eps), ...
    'LineWidth', 1.5);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(Y_after_mag(2:end) + eps), ...
    'LineWidth', 1.5);

plot( ...
    range_axis(2:end), ...
    20*log10(Y_B_after_mag(2:end) + eps), ...
    '--', ...
    'LineWidth', 1.5);

xlabel("Apparent Range (m)");
ylabel("Magnitude (dB)");

legend( ...
    "A + B before cancellation", ...
    sprintf("A + B after %d phases", num_phases), ...
    "B-only reference");

title("Recovery of Real Target Beneath H2");

grid on;
xlim([0 2]);


% =========================================================
% ORIGINAL PLOT 2
% SINGLE TARGET
% =========================================================

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(Y_B_before_mag(2:end) + eps), ...
    'LineWidth', 1.2);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(Y_B_after_mag(2:end) + eps), ...
    'LineWidth', 1.2);

xlabel("Apparent Range (m)");
ylabel("Magnitude (dB)");

legend( ...
    "1 phase", ...
    sprintf("%d phases", num_phases));

title(sprintf( ...
    "Single Target: 1 Phase vs %d Phases", ...
    num_phases));

grid on;
xlim([0 2]);


% =========================================================
% PLOT 3
% QUANTIZED DAC VS IDEAL DAC
% =========================================================

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(Y_B_before_mag(2:end) + eps), ...
    'LineWidth', 1.2);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(Y_B_ideal_before_mag(2:end) + eps), ...
    '--', ...
    'LineWidth', 1.2);

xlabel("Apparent Range (m)");
ylabel("Magnitude (dB)");

legend( ...
    "Quantized DAC", ...
    "Ideal DAC");

title("Single Target: Quantized DAC vs Ideal DAC");

grid on;
xlim([0 2]);


% =========================================================
% PLOT 4
% DAC ERROR FLOOR DIRECTLY
% =========================================================

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(E_B_before_mag(2:end) + eps), ...
    'LineWidth', 1.2);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(E_B_after_mag(2:end) + eps), ...
    'LineWidth', 1.2);

xlabel("Apparent Range (m)");
ylabel("DAC Error Magnitude (dB)");

legend( ...
    "1 phase", ...
    sprintf("%d phases", num_phases));

title(sprintf( ...
    "Extracted DAC Error: 1 Phase vs %d Phases", ...
    num_phases));

grid on;
xlim([0 2]);


% =========================================================
% QUANTITATIVE FLOOR COMPARISON
% =========================================================

analysis_mask = ...
    (range_axis > 0) & ...
    (range_axis <= 2);

ideal_search = Y_B_ideal_after_mag;
ideal_search(~analysis_mask) = 0;

[~, target_index] = max(ideal_search);

guard_bins = 2;

lower_index = ...
    max(1, target_index - guard_bins);

upper_index = ...
    min(num_positive, target_index + guard_bins);

floor_mask = analysis_mask;

floor_mask(lower_index:upper_index) = false;


% RMS DAC-error floor
rms_1_phase = sqrt( ...
    mean( ...
    abs(E_B_before_positive(floor_mask)).^2));

rms_multi_phase = sqrt( ...
    mean( ...
    abs(E_B_after_positive(floor_mask)).^2));

rms_reduction_db = ...
    20 * log10( ...
    (rms_1_phase + eps) / ...
    (rms_multi_phase + eps));


% Median DAC-error floor
median_1_phase = median( ...
    abs(E_B_before_positive(floor_mask)));

median_multi_phase = median( ...
    abs(E_B_after_positive(floor_mask)));

median_reduction_db = ...
    20 * log10( ...
    (median_1_phase + eps) / ...
    (median_multi_phase + eps));


% =========================================================
% PRINT RESULTS
% =========================================================

fprintf("\n");
fprintf("DAC ERROR ANALYSIS\n");
fprintf("------------------\n");

fprintf( ...
    "Number of phase states:      %d\n", ...
    num_phases);

fprintf("\n");

fprintf( ...
    "RMS error floor, 1 phase:    %.6e\n", ...
    rms_1_phase);

fprintf( ...
    "RMS error floor, %d phases:  %.6e\n", ...
    num_phases, rms_multi_phase);

fprintf( ...
    "RMS floor reduction:         %.2f dB\n", ...
    rms_reduction_db);

fprintf("\n");

fprintf( ...
    "Median error floor, 1 phase:    %.6e\n", ...
    median_1_phase);

fprintf( ...
    "Median error floor, %d phases:  %.6e\n", ...
    num_phases, median_multi_phase);

fprintf( ...
    "Median floor reduction:         %.2f dB\n", ...
    median_reduction_db);

fprintf("\n");

end