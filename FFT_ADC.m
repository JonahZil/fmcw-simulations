function FFT_ADC( ...
    ADC_clean, ADC_noisy, ...
    phase_set, adc_fs, c, S)

num_phases = length(phase_set);
N = size(ADC_clean, 2);

w = hann(N).';
w = w / mean(w);

Y_clean = zeros(num_phases, N);
Y_noisy = zeros(num_phases, N);

for k = 1:num_phases

    Y_clean(k, :) = ...
        fft(ADC_clean(k, :) .* w);

    Y_noisy(k, :) = ...
        fft(ADC_noisy(k, :) .* w);

end

Y_clean_16 = zeros(1, N);
Y_noisy_16 = zeros(1, N);

for k = 1:num_phases

    correction = exp(1j * phase_set(k));

    Y_clean_16 = ...
        Y_clean_16 + Y_clean(k, :) * correction;

    Y_noisy_16 = ...
        Y_noisy_16 + Y_noisy(k, :) * correction;

end

Y_clean_16 = Y_clean_16 / num_phases;
Y_noisy_16 = Y_noisy_16 / num_phases;

Y_noisy_1 = Y_noisy(1, :);

num_positive = floor(N/2) + 1;

Y_clean_16 = Y_clean_16(1:num_positive);
Y_noisy_16 = Y_noisy_16(1:num_positive);
Y_noisy_1 = Y_noisy_1(1:num_positive);

f = (0:num_positive-1) * adc_fs / N;
range_axis = c * f / (2*S);

analysis_mask = ...
    (range_axis > 0) & ...
    (range_axis <= 2);

clean_search = abs(Y_clean_16);
clean_search(~analysis_mask) = 0;

[target_peak, target_index] = max(clean_search);

guard_bins = 4;

lower_index = max(1, target_index - guard_bins);
upper_index = min(num_positive, target_index + guard_bins);

floor_mask = analysis_mask;
floor_mask(lower_index:upper_index) = false;

clean_floor = sqrt( ...
    mean(abs(Y_clean_16(floor_mask)).^2));

noisy_1_floor = sqrt( ...
    mean(abs(Y_noisy_1(floor_mask)).^2));

noisy_16_floor = sqrt( ...
    mean(abs(Y_noisy_16(floor_mask)).^2));

clean_floor_db = ...
    20*log10(clean_floor / target_peak);

noisy_1_floor_db = ...
    20*log10(noisy_1_floor / target_peak);

noisy_16_floor_db = ...
    20*log10(noisy_16_floor / target_peak);

phase_reduction_db = ...
    noisy_1_floor_db - noisy_16_floor_db;

thermal_increase_db = ...
    noisy_16_floor_db - clean_floor_db;

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_clean_16(2:end)) / target_peak + eps), ...
    'LineWidth', 1.2);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_noisy_1(2:end)) / target_peak + eps), ...
    'LineWidth', 1.2);

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_noisy_16(2:end)) / target_peak + eps), ...
    'LineWidth', 1.2);

xlabel("Apparent Range (m)");
ylabel("Magnitude Relative to Target (dB)");

legend( ...
    "No thermal noise, 16 phases", ...
    "Thermal noise, 1 phase", ...
    "Thermal noise, 16 phases");

title("Thermal Noise: 1 Phase vs 16 Phases");

grid on;
xlim([0 2]);

fprintf("\n");
fprintf("THERMAL NOISE ANALYSIS\n");
fprintf("----------------------\n");

fprintf( ...
    "Clean 16-phase floor:    %.2f dB\n", ...
    clean_floor_db);

fprintf( ...
    "Noisy 1-phase floor:     %.2f dB\n", ...
    noisy_1_floor_db);

fprintf( ...
    "Noisy 16-phase floor:    %.2f dB\n", ...
    noisy_16_floor_db);

fprintf( ...
    "1 to 16 phase reduction: %.2f dB\n", ...
    phase_reduction_db);

fprintf( ...
    "Thermal-noise increase:  %.2f dB\n", ...
    thermal_increase_db);

end