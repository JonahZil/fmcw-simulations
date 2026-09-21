function FFT_ADC( ...
    ADC_all, ADC_A, ADC_B, ADC_B_reference, ...
    phase_set, adc_fs, c, S, ...
    target_B_range, phase_error_std_deg)

num_phases = length(phase_set);
N = size(ADC_all, 2);

w = hann(N).';
w = w / mean(w);

Y_all = zeros(num_phases, N);
Y_A = zeros(num_phases, N);
Y_B = zeros(num_phases, N);
Y_B_reference = zeros(num_phases, N);

for k = 1:num_phases

    Y_all(k, :) = fft(ADC_all(k, :) .* w);
    Y_A(k, :) = fft(ADC_A(k, :) .* w);
    Y_B(k, :) = fft(ADC_B(k, :) .* w);
    Y_B_reference(k, :) = fft(ADC_B_reference(k, :) .* w);

end

Y_A_before = Y_A(1, :);

Y_all_after = zeros(1, N);
Y_A_after = zeros(1, N);
Y_B_after = zeros(1, N);
Y_B_reference_after = zeros(1, N);

for k = 1:num_phases

    correction = exp(1j * phase_set(k));

    Y_all_after = ...
        Y_all_after + Y_all(k, :) * correction;

    Y_A_after = ...
        Y_A_after + Y_A(k, :) * correction;

    Y_B_after = ...
        Y_B_after + Y_B(k, :) * correction;

    Y_B_reference_after = ...
        Y_B_reference_after ...
        + Y_B_reference(k, :) * correction;

end

Y_all_after = Y_all_after / num_phases;
Y_A_after = Y_A_after / num_phases;
Y_B_after = Y_B_after / num_phases;
Y_B_reference_after = Y_B_reference_after / num_phases;

num_positive = floor(N/2) + 1;

Y_A_before = Y_A_before(1:num_positive);
Y_all_after = Y_all_after(1:num_positive);
Y_A_after = Y_A_after(1:num_positive);
Y_B_after = Y_B_after(1:num_positive);
Y_B_reference_after = Y_B_reference_after(1:num_positive);

f = (0:num_positive-1) * adc_fs / N;
range_axis = c * f / (2*S);

[~, target_bin] = ...
    min(abs(range_axis - target_B_range));

H2_before = abs(Y_A_before(target_bin));
H2_after = abs(Y_A_after(target_bin));

harmonic_suppression_db = ...
    20*log10((H2_before + eps) / (H2_after + eps));

B_phase_error = abs(Y_B_after(target_bin));
B_reference = abs(Y_B_reference_after(target_bin));

target_retention_db = ...
    20*log10((B_phase_error + eps) / (B_reference + eps));

AB_recovered = abs(Y_all_after(target_bin));

recovery_error_db = ...
    20*log10((AB_recovered + eps) / (B_reference + eps));

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_A_before(2:end)) + eps), ...
    'LineWidth', 1.2);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_A_after(2:end)) + eps), ...
    'LineWidth', 1.2);

xlabel("Apparent Range (m)");
ylabel("Magnitude (dB)");

legend( ...
    "Harmonic before cancellation", ...
    "Harmonic after cancellation");

title(sprintf( ...
    "Harmonic Suppression with %.1f Degree Phase Error", ...
    phase_error_std_deg));

grid on;
xlim([0 2]);

figure;

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_all_after(2:end)) + eps), ...
    'LineWidth', 1.2);

hold on;

plot( ...
    range_axis(2:end), ...
    20*log10(abs(Y_B_reference_after(2:end)) + eps), ...
    '--', ...
    'LineWidth', 1.2);

xlabel("Apparent Range (m)");
ylabel("Magnitude (dB)");

legend( ...
    "A + B after cancellation", ...
    "Ideal-phase B reference");

title(sprintf( ...
    "Target Recovery with %.1f Degree Phase Error", ...
    phase_error_std_deg));

grid on;
xlim([0 2]);

fprintf("\n");
fprintf("PHASE ERROR ANALYSIS\n");
fprintf("--------------------\n");

fprintf( ...
    "Phase error std:       %.2f degrees\n", ...
    phase_error_std_deg);

fprintf( ...
    "H2 suppression:        %.2f dB\n", ...
    harmonic_suppression_db);

fprintf( ...
    "Real target retention: %.2f dB\n", ...
    target_retention_db);

fprintf( ...
    "Recovered target error: %.2f dB\n", ...
    recovery_error_db);

end