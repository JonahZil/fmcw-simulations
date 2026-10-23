
function FFT_ADC(ADC_all, ADC_ref, phase_set, adc_fs, c, S, target_range)

    num_phases = size(ADC_all, 1);
    N = size(ADC_all, 2);

    % Hann window
    w = hann(N).';
    w = w / mean(w);

    % FFT of nonlinear and linear-reference signals
    Y_all = zeros(num_phases, N);
    Y_ref = zeros(num_phases, N);

    for k = 1:num_phases
        Y_all(k, :) = fft(ADC_all(k, :) .* w);
        Y_ref(k, :) = fft(ADC_ref(k, :) .* w);
    end

    % ========================================================
    % Before phase diversity
    % ========================================================

    Y_before = Y_all(1, :);
    Y_ref_before = Y_ref(1, :);

    % ========================================================
    % Phase compensation and coherent averaging
    % ========================================================

    Y_after = zeros(1, N);
    Y_ref_after = zeros(1, N);

    for k = 1:num_phases

        correction = exp(1j * phase_set(k));

        Y_after = Y_after + Y_all(k, :) * correction;
        Y_ref_after = Y_ref_after + Y_ref(k, :) * correction;

    end

    Y_after = Y_after / num_phases;
    Y_ref_after = Y_ref_after / num_phases;

    % ========================================================
    % Isolate nonlinear distortion
    % ========================================================

    D_before = Y_before - Y_ref_before;
    D_after = Y_after - Y_ref_after;

    % ========================================================
    % Positive-frequency spectrum
    % ========================================================

    num_positive = floor(N/2) + 1;

    Y_before = Y_before(1:num_positive);
    Y_after = Y_after(1:num_positive);

    Y_ref_before = Y_ref_before(1:num_positive);
    Y_ref_after = Y_ref_after(1:num_positive);

    D_before = D_before(1:num_positive);
    D_after = D_after(1:num_positive);

    f = (0:num_positive-1) * adc_fs / N;
    range_axis = c * f / (2*S);

    % ========================================================
    % Locate fundamental and second harmonic
    % ========================================================

    [~, idx_H1] = min(abs(range_axis - target_range));
    [~, idx_H2] = min(abs(range_axis - 2*target_range));

    % ========================================================
    % Quantitative metrics
    % ========================================================

    % Second-harmonic suppression (isolated distortion)
    H2_suppression_db = 20*log10( ...
        (abs(D_before(idx_H2)) + eps) / ...
        (abs(D_after(idx_H2)) + eps));

    % Fundamental amplitude retention
    H1_retention_db = 20*log10( ...
        (abs(Y_after(idx_H1)) + eps) / ...
        (abs(Y_before(idx_H1)) + eps));

    % Fundamental phase change
    H1_phase_error_deg = rad2deg( ...
        angle(Y_after(idx_H1) * conj(Y_before(idx_H1))));

    % H2 distortion relative to H1
    H2_relative_before_db = 20*log10( ...
        (abs(D_before(idx_H2)) + eps) / ...
        (abs(Y_ref_before(idx_H1)) + eps));

    H2_relative_after_db = 20*log10( ...
        (abs(D_after(idx_H2)) + eps) / ...
        (abs(Y_ref_after(idx_H1)) + eps));

    % ========================================================
    % Print results
    % ========================================================

    fprintf('\n');
    fprintf('PHASE-DIVERSITY PERFORMANCE\n');
    fprintf('---------------------------\n');

    fprintf('H1 range bin:         %.4f m\n', range_axis(idx_H1));
    fprintf('H2 range bin:         %.4f m\n', range_axis(idx_H2));

    fprintf('H2 suppression:       %.2f dB\n', H2_suppression_db);
    fprintf('H1 amplitude change:  %.4f dB\n', H1_retention_db);
    fprintf('H1 phase change:      %.4f degrees\n', H1_phase_error_deg);

    fprintf('H2/H1 before:         %.2f dBc\n', H2_relative_before_db);
    fprintf('H2/H1 after:          %.2f dBc\n', H2_relative_after_db);

    % ========================================================
    % Original spectrum plot
    % ========================================================

    figure;

    plot(range_axis(2:end), ...
        20*log10(abs(Y_before(2:end)) + eps), ...
        'LineWidth', 1.5);

    hold on;

    plot(range_axis(2:end), ...
        20*log10(abs(Y_after(2:end)) + eps), ...
        'LineWidth', 1.5);

    xlabel("Apparent Range (m)");
    ylabel("Magnitude (dB)");

    legend("Before phase cycling", ...
           "After phase cycling");

    title("Harmonic Cancellation Using Start-Phase Modulation");

    grid on;
    xlim([0 2]);

    % ========================================================
    % Isolated nonlinear distortion plot
    % ========================================================

    figure;

    plot(range_axis(2:end), ...
        20*log10(abs(D_before(2:end)) + eps), ...
        'LineWidth', 1.5);

    hold on;

    plot(range_axis(2:end), ...
        20*log10(abs(D_after(2:end)) + eps), ...
        'LineWidth', 1.5);

    xlabel("Apparent Range (m)");
    ylabel("Distortion Magnitude (dB)");

    legend("Distortion before cancellation", ...
           "Distortion after cancellation");

    title("Isolated Nonlinear Distortion");

    grid on;
    xlim([0 2]);

end
