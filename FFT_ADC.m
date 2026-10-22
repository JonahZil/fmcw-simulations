function FFT_ADC(ADC_all, phase_set, adc_fs, c, S, target_range)
    num_phases = size(ADC_all, 1);
    N = size(ADC_all, 2);

    w = hann(N).';
    w = w / mean(w);

    % Store complex FFT for every transmitted phase
    Y_all = zeros(num_phases, N);

    for k = 1:num_phases

        Y_all(k, :) = fft(ADC_all(k, :) .* w);

    end

    % ========================================================
    % Before cancellation
    % ========================================================

    % Use phase = 0 capture as the original spectrum
    Y_before = Y_all(1, :);

    % ========================================================
    % Digital phase demodulation
    % ========================================================

    Y_after = zeros(1, N);

    for k = 1:num_phases

        phi = phase_set(k);

        % Fundamental rotates approximately as exp(-j*phi).
        %
        % Multiplying by exp(+j*phi) puts the
        % fundamental from every capture back in phase.
        %
        % The second harmonic rotates as exp(-j*2*phi),
        % so after correction it still rotates as
        % exp(-j*phi).

        Y_demod = Y_all(k, :) * exp(1j * phi);

        Y_after = Y_after + Y_demod;

    end

    % Coherent average
    Y_after = Y_after / num_phases;

    % ========================================================
    % Positive frequencies
    % ========================================================

    num_positive = floor(N/2) + 1;

    Y_before = abs(Y_before(1:num_positive));
    Y_after = abs(Y_after(1:num_positive));

    f = (0:num_positive-1) * adc_fs / N;

    range_axis = c * f / (2*S);
    
    % ========================================================
    % Positive frequencies
    % ========================================================

    num_positive = floor(N/2) + 1;

    % Preserve complex FFT values for metric calculations
    Y_before = Y_before(1:num_positive);
    Y_after = Y_after(1:num_positive);

    f = (0:num_positive-1) * adc_fs / N;
    range_axis = c * f / (2*S);

    % ========================================================
    % Quantitative performance metrics
    % ========================================================

    % Fundamental and second-harmonic apparent ranges
    R_H1 = target_range;
    R_H2 = 2 * target_range;

    % Find nearest FFT bins
    [~, idx_H1] = min(abs(range_axis - R_H1));
    [~, idx_H2] = min(abs(range_axis - R_H2));

    % Complex FFT coefficients
    H1_before = Y_before(idx_H1);
    H1_after = Y_after(idx_H1);

    H2_before = Y_before(idx_H2);
    H2_after = Y_after(idx_H2);

    % 1. Harmonic suppression
    H2_suppression_db = 20*log10( ...
        (abs(H2_before) + eps) / ...
        (abs(H2_after) + eps));

    % 2. Fundamental amplitude retention
    H1_retention_db = 20*log10( ...
        (abs(H1_after) + eps) / ...
        (abs(H1_before) + eps));

    % 3. Harmonic relative to fundamental
    H2_relative_before_db = 20*log10( ...
        (abs(H2_before) + eps) / ...
        (abs(H1_before) + eps));

    H2_relative_after_db = 20*log10( ...
        (abs(H2_after) + eps) / ...
        (abs(H1_after) + eps));

    % 4. Fundamental phase error
    H1_phase_error_deg = rad2deg( ...
        angle(H1_after * conj(H1_before)));

    % ========================================================
    % Print results
    % ========================================================

    fprintf('\n');
    fprintf('PHASE-DIVERSITY PERFORMANCE\n');
    fprintf('---------------------------\n');

    fprintf('H1 range bin:          %.4f m\n', range_axis(idx_H1));
    fprintf('H2 range bin:          %.4f m\n', range_axis(idx_H2));

    fprintf('H2 suppression:        %.2f dB\n', H2_suppression_db);
    fprintf('H1 amplitude change:   %.4f dB\n', H1_retention_db);

    fprintf('H2/H1 before:          %.2f dBc\n', H2_relative_before_db);
    fprintf('H2/H1 after:           %.2f dBc\n', H2_relative_after_db);

    fprintf('H1 phase error:        %.4f degrees\n', H1_phase_error_deg);

    % Convert to magnitude for plotting
    Y_before = abs(Y_before);
    Y_after = abs(Y_after);

    % ========================================================
    % Plot
    % ========================================================

    figure;

    plot( ...
        range_axis(2:end), ...
        20*log10(Y_before(2:end) + eps), ...
        'LineWidth', 1.5);

    hold on;

    plot( ...
        range_axis(2:end), ...
        20*log10(Y_after(2:end) + eps), ...
        'LineWidth', 1.5);

    xlabel("Apparent Range (m)");
    ylabel("Magnitude (dB)");

    legend( ...
        "Before phase cycling", ...
        "After phase cycling");

    title("Harmonic Cancellation Using Start-Phase Modulation");

    grid on;
    xlim([0 2]);

end