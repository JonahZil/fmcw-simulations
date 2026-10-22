
function FFT_ADC(ADC_all, ADC_B_only, phase_set, adc_fs, c, S, ...
    ADC_linear_all, ADC_linear_B, target_A_range, target_B_range)

    num_phases = size(ADC_all, 1);
    N = size(ADC_all, 2);

    w = hann(N).';
    w = w / mean(w);

    Y_all = zeros(num_phases, N);
    Y_B = zeros(num_phases, N);
    Y_linear_all = zeros(num_phases, N);
    Y_linear_B = zeros(num_phases, N);

    for k = 1:num_phases

        Y_all(k, :) = ...
            fft(ADC_all(k, :) .* w);

        Y_B(k, :) = ...
            fft(ADC_B_only(k, :) .* w);

        Y_linear_all(k, :) = ...
            fft(ADC_linear_all(k, :) .* w);

        Y_linear_B(k, :) = ...
            fft(ADC_linear_B(k, :) .* w);

    end
    
    % before demodulation
    Y_before = Y_all(1, :);
    
    % after demodulation
    Y_after = zeros(1, N);

    % b reference
    Y_B_after = zeros(1, N);
    Y_linear_after = zeros(1, N);
    Y_linear_B_after = zeros(1, N);

    for k = 1:num_phases

        phi = phase_set(k);

        correction = exp(1j * phi);

        Y_after = Y_after ...
            + Y_all(k, :) * correction;

        Y_B_after = Y_B_after ...
            + Y_B(k, :) * correction;

        Y_linear_after = Y_linear_after ...
            + Y_linear_all(k, :) * correction;

        Y_linear_B_after = Y_linear_B_after ...
            + Y_linear_B(k, :) * correction;

    end

    Y_after = Y_after / num_phases;
    Y_B_after = Y_B_after / num_phases;
    Y_linear_after = Y_linear_after / num_phases;
    Y_linear_B_after = Y_linear_B_after / num_phases;

    % Quantitative measurements at target B (the H2 range of A).
    % Complex spectra must be used here, before taking magnitudes.
    range_resolution = c * adc_fs / (2*S*N);
    bin_A = round(target_A_range / range_resolution) + 1;
    bin_B = round(target_B_range / range_resolution) + 1;

    % Nonlinear distortion, with the linear A+B reference subtracted.
    distortion_before = Y_all(1, bin_B) - Y_linear_all(1, bin_B);
    distortion_after = Y_after(bin_B) - Y_linear_after(bin_B);
    h2_suppression_db = 20*log10( ...
        (abs(distortion_before) + eps) / (abs(distortion_after) + eps));

    % Subtract linear target A to isolate the recovered B contribution.
    linear_A_before = Y_linear_all(1, bin_B) - Y_linear_B(1, bin_B);
    linear_A_after = Y_linear_after(bin_B) - Y_linear_B_after(bin_B);
    B_recovered_before = Y_all(1, bin_B) - linear_A_before;
    B_recovered_after = Y_after(bin_B) - linear_A_after;
    B_reference_before = Y_linear_B(1, bin_B);
    B_reference_after = Y_linear_B_after(bin_B);

    recovery_error_before_db = 20*log10( ...
        (abs(B_recovered_before - B_reference_before) + eps) / ...
        (abs(B_reference_before) + eps));
    recovery_error_after_db = 20*log10( ...
        (abs(B_recovered_after - B_reference_after) + eps) / ...
        (abs(B_reference_after) + eps));
    magnitude_error_before_db = 20*log10( ...
        (abs(B_recovered_before) + eps) / (abs(B_reference_before) + eps));
    magnitude_error_after_db = 20*log10( ...
        (abs(B_recovered_after) + eps) / (abs(B_reference_after) + eps));

    fprintf("Target A: %.3f m (FFT bin %d)\n", target_A_range, bin_A);
    fprintf("Target B: %.3f m (FFT bin %d, %.3f m apparent)\n", ...
        target_B_range, bin_B, (bin_B-1)*range_resolution);
    fprintf("H2 distortion at B: before %.3e, after %.3e\n", ...
        abs(distortion_before), abs(distortion_after));
    fprintf("H2 suppression at B bin: %.2f dB\n", h2_suppression_db);
    fprintf("B complex recovery error: before %.2f dB, after %.2f dB (relative to B)\n", ...
        recovery_error_before_db, recovery_error_after_db);
    fprintf("B magnitude bias: before %+.2f dB, after %+.2f dB\n", ...
        magnitude_error_before_db, magnitude_error_after_db);
    
    num_positive = floor(N/2) + 1;

    Y_before = abs(Y_before(1:num_positive));
    Y_after = abs(Y_after(1:num_positive));
    Y_B_after = abs(Y_B_after(1:num_positive));

    f = (0:num_positive-1) * adc_fs / N;

    range_axis = c * f / (2*S);
    
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

    plot( ...
        range_axis(2:end), ...
        20*log10(Y_B_after(2:end) + eps), ...
        '--', ...
        'LineWidth', 1.5);

    xlabel("Apparent Range (m)");
    ylabel("Magnitude (dB)");

    legend( ...
        "A + B before cancellation", ...
        "A + B after cancellation", ...
        "B-only reference");

    title("Recovery of Real Target Beneath H2");

    grid on;
    xlim([0 2]);

end
