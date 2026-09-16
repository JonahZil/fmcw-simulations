function FFT_ADC(ADC_all, ADC_B_only, phase_set, adc_fs, c, S)

    num_phases = size(ADC_all, 1);
    N = size(ADC_all, 2);

    w = hann(N).';
    w = w / mean(w);

    Y_all = zeros(num_phases, N);
    Y_B = zeros(num_phases, N);

    for k = 1:num_phases

        Y_all(k, :) = ...
            fft(ADC_all(k, :) .* w);

        Y_B(k, :) = ...
            fft(ADC_B_only(k, :) .* w);

    end
    
    % before demodulation
    Y_before = Y_all(1, :);
    
    % after demodulation
    Y_after = zeros(1, N);

    % b reference
    Y_B_after = zeros(1, N);

    for k = 1:num_phases

        phi = phase_set(k);

        correction = exp(1j * phi);

        Y_after = Y_after ...
            + Y_all(k, :) * correction;

        Y_B_after = Y_B_after ...
            + Y_B(k, :) * correction;

    end

    Y_after = Y_after / num_phases;
    Y_B_after = Y_B_after / num_phases;
    
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