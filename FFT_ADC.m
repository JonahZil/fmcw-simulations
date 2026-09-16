function FFT_ADC(ADC_all, phase_set, adc_fs, c, S)

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