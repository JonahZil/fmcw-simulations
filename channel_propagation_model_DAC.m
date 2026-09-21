function rx = channel_propagation_model_DAC( ...
    t, c, lambda, antenna_gain, ...
    f_start, S, phase_start, range, rcs, dac_step)

    delay = 2*range/c;

    amplitude = sqrt( ...
        (rcs * lambda^2 * antenna_gain^2) / ...
        (4*pi^3*range^4));

    rx = zeros(size(t));

    valid = t >= delay;
    delayed_t = t(valid) - delay;

    tx = Chirp_Gen( ...
        delayed_t, f_start, S, phase_start);

    tx = DAC_Model(tx, dac_step);

    rx(valid) = amplitude * tx;

end