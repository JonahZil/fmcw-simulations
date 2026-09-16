function chirp = Chirp_Gen(t, f_start, S, phase_start)

    phase = 2*pi*(f_start*t + 0.5*S*t.^2) + phase_start;
    chirp = sin(phase);

end