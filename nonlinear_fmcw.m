chirp_times = (100:10:200) * 1e-6;

az_A = deg2rad(20);
el_A = deg2rad(0);

uAx = cos(el_A) * sin(az_A);
uAy = sin(el_A);

delta_r_A = rx_pos(:,1) * uAx + ...
            rx_pos(:,2) * uAy;

tau_A = (2 * range_A - delta_r_A) / c;

H2_A = zeros(length(chirp_times), num_rx);

for k = 1:length(chirp_times)

    t_chirp = chirp_times(k);
    S = B / t_chirp;

    t = (0:N-1) / f_sampling;

    tx = exp(1j * ...
        (2*pi*f_start*t + pi*S*t.^2));

    for m = 1:num_rx

        tau = tau_A(m);

        rx = exp(1j * ...
            (2*pi*f_start*(t - tau) + ...
            pi*S*(t - tau).^2));

        beat = tx .* conj(rx);

        h2 = a2 * beat.^2;

        f_h2 = 2 * S * tau;

        H2_A(k,m) = sum(h2 .* ...
            exp(-1j * 2*pi*f_h2*t));

    end
end

% Normalize H2 phase
H2_A_phase = H2_A ./ abs(H2_A);
h2_vec = H2_A_phase(:);

% Target B: real reflector at twice the range
range_B = 2 * range_A;

azimuths = -60:1:60;
elevations = -60:1:60;

rho_map = zeros(length(elevations), length(azimuths));

for i_el = 1:length(elevations)

    el_B = deg2rad(elevations(i_el));

    for i_az = 1:length(azimuths)

        az_B = deg2rad(azimuths(i_az));

        uBx = cos(el_B) * sin(az_B);
        uBy = sin(el_B);

        delta_r_B = rx_pos(:,1) * uBx + ...
                    rx_pos(:,2) * uBy;

        tau_B = (2 * range_B - delta_r_B) / c;

        H1_B = zeros(length(chirp_times), num_rx);

        for k = 1:length(chirp_times)

            t_chirp = chirp_times(k);
            S = B / t_chirp;

            t = (0:N-1) / f_sampling;

            tx = exp(1j * ...
                (2*pi*f_start*t + pi*S*t.^2));

            for m = 1:num_rx

                tau = tau_B(m);

                rx = exp(1j * ...
                    (2*pi*f_start*(t - tau) + ...
                    pi*S*(t - tau).^2));

                beat = tx .* conj(rx);

                h1 = a1 * beat;

                f_h1 = S * tau;

                H1_B(k,m) = sum(h1 .* ...
                    exp(-1j * 2*pi*f_h1*t));

            end
        end

        H1_B_phase = H1_B ./ abs(H1_B);
        h1_vec = H1_B_phase(:);

        rho_map(i_el,i_az) = ...
            abs(h2_vec' * h1_vec) / ...
            (norm(h2_vec) * norm(h1_vec));

    end
end