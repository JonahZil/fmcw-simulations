% spatial phases
H1_relative = H1 ./ H1(:,1);
H2_relative = H2 ./ H2(:,1);

fprintf('H1 RX phase relative to RX1 at 100 us:\n');
disp(rad2deg(angle(H1_relative(1,:))));

fprintf('H2 RX phase relative to RX1 at 100 us:\n');
disp(rad2deg(angle(H2_relative(1,:))));

% correlation

H1_phase = H1 ./ abs(H1);
H2_phase = H2 ./ abs(H2);

h1_vec = H1_phase(:);
h2_vec = H2_phase(:);

rho = abs(h1_vec' * h2_vec) / (norm(h1_vec) * norm(h2_vec));

fprintf('H1/H2 joint chirp-RX correlation: %.10f\n', rho);