% parameters
c = 299792458;

f_start = 58.1e9;
f_stop = 63.1e9;

B = f_stop - f_start;

t_chirp = 100e-6;
S = B/t_chirp;
f_sampling = 2e6;
N = 128;

a1 = 1;
a2 = 0.3;

f_c = (f_start + f_stop) / 2;
lambda = c / f_c;

d = lambda / 2;

% antenna array
rx_pos = [
    0  0;
    d  0;
    0  d
];

num_rx = size(rx_pos, 1);

% coordinates
az = deg2rad(45);
el = deg2rad(0);
range_A = 0.5;