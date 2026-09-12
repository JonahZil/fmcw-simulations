% ambiguity map

figure;

imagesc(azimuths, elevations, rho_map);

axis xy;
axis square;

colorbar;
caxis([0 1]);

xlabel('Target B Azimuth (deg)');
ylabel('Target B Elevation (deg)');

title(sprintf( ...
    'H2_A / H1_B Ambiguity: A = (%.0f deg, %.0f deg)', ...
    rad2deg(az_A), rad2deg(el_A)));

hold on;

contour(azimuths, elevations, rho_map, ...
    [0.9 0.95 0.99], 'k', 'ShowText', 'on');

hold off;