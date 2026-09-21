function n = Thermal_Noise(N, adc_fs, bandwidth, rms_noise)

    [b, a] = butter(2, bandwidth/(adc_fs/2), 'low');
    
    n = filter(b, a, randn(1, N));
    n = n - mean(n);
    n = rms_noise * n / std(n);

end