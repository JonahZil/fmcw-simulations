function y = If_Amp_LowPass_Filter(dt, x, gain)

    order = 2;
    cutoff = 800e3;
    
    fs = 1/dt;
    [b, a] = butter(order, cutoff/(fs/2), 'low');
    
    y = gain * filter(b, a, x);

end