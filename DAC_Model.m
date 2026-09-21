function [y, error] = DAC_Model(x, step)

    y = step * round(x / step);
    
    if nargout > 1
        error = x - y;
    end

end