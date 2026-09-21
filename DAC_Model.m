function [y, error] = DAC_Model(x, step)

    if step == 0
        y = x;
    else
        y = step * round(x / step);
    end

    if nargout > 1
        error = x - y;
    end

end