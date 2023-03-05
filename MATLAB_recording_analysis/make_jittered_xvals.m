function [out_x] = make_jittered_xvals(indata, x_range)

% x_range: [x_min, x_max]


out_x = (x_range(2)-x_range(1)).*rand(numel(indata),1) + x_range(1);

end % of function