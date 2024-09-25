function [pdc] = PDC_2(X, pmax, freq)
%PDC Computes the Partial Directed Coherence (PDC) between signals using an AR model.
%
%   [pdc] = PDC(X, pmax, freq) computes the PDC between signals in matrix X
%   using an autoregressive (AR) model with maximum model order pmax. The PDC
%   is computed at frequencies specified by the vector freq.
%
%   Inputs:
%       X: Input data matrix where each column represents a different time series or channel.
%       pmax: Maximum model order for the AR model.
%       freq: Frequency vector for which PDC will be computed.
%
%   Output:
%       pdc: Partial Directed Coherence matrix containing PDC values in form 
%       frequency x channel i x channel j, where causality is from channel i to channel j

% The ARfit package for computing AR models is available here: https://climate-dynamics.org/software/#arfit
% it is required in order for this function to run

%----------------------------------------------------------------------------------------------------------
% Get the number of channels
Nc = size(X, 2);

% Initialize the PDC matrix
pdc = zeros(length(freq), Nc, Nc);

% Fit an AR model to the input data
[w, B, Sig] = arfit(X, 1, pmax, 'fpe', 'zero');

% Handle NaN (Not a Number) case if model fitting fails
if isnan(Sig)
    pdc = NaN;
    return
end

% Convert AR coefficients to the frequency domain
Bf = getBf(Nc, B, freq);

% Compute PDC for each channel pair at each frequency
for xx = 1:Nc
    for yy = 1:Nc
        if (yy ~= xx)
            for ff = 1:length(freq)
                % Compute PDC using the spectral matrices
                pdc(ff, yy, xx) = abs(Bf(xx, yy, ff)) / sqrt(Bf(:, yy, ff)' * Bf(:, yy, ff));
            end
        end
    end
end

end

function [Bf] = getBf(m, B, freq)
%   Computes the frequency response of the AR coefficients.
%
%   [Bf] = getBf(m, B, freq) computes the frequency response of the AR coefficients
%   B for each channel at frequencies specified by the vector freq.
%
%   Inputs:
%       m: Number of channels.
%       B: AR coefficients.
%       freq: Frequency vector.
%
%   Output:
%       Bf: Frequency response of the AR coefficients.

% Compute the number of parameters per channel
popt = size(B, 2) / m;

% Initialize Bf matrix
Bf = zeros(m, m, length(freq));

% Initialize identity matrix
A0 = eye(m);

% Compute frequency response for each frequency
for g = 1:length(freq)
    for ch = 1:m
        Bj = B(:, [ch:m:m * popt]);
        for j = 1:popt
            Bf(:, ch, g) = Bf(:, ch, g) + Bj(:, j) .* exp(1i * j * freq(g));
        end
    end
    % Compute the inverse frequency response matrix
    Bf(:, :, g) = A0 - Bf(:, :, g);
end

end