function [GC_xy, GC_yx] = GrangerCausality(X, Y, max_order)
%   Computes Granger causality between two time series.
%
%   [GC_xy, GC_yx] = computeGrangerCausality(X, Y, max_order) computes the
%   Granger causality between time series X and Y up to the specified maximum
%   model order using a multivariate autoregressive (MVAR) model. The function
%   also returns the automatically determined model order.
%
%   Inputs:
%       X: Time series data for the first variable (channel).
%       Y: Time series data for the second variable (channel).
%       max_order: Maximum model order to consider.
%
%   Outputs:
%       GC_xy: Granger causality from X to Y.
%       GC_yx: Granger causality from Y to X.

% Initialize variables
num_samples = min(length(X), length(Y));
aic_values = zeros(max_order, 1);

warning('off', 'MATLAB:nearlySingularMatrix');

% Concatenate the two time series data
data = [X(1:num_samples), Y(1:num_samples)];

% Determine the optimal model order using AIC
for order = 1:max_order
    % Fit VAR model with current order
    [coefficients, ~] = estimateVAR(data, order);
    
    % Compute residuals
    residuals = computeResiduals(data, coefficients, order);
    
    % Compute residual covariance matrix
    residual_covariance = cov(residuals);
    
    % Compute AIC
    aic_values(order) = log(det(residual_covariance)) + 2 * order * size(data, 2) / num_samples;
end

% Find the model order with the minimum AIC
[~, order] = min(aic_values);

% Fit the MVAR model with the determined order
[coefficients, ~] = estimateVAR(data, order);
residuals = computeResiduals(data, coefficients, order);
residual_covariance = cov(residuals);


% Compute Granger causality from X to Y
residuals_lagged_X = lagmatrix(data(:, 1), 1:order); % Lagged values of X
num_lags_X = min(size(residuals, 1), size(residuals_lagged_X, 1));
residuals_lagged_X = residuals_lagged_X(order+1:end, :);
residuals_lagged_X = residuals_lagged_X(1:num_lags_X, :);  % Ensure consistent number of rows
residual_covariance_X = cov([residuals, residuals_lagged_X]);

GC_xy = real(log(det(residual_covariance) / det(residual_covariance_X)));

% Compute Granger causality from Y to X
residuals_lagged_Y = lagmatrix(data(:, 2), 1:order);
num_lags_Y = min(size(residuals, 1), size(residuals_lagged_Y, 1));
residuals_lagged_Y = residuals_lagged_Y(order+1:end, :);
residuals_lagged_Y = residuals_lagged_Y(1:num_lags_Y, :);  % Ensure consistent number of rows
residual_covariance_Y = cov([residuals, residuals_lagged_Y]);

GC_yx = real(log(det(residual_covariance) / det(residual_covariance_Y)));
end

function [coefficients, constant] = estimateVAR(data, order)
% ESTIMATEVAR Estimates VAR model coefficients using least squares.
%
%   [coefficients, constant] = estimateVAR(data, order) fits a VAR model
%   of specified order to the input data using least squares estimation.
%
%   Inputs:
%       data: Input data matrix where each column represents a variable.
%       order: Model order for the VAR model.
%
%   Outputs:
%       coefficients: Coefficients matrix of the VAR model.
%       constant: Constant vector of the VAR model.

num_vars = size(data, 2);
num_samples = size(data, 1) - order;

% Construct lagged data matrix
X = zeros(num_samples, order * num_vars);
for lag = 1:order
    X(:, (lag - 1) * num_vars + 1 : lag * num_vars) = data(order - lag + 1 : end - lag, :);
end

% Target data (response)
Y = data(order + 1 : end, :);

% Estimate coefficients using least squares
coefficients = (X' * X) \ (X' * Y);

% Extract constant term
constant = mean(Y - X * coefficients, 1);

end

function residuals = computeResiduals(data, coefficients, order)
% COMPUTERESIDUALS Computes residuals of a VAR model.
%
%   residuals = computeResiduals(data, coefficients, order) computes the
%   residuals of a VAR model given the input data and model coefficients.
%
%   Inputs:
%       data: Input data matrix where each column represents a variable.
%       coefficients: Coefficients matrix of the VAR model.
%       order: Model order for the VAR model.
%
%   Output:
%       residuals: Residuals matrix of the VAR model.

num_vars = size(data, 2);
num_samples = size(data, 1) - order;

% Construct lagged data matrix
X = zeros(num_samples, order * num_vars);
for lag = 1:order
    X(:, (lag - 1) * num_vars + 1 : lag * num_vars) = data(order - lag + 1 : end - lag, :);
end

% Target data (response)
Y = data(order + 1 : end, :);

% Compute residuals
residuals = Y - X * coefficients;

end