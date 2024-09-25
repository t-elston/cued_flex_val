function lagged_corr = calculate_lagged_cross_correlation(signal1, signal2, max_lag)
    % Inputs:
    % signal1, signal2: Time series data for brain areas (column vectors)
    % max_lag: Maximum lag to consider
    
    % Calculate the cross-correlation at different lags
    [corr_values, lags] = xcorr(signal1, signal2, max_lag, 'coeff');
    
    % Find the lag with the maximum correlation
    [~, max_corr_idx] = max(corr_values);
    optimal_lag = lags(max_corr_idx);
    
%     % Plot the cross-correlation function
%     figure;
%     stem(lags, corr_values);
%     xlabel('Lag');
%     ylabel('Cross-Correlation');
%     title('Lagged Cross-Correlation');
%     grid on;
%     
    % Return the lag with the maximum correlation
    lagged_corr = optimal_lag;
end