function mean_angle = circular_mean(angles)
    % Convert angles to complex numbers on the unit circle
    z = exp(1i * angles);
    
    % Calculate the mean angle using the complex number representation
    mean_z = mean(z);
    
    % Convert the mean complex number back to an angle
    mean_angle = angle(mean_z);
    
    % Ensure the result is between 0 and 2*pi (0 to 360 degrees)
    if mean_angle < 0
        mean_angle = 2*pi + mean_angle;
    end
end