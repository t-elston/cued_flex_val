function reconstructed_signal = reconstruct_from_phase_amplitude(phase, amplitude)
    
% reconstruct a waveform from instantaneous phase and amplitude estimates derived from the Hilbert transform

% Ensure that the input phase is in radians
    phase = mod(phase, 2*pi); % Make sure phase is in [0, 2*pi]

    % Reconstruct the complex signal from phase and amplitude
    complex_signal = amplitude .* exp(1i * phase);

    % Obtain the real part of the complex signal as the reconstructed waveform
    reconstructed_signal = real(complex_signal);
end