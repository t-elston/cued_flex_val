% Generate synthetic data for two signals
fs = 1000; % Sampling frequency (Hz)
t = 0:1/fs:5; % Time vector (5 seconds)
f1 = 4; % Frequency of signal 1 (Hz)
f2 = 4; % Frequency of signal 2 (Hz)
x1 = sin(2*pi*f1*t); % Signal 1 (sine wave)
x2 = circshift(x1, fs/2);

% Add noise to the signals
noise_level = .1; % Noise level
x1 = x1 + noise_level*randn(size(t));
x2 = x2 + .5*randn(size(t));

% Compute the GPDC between the two signals
order = 10; % Model order for AR estimation
freq = 1:50; % Frequency range (up to Nyquist frequency)

X = [x1', x2'];

pdc = PDC_2(X, order, freq);

[GC_xy, GC_yx, order] = GrangerCausality2(x1', x2', order)

figure;
hold on
plot(pdc(:,1,2))
plot(pdc(:,2,1))