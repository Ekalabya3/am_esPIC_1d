% FFT ANALYSIS
% This script performs FFT analysis on the electric field data at a particular point.

% Load the data from an Excel file
data = readtable('C:\Users\abhil\Desktop\Images\1d EQ\total potential eq 1.xlsx');  % Update the file path if necessary
timestep = data.TIME;        % Assuming the column name for time is 'TIME'
electricField = data.E_FIELD;  % Assuming the column name for electric field is 'E_FIELD'

% Calculate delta_t from the timestep data
delta_t = mean(diff(timestep));  % Computational timestep in seconds

% Sampling frequency and interval
Fs = 1 / delta_t;  % Sampling frequency in Hz
T = delta_t;       % Sampling period in seconds

% Length of the signal
L = length(electricField);

% Perform FFT
Y = fft(electricField);

% Compute the two-sided spectrum P2 and the single-sided spectrum P1
P2 = abs(Y / L);               % Two-sided spectrum
P1 = P2(1:L/2+1);              % Single-sided spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Correct the amplitude

% Define the frequency domain f
f = Fs * (0:(L/2)) / L;

% Plot the single-sided amplitude spectrum
figure;
plot(f, P1);
title('Single-sided Amplitude Spectrum of Potential');
xlabel('Frequency (Arbitrary Unit)');
ylabel('|P1(f)|');

% Save the figure (optional)
% Uncomment the following line if you want to save the plot
% saveas(gcf, 'fft_spectrum.png');

