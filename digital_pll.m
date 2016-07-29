clearvars clc
% Digital PLL Implementation
%
% Experiment with the loop filter constants Kp and Ki to verify classical
% observations about PLL, such that frequency offsets can not be corrected
% if Ki=0;

nIterations = 2e3;
fs          = 1e6;   % Nominal sampling frequency
y_ppm       = 50;    % Fractional frequency offset in ppm
f0          = 1e3;   % Nominal clock frequency
pn_var      = 1e-9;  % Phase noise variance
Kp          = 0.05;  % Proportional Constant
Ki          = 0.01;  % Integral Constant

%% Derived parameters
y        = y_ppm * 1e-6;  % Fractional frequency offset in ppm
f_offset = y * f0;        % Absolute freq error (in Hz)

% Nominal Phase Increment (used in the loop phase accumulator)
delta_phi_c = 2*pi*f0/fs;

%% Generate input signal
% The input signal obviously does not depend on the loop processing, so we
% can generate it before simulating the loop.
%
% The input here is considered to be a pure complex sinusoid (exponential)
% slightly corrupted by phase noise and eventually with frequency offset.
% We can simulate such an exponential by considering that in each clock
% cycle (sample period), the exponential grows by the following increment:
delta_phi_in = 2*pi*(f0 + f_offset)/fs;
% which when combined to the following phase noise random sequence:
phase_noise  = sqrt(pn_var)*randn(nIterations, 1);
% results in the actual instantanous phase:
phi_in       = (0:nIterations-1).' * delta_phi_in + phase_noise;
% Finally, generate the exponential:
s_in         = exp(1j * phi_in);

%% Loop
% The loop has three main components:
%   1) Direct Digital Synthesizer (DDS): iteratively generates a a complex
%   sinusoid by computing e(j*phi_loop), where "phi_loop" is its
%   instantaneous phase that continuously grows based on a configurable
%   increment value. The angle increment added to the DDS (accumulated in
%   "phi_loop", i.e. its "phase accumulator") at each clock cycle (sample
%   period) is the nominal increment 2*pi*f0/fs plus a varying correction
%   term. The latter term should vary considerably until the point in which
%   the DDS output is synchronized to the input, namely until the loop
%   "locks".
%
%   2) Phase Detector: extracts the phase error from the input complex
%   exponential and the loop complex exponential
%
%   3) Loop Filter: filters the phase error such that, after sufficient
%   iterations, it is driven to zero. The output of the loop filter (called
%   filtered phase error) is the correction term used in the DDS.

% Preallocate
dds_loop           = zeros(nIterations, 1); % DDS Complex Value
phi_loop           = zeros(nIterations, 1); % DDS Phase Accumulator
dds_mult           = zeros(nIterations, 1); % Conjugate Product
phi_error          = zeros(nIterations, 1); % Phase Error
phi_error_filtered = zeros(nIterations, 1); % Filtered Phase Error

% Initialize the loop DDS to a random phase
phi_loop(1) = sqrt(pi)*randn;

for i = 1:nIterations

%% Loop DDS:
dds_loop(i) = exp(1j*phi_loop(i));

%% Phase Detector
% Phase error is obtained by computing the angle of the conjugate product
%   Input                 : e^(j*phi_in)
%   Loop                  : e^(j*phi_loop)
%   Conjugate Product     : e^[j*(phi_in - phi_loop)]
%   ---------------------------------------------------
%   Conj. Product Angle   : phi_in - phi_loop

% Multiply the input complex exponential by the conjugate of the loop DDS:
dds_mult(i) = s_in(i) * conj(dds_loop(i));

% Phase error:
phi_error(i) = angle(dds_mult(i));

%% Loop Filter
% The loop filter consists of a Proportional+Integral (PI) controller

% Proportional term
proportional_out = phi_error(i)*Kp;

% Integral term
integral_out     = phi_error(i)*Ki + integral_out;

% PI Filter output:
phi_error_filtered(i) = proportional_out + integral_out;

%% Update the phase accumulator of the loop DDS
% The angle of the DDS complex exponential in the next clock cycle results
% from the sum between the nominal phase increment and the filtered phase
% error term:
phi_loop(i+1) = phi_loop(i) + delta_phi_c + phi_error_filtered(i); 

end

%% Performance
% Input vs. Output Instantaneous Phase
figure
plot(phi_in)
hold on
plot(phi_loop, 'r')
title('Instantaneous Phase')
xlabel('Sample')
ylabel('Phase (rad)')
legend('Input', 'Output')

% Filtered Phase Error
figure
plot(phi_error_filtered)
title('Loop Filter Output')
xlabel('Sample')
ylabel('Error (rad)')

% Input vs. Output Sinusoids
figure
subplot(121)
plot(real(s_in))
hold on
plot(real(dds_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Real part')
legend('Input', 'Output')

subplot(122)
plot(imag(s_in))
hold on
plot(imag(dds_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Imaginary part')
legend('Input', 'Output')

