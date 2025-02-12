% Function for initial loading of particles
function [position, velocity] = InitialLoading(N_particles, V0a, V0b, Vth, Amplitude, Mode, L)
    position = transpose(linspace(0, L - L/N_particles, N_particles));

    %velocity = Vth*randn(N_particles, 1) + V0; % Maxwellian distribution of velocities


    % for electron twostream
    velocity_1 = Vth*randn(N_particles/2, 1) + V0a;
    velocity_2 = Vth*randn(N_particles/2, 1) + V0b;



    velocity = [velocity_1;velocity_2];
    % velocity = velocity(randperm(length(velocity)));
    %velocity = -1 + 2 * rand(N_particles, 1); 

    % V_min = -1.0; % Minimum velocity
    % V_max = 1.0; % Maximum velocity
    % % Generate the uniform velocity distribution
    % velocity = V_min + (V_max - V_min) * rand(N_particles, 1);

    if Amplitude ~= 0
        position = position + Amplitude*cos(2*pi*Mode*position/L); % Perturbation
    end
end

