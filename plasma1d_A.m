clear

% Parameters
L = 128;
dt = 0.05; 
Nt = 2000; 
Ng = 512; 

Num_beams = 2;     % have to remove beam

% Beam 1
N1 = 50000; 
V01a = 10;      %----
V01b = 10;      %----
Vth1 = 1; 
QM1 = -1; 
m1   = 1;
XP1 = 0; 
Mode1 = 0; 
WP = 1;
% Beam 2
N2 = 50000; 
V02a =  -10;        %----
V02b =  -10;        %----
Vth2 = 1; 
QM2 = -1;
m2  = 1;
XP2 = 0; 
Mode2 = 0;
% Ions in the background
IB = 100000;%0*20000;

% Size of the cell
dx = L/(Ng);
time = 0:dt:(Nt-1)*dt;

% Methods
Motion_method = 'Leapfrog';  % 'Leapfrog' 'Runge Kutta (RK4)' 'Euler method'
Field_method = 'Fast Fourier Transform (FFT)';  % 'Finite Difference Method' 'Fast Fourier Transform (FFT)' 'Direct Integration'
Interpolation_method = 'Cloud in Cell (CIC)'; % 'Nearest Grid Point (NGP)'

% Charge
[Q1, Q2, rho_back] = Charge(QM1, QM2, IB, N1, N2, L, WP);     % no idea why
%rho_back = IB/Ng;
% Initial loading
[xp1, vp1] = InitialLoading(N1, V01a, V01b, Vth1, XP1, Mode1, L);
[xp2, vp2] = InitialLoading(N2, V02a, V02b, Vth2, XP2, Mode2, L);

% figure;
% scatter(xp1, vp1,'red','.');
% hold on;
% scatter(xp2, vp2,'blue', '.');
% xlabel('Position');
% ylabel('Velocity');
% 
% 
% figure;
% histogram(vp1,50,"Normalization","count");
% hold on;
% histogram(vp2,50,"Normalization","count");
% xlabel('Velocity');
% ylabel('Number of Particles');

% figure;
% histogram([vp1,vp2], 100, 'Normalization', 'count', 'FaceColor', 'b', 'EdgeColor', 'b');
% xlabel('$v$ (Velocity)', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('Number of particles', 'Interpreter', 'latex', 'FontSize', 12);
% pause(0.00001);
% filename = sprintf('histogram_plot.png');
% % Save the figure as a PNG file
% saveas(gcf, filename);
% % Close the figure to save memory
% close(gcf);

%save_2

% Auxiliarity vectors
p1 = AuxVector(N1); p2 = AuxVector(N2);

% Poisson's equation preparation
un = ones(Ng-1, 1); % Ng-1 * 1
Ax = spdiags([un -2*un un], [-1 0 1], Ng-1, Ng-1); % Matrix that gives the indices of the Poisson equation Ng-1 * Ng-1.
kap = (2*pi/L)*(-Ng/2:Ng/2-1);
kap = fftshift(kap');
kap(1) = 1;

% Others
mat1 = 0; 
mat2 = 0; 
Eg = 0; 
Phi = 0;

disp = [];



%figure;
% Computational cycle
for it = 1:Nt

    mom(it) = (Q1/QM1)*(sum(vp1)) + (Q2/QM2)*(sum(vp2));    %Q1/  Q2/ 
    E_kin(it) = 0.5*abs(Q1/QM1)*(sum(vp1.^2)) + 0.5*abs(Q2/QM2)*(sum(vp2.^2)); %Q1/ Q2/
    E_pot(it) = 0.5*sum(Eg.^2)*dx;
    E_tot(it) = sum(Eg);
    FLE(it)   = Eg(1,1);
    FLP(it)   = Phi(1,1);
    phi_total(it) = sum(Phi); 



    switch Motion_method
        case 'Leapfrog'
            vp1 = MotionV(vp1, QM1, mat1, Eg, N1, Motion_method, dt,m1);
            vp2 = MotionV(vp2, QM2, mat2, Eg, N2, Motion_method, dt,m2);
    end

    % Updating positions
    xp1 = MotionX(xp1, vp1, Motion_method, dt);
    xp2 = MotionX(xp2, vp2, Motion_method, dt);

    % Periodic boundary conditions
    xp1 = PeriodicBC(xp1, 0, L);
    xp2 = PeriodicBC(xp2, 0, L);

    % Interpolation functions
    mat1 = InterpolationF(Interpolation_method, dx, Ng, xp1, N1, p1);
    mat2 = InterpolationF(Interpolation_method, dx, Ng, xp2, N2, p2);

    % Charge density
    rho1 = Charge_density(Q1, mat1, dx);
    rho2 = Charge_density(Q2, mat2, dx);
    rhot = rho1 + rho2 + rho_back;

    % Field equations
    [Phi, Eg] = Field(Field_method, rhot, Ng, dx, L, Ax, kap);

    % Updating velocity
    vp1 = MotionV(vp1, QM1, mat1, Eg, N1, Motion_method, dt,m1);
    vp2 = MotionV(vp2, QM2, mat2, Eg, N2, Motion_method, dt,m2);

    % Append the new array to the resultArray
    % disp = [disp, Eg];

    iteration = it

    %save_1;

    % if rem(it,100)==0
    scatter(xp1, vp1, 5, 'red', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); % Smaller, semi-transparent red points
    hold on
    scatter(xp2, vp2, 5, 'blue', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); % Smaller, semi-transparent blue points
    xlabel('$x$ (Position)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$v$ (Velocity)', 'Interpreter', 'latex', 'FontSize', 12);
    pause(0.00001)
    hold off
    %end

    % 
    
    % histogram(vp1,50,"Normalization","count");
    % hold on
    % histogram(vp2,50,"Normalization","count");
    % pause(0.00001);
    % hold off

    % if rem(it,2500)==0
    % histogram([vp1,vp2], 100, 'Normalization', 'count', 'FaceColor', 'b', 'EdgeColor', 'b'); % Red color for vp1
    % % hold on
    % % histogram(vp2, 50, 'Normalization', 'count', 'FaceColor', 'b', 'EdgeColor', 'b'); % Blue color for vp2
    % % Add axis labels with LaTeX-style formatting
    % xlabel('$v$ (Velocity)', 'Interpreter', 'latex', 'FontSize', 12);
    % ylabel('Number of particles', 'Interpreter', 'latex', 'FontSize', 12);
    % % Optional: Add title
    % % title('Particle Velocity Distribution', 'Interpreter', 'latex', 'FontSize', 14);
    % pause(0.001);
    % % hold off
    % filename = sprintf('histogram_plot_%d.png', it);
    % % Save the figure as a PNG file
    % saveas(gcf, filename);
    % % Close the figure to save memory
    % close(gcf);
    % end 

    % Plotting and saving every 100 iterations
    % if rem(it,10)==0
    %     figure; % Create a new figure
    %     scatter(xp1, vp1, 'red', '.');
    %     hold on;
    %     scatter(xp2, vp2, 'blue', '.');
    %     hold off;
    % 
    %     % Construct the filename
    %     filename = sprintf('scatter_plot_%d.png', it);
    % 
    %     % Save the figure as a PNG file
    %     saveas(gcf, filename);
    % 
    %     % Close the figure to save memory
    %     close(gcf);
    % end 
    

end

% DR = (fftshift(fft2(disp))'); %fftshift
% 
% % Frequency and wavenumber axes
% K = ((2*pi/L)*(-Ng/2:Ng/2-1));   %- floor(Ng/2)
% omega = (((0:Nt-1)- floor(Nt/2)) * (2 * pi / (Nt * dt)));  %  - floor(Nt/2)





% figure;
% imagesc(disp);
% colorbar;
% xlabel('Time');
% ylabel('Position');

% figure;
% imagesc(K, flip(omega), abs(DR));   %image
% colorbar;
% xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 12);
% %title('1D Electron Dispersion Relation')
% grid on;
% % Construct the filename
% filename = sprintf('dispersion_plot_%d.png', it);
% 
% % Save the figure as a PNG file
% saveas(gcf, filename);

% figure;
% plot(time, E_tot/sqrt(max(E_tot)^2));
% xlabel('$t$ (Timestep)', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('Total Electrostatic field', 'Interpreter', 'latex', 'FontSize', 12);
% %title('1D Electron Dispersion Relation')
% %grid on;
% % Construct the filename
% filename = sprintf('E_field_plot.png');
% 
% % Save the figure as a PNG file
% saveas(gcf, filename);
% xlabel('Time')
% ylabel('Total Electric Field')
% 
% figure;
% plot(time, mom)
% xlabel('Time')
% ylabel('Momentum')
% 
% figure;
% plot(time, E_kin, 'b')
% xlabel('Time')
% ylabel('Kinetic Energy')
% 
% figure;
% plot(time, E_pot, 'r')
% xlabel('Time')
% ylabel('Potential Energy')
% 
figure;
hold on
plot(time, E_kin, 'b');
plot(time, E_pot, 'r');
plot(time, E_kin + E_pot, 'k');
xlabel('$t$ (Timestep)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 12);
legend('show'); % Display the legend
legend('Location', 'best'); % Position the legend in the best location

%title('1D Electron Dispersion Relation')
%grid on;
% Construct the filename
filename = sprintf('energy_plot.png');

% Save the figure as a PNG file
saveas(gcf, filename);
% xlabel('Time')
% ylabel('Energy')

% figure;
% plot(time, FLE)
% xlabel('Time')
% ylabel('Electric Field at point (1,1)')
% 
% figure;
% plot(time, FLP)
% xlabel('Time')
% ylabel('Electrostatic Potential at point (1,1)')
% 
% figure;
% plot(time, phi_total)
% xlabel('Time')
% ylabel('Total Potential')
% 
% figure;
% scatter(xp1, vp1,'red','.');
% hold on;
% scatter(xp2, vp2,'blue', '.');
% xlabel('Position');
% ylabel('Velocity');
% 
% figure;
% histogram(vp1,50,"Normalization","count");
% hold on;
% histogram(vp2,50,"Normalization","count");
% xlabel('Velocity');
% ylabel('Number of Particles');