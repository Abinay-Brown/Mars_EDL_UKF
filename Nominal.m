%% Mars Atmospheric Reentry Trajectory
clear; clc;
format longg;
%% Mars Constants
Rm = 3376.2 * 1000;     % Mars Radius (m)
mu = 4.2828 * 10^13;    % Mars Gravitational Parameter (m^3/s^2)
rho0 = 0.02;            % Mars Sea-Level Density kg/m^3
H = 11.1;               % Mars scale height (km)
    
%% Reentry Vehicle Parameters
mass   = 3520;           % Mass in kg
radius = 4.5;          % Shell Radius
S      = pi*radius^2;   % Shell Area   
Cd     = 1.7;           % Shell Drag Coefficient
LD     = 0.3;          % Vehicle Light Drag Ratio
sig    = 2.5 * (pi/180);            % Bank Angle

%% Initial Conditions
V0 = 5.5 * 1000;
X0 = 0   * (pi/180);
y0 = -12 * (pi/180);
h0 = 135 * 1000;
state0 = [V0, X0, y0, h0, 0, 0, h0];
dt = 0.001;
tspan = 0:dt:1310;
params = [Rm, mu, rho0, H, mass, S, Cd, LD, sig, dt];


%% Simulation Run
options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
func = @(t, y) ReentryDynamics(t, y, params);
[t, nominal] = ode45(func, tspan, state0, options);
nominal = nominal(nominal(:,4)>0, :);
t = t(nominal(:,4)>0);
plot3(nominal(:, 5) / 1000, nominal(:, 6)  / 1000, nominal(:, 7) / 1000);
axis square;
grid on;
save('Nominal.mat', "state0", "params", "t" ,"nominal");
