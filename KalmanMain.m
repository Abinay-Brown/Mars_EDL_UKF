%% Mars Atmospheric Reentry Trajectory
clear; clc;
% format longg;
format short;
load Nominal.mat;
load SensorData.mat;

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
V0 = 5.4 * 1000;
X0 = 0   * (pi/180);
y0 = -11.5 * (pi/180);
h0 = 135 * 1000;
state0 = [V0, X0, y0, h0];
dt =  0.001;
freq =50;
T = 1/freq;
tspan = 0:dt:T;

params = [Rm, mu, rho0, H, mass, S, Cd, LD, sig, dt];
xhatukf = state0';
xhatukf_Array = [state0];


Q = diag([0, 0, 0, 0]);
%% Process Noise Covariance
P = eye(4); 
P(1, 1) = 1; P(2, 2) = 0.1; P(3, 3) = 0.1; P(4, 4) = 1;
Pukf = P;
W = ones(8, 1)/8; % UKF Weights



%% Simulation Run
options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
iter = 1;
tfinal = t(end);

for t = T : T : tfinal
    
    R = [sigVx_sq(iter), 0, 0, 0; 
         0, sigVy_sq(iter), 0, 0
         0, 0, sigVz_sq(iter), 0
         0, 0, 0, sigH_sq]; 

    z = SensorObs(:, iter);

    % Generate the UKF sigma points
    [root, ~] = chol(4*Pukf);
    xbreve(:, 1:4) = xhatukf + root(1:4,:)';
    xbreve(:, 5:8) = xhatukf - root(1:4,:)';
    
    % UKF time update
    for i = 1 : 8
        ODE_UKF = @(t, y) UKF(t, y, params);
        [t_UKF, y_UKF] = ode45(ODE_UKF, tspan, xbreve(:, i)', options);
        xbreve(:, i) = y_UKF(end, 1:4)';
    end

    xhatukf = zeros(4,1);
    
    for i = 1 : 8
        xhatukf = xhatukf + W(i) * xbreve(:,i);
    end

    Pukf = zeros(4,4);
    for i = 1 : 8
        Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
    end
    Pukf = Pukf + Q; 
    
%     if sum(sum(isinf(Pukf))) + sum(sum(isnan(Pukf))) > 0
%         Pukf = P0; % reinitialize the UKF covariance
%     end

    % UKF measurement update
    zukf = zeros(4, 8);
    zhat = 0;
    for i = 1 : 8
        Vval = xbreve(1, i);
        Xval = xbreve(2, i);
        yval = xbreve(3, i);
        hval = xbreve(4, i);

        zukf(1:4, i) = [Vval*cos(yval)*cos(Xval); Vval*cos(yval)*sin(Xval); -Vval*sin(yval); hval];
        zhat = zhat + W(i) * zukf(:,i); % Eq. 14.62 
    end

     
    Py = 0;
    Pxy = zeros(4,1);
    for i = 1 : 8
        Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';  %  Eq. 14.63
        Pxy = Pxy + W(i) * (xbreve(:,i) - xhatukf) * (zukf(:,i) - zhat)'; % Eq. 14.64
    end
    Py = Py + R;
    Kukf = Pxy / Py;
    xhatukf = xhatukf + Kukf * (z - zhat);
    xhatukf_Array = [xhatukf_Array; xhatukf'];
    
    Pukf = Pukf - Kukf * Py * Kukf';
    

    iter = iter + 1;
    disp(iter);
end

Vukf = xhatukf_Array(:, 1);
Xukf = xhatukf_Array(:, 2);
yukf = xhatukf_Array(:, 3);
hukf = xhatukf_Array(:, 4);

Vxukf = Vukf.*cos(yukf).*cos(Xukf);
Vyukf = Vukf.*cos(yukf).*sin(Xukf);
Vzukf = -Vukf.*sin(yukf);

tint = 0:T:tfinal;

Posx = 0 + cumtrapz(tint, Vxukf);

Posy = 0 + cumtrapz(tint, Vyukf);

Posz = h0- cumtrapz(tint, Vzukf);

plot3(Posx/1000, Posy/1000, hukf/1000);
hold on;
plot3(nominal(:, 5) / 1000, nominal(:, 6)  / 1000, nominal(:, 7) / 1000);
axis square;
grid on;
[Posx(end)/1000, Posy(end)/1000, hukf(end)/1000];