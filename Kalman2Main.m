%% Mars Atmospheric Reentry Trajectory
clear; clc;
% format longg;
format short;
load Nominal.mat;
load SensorData.mat;
rng("default");
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
den0 = 1.0449e-07;
state0 = [V0, X0, y0, h0, den0];
dt =  0.001;
freq =50;
T = 1/freq;
tspan = 0:dt:T;

params = [Rm, mu, rho0, H, mass, S, Cd, LD, sig, dt];
xhatukf = state0';
xhatukf_Array = [state0];


Q = diag([0.005, 0.0001, 0.0001, 0.0001, 10^-7]);

% Q = diag([0.00, 0.000, 0.000, 0.000, 10^-7]);
%% Process Noise Covariance
P = eye(5); 
P(1, 1) = 1; P(2, 2) = 0.1; P(3, 3) = 0.1; P(4, 4) = 1; P(5, 5) = 10^-10;
Pukf = P;
W = ones(10, 1)/10; % UKF Weights
p11_array = [P(1, 1)];
p22_array = [P(2, 2)];
p33_array = [P(3, 3)];
p44_array = [P(4, 4)];
p55_array = [P(5, 5)];


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
    [root, ~] = chol(5*Pukf);
    xbreve(:, 1:5) = xhatukf + root(1:5,:)';
    xbreve(:, 6:10) = xhatukf - root(1:5,:)';
    
    % UKF time update
    for i = 1 : 10
        ODE_UKF = @(t, y) UKF(t, y, params);
        [t_UKF, y_UKF] = ode45(ODE_UKF, tspan, xbreve(:, i)', options);
        xbreve(:, i) = y_UKF(end, 1:5)';
    end

    xhatukf = zeros(5,1);
    
    for i = 1 : 10
        xhatukf = xhatukf + W(i) * xbreve(:,i);
    end

    Pukf = zeros(5,5);
    for i = 1 : 10
        Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
    end
    
    Pukf = Pukf + Q; 
    
%     if sum(sum(isinf(Pukf))) + sum(sum(isnan(Pukf))) > 0
%         Pukf = P0; % reinitialize the UKF covariance
%     end

    % UKF measurement update
    zukf = zeros(4, 10);
    zhat = 0;
    for i = 1 : 10
        Vval = xbreve(1, i);
        Xval = xbreve(2, i);
        yval = xbreve(3, i);
        hval = xbreve(4, i);

        zukf(1:4, i) = [Vval*cos(yval)*cos(Xval); Vval*cos(yval)*sin(Xval); -Vval*sin(yval); hval];
        zhat = zhat + W(i) * zukf(:,i); % Eq. 14.62 
    end

     
    Py = 0;
    Pxy = zeros(5,1);
    for i = 1 : 10
        Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';  %  Eq. 14.63
        Pxy = Pxy + W(i) * (xbreve(:,i) - xhatukf) * (zukf(:,i) - zhat)'; % Eq. 14.64
    end
    Py = Py + R;
    Kukf = Pxy / Py;
    xhatukf = xhatukf + Kukf * (z - zhat);
    xhatukf_Array = [xhatukf_Array; xhatukf'];
    p11_array = [p11_array; Pukf(1, 1)];
    p22_array = [p22_array; Pukf(2, 2)];
    p33_array = [p33_array; Pukf(3, 3)];
    p44_array = [p44_array; Pukf(4, 4)];
    p55_array = [p55_array; Pukf(5, 5)];
    Pukf = Pukf - Kukf * Py * Kukf';
    

    iter = iter + 1;
    fprintf("Running UKF time step %0.6f of %0.1f sec \n", t, tfinal);
end

xhatukf_Array = xhatukf_Array(xhatukf_Array(:, 4) >0, :)

Vukf = xhatukf_Array(:, 1);
Xukf = xhatukf_Array(:, 2);
yukf = xhatukf_Array(:, 3);
hukf = xhatukf_Array(:, 4);

Vxukf = Vukf.*cos(yukf).*cos(Xukf);
Vyukf = Vukf.*cos(yukf).*sin(Xukf);
Vzukf = -Vukf.*sin(yukf);

tint = 0:T:tfinal;
tint = tint(xhatukf_Array(:, 4) >0)
Posx = 0 + cumtrapz(tint, Vxukf);

Posy = 0 + cumtrapz(tint, Vyukf);

Posz = h0- cumtrapz(tint, Vzukf);

plot3(Posx/1000, Posy/1000, hukf/1000);
hold on;
plot3(nominal(:, 5) / 1000, nominal(:, 6)  / 1000, nominal(:, 7) / 1000);
axis square;
grid on;
[Posx(end)/1000, Posy(end)/1000, hukf(end)/1000];

save('UKF.mat', 'xhatukf_Array', 'tint', 'Posx', 'Posy', 'Posz', 'Vxukf',...
    'Vyukf', 'Vzukf', 'Vukf', 'Xukf', 'yukf', 'hukf',...
    'p11_array', 'p22_array', 'p33_array', 'p44_array', 'p55_array')
