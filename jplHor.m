%----------------------
% Daniel Newallis
% Orbit Propogation of ARGO
% Gravational Parameters from JPL
%----------------------
% ----------------------
% Program setup
%-----------------------
clear
close all
clc
% ----------------------
% Time Window Setup
%-----------------------
% Depature Window (Gregorian dates)
start1 = datetime(2042, 01, 01);     % start of depature window
end1   = datetime(2042, 01, 03);     % end of depature window
gd1 = (start1 : minutes(1) : end1)';    % generates array with timesteps of 1 minute
jd1 = juliandate(gd1);               % Convert Gregian dates to Julian dates

% ----------------------
% Retrieve Planetary Ephemeris Data (SCI)
% ----------------------
[r_mars_vec, v_mars_vec] = get_jpl_horizons('499', start1, end1); % Mars (499)

% ----------------------
% Initial Information
%-----------------------
mu  = 1.32712440041279e11;             % Grav param of Sun (km^3/s^2)
mu_mars = 4.2828375214e4;              % Grav param of Mars (km^3/s^2)
R_mars = 3389.5;                       % Mars Radius (km)
R_sun = 696300;                        % Sun Radius (km)

% Spacecraft Orbit Parameters
ecc = 0.75;                             % eccentricity
inc = 0;                               % deg
alt = 25000;                           % km

% Orbit Geometry & Initial State
r_sc_p = R_mars + alt;                         % Radius at Apoapsis (km)
a      = r_sc_p / (1 + ecc);                   % Semi-Major Axis of S/C Orbit (km)
v_p    = sqrt(mu_mars * ((2/r_sc_p) - (1/a))); % Vis-Viva Velocity at Periapsis (km/s)
n      = sqrt(mu_mars / a^3);                  % Mean Motion (rad/s)
T_period = 2*pi*sqrt(a^3/mu_mars);             % Orbital Period (s)

% Initial State Vectors (MCI Frame)
r_p_MCI = [r_sc_p; 0; 0];               % Position Vector [x, y, z] (km)
v_p_MCI = [0; v_p; 0];                  % Velocity Vector [vx, vy, vz] (km/s)

% ----------------------
% Propagate Orbit
% ----------------------
% Define time vector for Spacecraft
t_total = seconds(end1 - start1);           % Total Duration in Seconds
dt_sc   = 60;                               % Time step: 60 Seconds
t_vec   = (0 : dt_sc : t_total)';           % Time Vector from 0 to end
r_sc_mars = zeros(length(t_vec), 3);        % Pre-allocate Output Matrix [Rx, Ry, Rz] in MCI Frame

for i = 1:length(t_vec)
    t_now = t_vec(i);
    M = n * t_now;    % M = M0 + n*t
    E = M;            % Initial Guess
    tol = 1e-12;     % Convergence tolerance
    max_iter = 50;
    
    for iter = 1:max_iter
        f_val = E - ecc*sin(E) - M;       % f(E)
        f_der = 1 - ecc*cos(E);           % f'(E)
        
        step = f_val / f_der;
        E_new = E - step;
        
        if abs(E_new - E) < tol
            E = E_new;
            break;
        end
        E = E_new;
    end
    
    % Calculate Position in Orbital Plane (MCI)
    x_orbit = a * (cos(E) - ecc);
    y_orbit = a * sqrt(1 - ecc^2) * sin(E);
    z_orbit = 0;
    
    % Store in matrix
    r_sc_mars(i, :) = [x_orbit, y_orbit, z_orbit];
end
% ----------------------
% Plotting to Verify (3D)
% ----------------------
figure('Color', 'w');
hold on;
grid on;
axis equal;

% Plot the Spacecraft Trajectory
plot3(r_sc_mars(:,1), r_sc_mars(:,2), r_sc_mars(:,3), 'b', 'LineWidth', 1.5);

% Plot Mars (As a 3D Sphere
[X, Y, Z] = sphere(50); 
X = X * R_mars;
Y = Y * R_mars;
Z = Z * R_mars;
mars_surf = surf(X, Y, Z);

% Label PLot
title('Spacecraft Trajectory (MCI Frame)', 'FontSize', 16);
xlabel('X (km)', 'FontSize', 14);
ylabel('Y (km)', 'FontSize', 14);
zlabel('Z (km)', 'FontSize', 14);

% Set viewing angle
view(45, 30);

% ----------------------
% Subplot
% ----------------------
r_mars_norms = vecnorm(r_mars_vec, 2, 2);
r_sc_norms   = vecnorm(r_sc_mars, 2, 2);

theta_e = acos( dot(r_mars_vec, r_sc_mars, 2) ./ (r_mars_norms .* r_sc_norms) );
theta_p = pi/2 - atan( r_sc_norms ./ (R_sun + R_mars) );

h_u = (R_mars .* r_mars_norms) ./ (R_sun - R_mars);
theta_u = atan(R_mars ./ h_u);
[S_percent, gd1] = Scalc(ecc,alt,R_mars,R_sun,mu_mars);

figure
subplot(2,2,1)
plot(gd1, theta_e)
ylabel("theta_e [rad]")
xlabel("time [days]")
ylim([0 pi]);
grid on
subplot(2,2,2)
plot(gd1, theta_u)
ylabel("theta_u [rad]")
xlabel("time [days]")
grid on
subplot(2,2,3)
plot(gd1, theta_p)
ylabel("theta_p [rad]")
xlabel("time [days]")
grid on
subplot(2,2,4)
plot(gd1, S_percent,'LineWidth',1.5);
ylabel("% Sunlight")
xlabel("time [days]")
ylim([0 100]);
grid on


function [S_percent, gd] = Scalc(e,alt,R_mars,R_sun,mu_mars)
    r_a = R_mars + alt;             % apoapsis radius (km)
    a  = r_a/(1+e);
    
    % Time Setup
    start = datetime(2042,01,01,0,0,0);
    END = datetime(2042,01,03,0,0,0);
    gd = (start:minutes(1):END)';
    N = length(gd);

    % Get Mars heliocentric position
    [r_mars_vec, ~] = get_jpl_horizons('499', start, END);

    % Safety check
    if length(r_mars_vec) ~= N
        error('Horizons time vector does not match gd1 length.')
    end

    % Preallocate
    S_percent = zeros(N,1);

    % Mean motion
    n = sqrt(mu_mars/a^3);

    for k = 1:N
    
        t = (k-1)*60;    % 60 seconds per step (1 minute resolution)
    
        % Mean anomaly
        M = n*t;
    
        % Solve Kepler's Equation (Newton iteration)
        E = M;
        for iter = 1:20
            E = E - (E - e*sin(E) - M)/(1 - e*cos(E));
        end
    
        % True anomaly
        nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    
        % Radius
        r = a*(1 - e*cos(E));
    
        % Spacecraft Mars-centered position
        r_sc_mars = [r*cos(nu), r*sin(nu), 0];
    
        % Spacecraft heliocentric position
        r_sc_sun = r_mars_vec(k,:) + r_sc_mars;
    
        % Distances
        d_sun  = norm(r_sc_sun);
        d_mars = norm(r_sc_mars);
    
        % Angular radii
        alpha_s = asin( min(1, max(-1, R_sun/d_sun)) );
        alpha_m = asin( min(1, max(-1, R_mars/d_mars)) );
    
        % Angular separation
        cos_theta = dot(-r_sc_sun, -r_sc_mars)/(d_sun*d_mars);
        cos_theta = min(1, max(-1, cos_theta));
        theta = acos(cos_theta);
    
        % Eclipse logic
        if theta >= alpha_s + alpha_m
            S_percent(k) = 100;
        
        elseif alpha_m >= alpha_s && theta <= (alpha_m - alpha_s)
            S_percent(k) = 0;
        
        else
            % Partial eclipse (safe overlap math)
            r1 = alpha_s;
            r2 = alpha_m;
            d  = theta;
        
            % Prevent divide-by-zero
            if d == 0
                S_percent(k) = 0;
            else
                arg1 = (d^2 + r1^2 - r2^2)/(2*d*r1);
                arg2 = (d^2 + r2^2 - r1^2)/(2*d*r2);
            
                arg1 = min(1,max(-1,arg1));
                arg2 = min(1,max(-1,arg2));
            
                part1 = r1^2 * acos(arg1);
                part2 = r2^2 * acos(arg2);
                part3 = 0.5*sqrt(max(0,(-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2)));
            
                overlap = part1 + part2 - part3;
                S_percent(k) = 100*(1 - overlap/(pi*r1^2));
            end
        end
    end
    
% Clean numerical noise
S_percent = real(S_percent);
S_percent(isnan(S_percent)) = 0;
S_percent = max(0,min(100,S_percent));
end

% ----------------------
% JPL HORIZONS API
% ----------------------
function [r_vec, v_vec] = get_jpl_horizons(body_id, t_start, t_end)
    % INPUTS
    % body_id: ('499' for Mars)
    % t_start: datetime object
    % t_end: datetime object
    base_url = 'https://ssd.jpl.nasa.gov/api/horizons.api';
    
    % Format dates for API (YYYY-MM-DD)
    s_str = datestr(t_start, 'yyyy-mm-dd');
    e_str = datestr(t_end, 'yyyy-mm-dd');
    
        % Fetch Data using explicit Name-Value pairs
        raw_text = webread(base_url, ...
            'format', 'text', ...
            'COMMAND', body_id, ...
            'OBJ_DATA', 'YES', ...
            'MAKE_EPHEM', 'YES', ...
            'EPHEM_TYPE', 'VECTORS', ...
            'CENTER', '@10', ...
            'START_TIME', s_str, ...
            'STOP_TIME', e_str, ...
            'STEP_SIZE', '1m', ...
            'CSV_FORMAT', 'YES');
            
        % Parse Text
        % Data is located between $$SOE and $$EOE markers
        soe_idx = strfind(raw_text, '$$SOE');
        eoe_idx = strfind(raw_text, '$$EOE');
        
        if isempty(soe_idx) || isempty(eoe_idx)
             error('No Ephemeris Data Found between markers.');
        end

        % Extract the CSV content segment
        data_str = raw_text(soe_idx+5 : eoe_idx-1);
        
        % Read the data using textscan
        % Format: JDTDB, Calendar, X, Y, Z, VX, VY, VZ
        C = textscan(data_str, '%f %s %f %f %f %f %f %f %[^\n]', ...
            'Delimiter', ',', 'MultipleDelimsAsOne', false);
        % Extract Position (km) and Velocity (km/s)
        % Columns: 3=X, 4=Y, 5=Z, 6=VX, 7=VY, 8=VZ
        r_vec = [C{3}, C{4}, C{5}];
        v_vec = [C{6}, C{7}, C{8}];
end
