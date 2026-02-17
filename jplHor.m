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
end1   = datetime(2042, 02, 01);     % end of depature window
gd1 = (start1 : days(1) : end1)';    % generates array with timesteps of 1 day
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
ecc = 0.8;                             % eccentricity
inc = 0;                               % deg
alt = 20000;                           % km

% Orbit Geometry & Initial State
r_sc_p = R_mars + alt;                         % Radius at Periapsis (km)
a      = r_sc_p / (1 - ecc);                   % Semi-Major Axis of S/C Orbit (km)
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

% Plotting to Verify
figure;
plot(r_sc_mars(:,1), r_sc_mars(:,2));
hold on;
viscircles([0,0], R_mars, 'Color', 'r'); % Draw Mars
axis equal; grid on;
title('Spacecraft Trajectory (MCI Frame)');
xlabel('X (km)'); ylabel('Y (km)');

theta_p = pi/2 - atan(norm(r_sc_mars)/(R_sun + R_mars)); %radians

subplot(2,2,3)
plot(t, theta_p)
xaxis("time [days]")
yaxis("theta_p [rad]")
grid on

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
            'STEP_SIZE', '1d', ...
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


