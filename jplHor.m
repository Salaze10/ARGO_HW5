%----------------------
% Daniel Newallis (Modified)
% PorkChop plot for Earth to Mars Transfer
% Lambert's Problem Solver Universal Variable Method
% Assume Prograde Solution
% Planetary Data from JPL Horizons API (Web)
% One Day Time Steps During Time Windows
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
start1 = datetime(2038, 01, 01);     % start of depature window
end1   = datetime(2041, 01, 01);     % end of depature window
gd1 = (start1 : days(1) : end1)';    % generates array with timesteps of 1 day

% Convert Gregian dates to Julian dates
jd1 = juliandate(gd1);

% ----------------------
% Retrieve Planetary Ephemeris Data
% ----------------------
fprintf('Fetching Mars Ephemeris from JPL Horizons...\n');
[r_mars_vec, v_mars_vec]   = get_jpl_horizons('499', start1, end1); % Mars (499)

% ----------------------
% Initial Information
%-----------------------
mu  = 1.32712440041279e11;             % Grav param of Sun (km^3/s^2)
mu_mars = 4.2828375214e4;              % Grav param of Mars (km^3/s^2)
R_mars = 3389.5;                       % Mars Radius (km)
R_sun = 696300;                        % Sun Radius (km)
%r_spacecraft                           % Radius of Spacecraft Orbit Altitude
%r_mars_vec                            % Distance from Sun to Mars


% ----------------------
% Main Code
% ----------------------

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
    
    try
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
end
