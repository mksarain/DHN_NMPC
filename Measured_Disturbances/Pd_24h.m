function Pd = Pd_24h(t)
%Pd_24h
%
% This function generates a 24-hour heat demand profile in Watts.
%
% The demand profile is intended for the small-scale district heating
% network model.
%
% Pd represents the thermal power extracted from the DHN at node 3.
%
% Input:
%   t  = simulation time in seconds.
%
% Output:
%   Pd = heat demand / thermal power extracted at node 3 in Watts.
%
% The profile includes:
%   1) A constant base demand.
%   2) A morning Gaussian-shaped demand peak.
%   3) An evening Gaussian-shaped demand peak.
%
% The function repeats the same demand profile every 24 hours.

% Declare persistent variables.
% Persistent variables keep their values between function calls.
% This avoids recalculating the demand profile grid every simulation step.
persistent tGrid PdGrid Ts

% Check whether the time grid has already been created.
% If tGrid is empty, this is the first time the function is being called.
if isempty(tGrid)

    % Define the sampling time for the demand profile.
    % Here, the demand profile is stored at 60-second intervals.
    Ts = 60;                  % [s]

    % Create a time grid from 0 seconds to 86400 seconds.
    % 86400 seconds = 24 hours.
    % This grid represents one complete daily demand profile.
    tGrid = 0:Ts:86400;       % [s]

    % Define the base heat demand.
    % This represents the minimum continuous thermal load.
    % 0.5e3 W = 500 W = 0.5 kW.
    Pd_base = 0.5e3;          % 0.5 kW base load

    % Define the additional morning peak demand.
    % This value is added on top of the base demand during the morning peak.
    % 1.0e3 W = 1000 W = 1.0 kW.
    Pd_morning_peak = 1.0e3;  % +1.0 kW

    % Define the additional evening peak demand.
    % This value is added on top of the base demand during the evening peak.
    % 1.5e3 W = 1500 W = 1.5 kW.
    Pd_evening_peak = 1.5e3;  % +1.5 kW

    % Define the center time of the morning demand peak.
    % 7.5 hours means 07:30.
    % Multiplication by 3600 converts hours into seconds.
    t_morning = 7.5*3600;     % 07:30

    % Define the center time of the evening demand peak.
    % 19.0 hours means 19:00.
    % Multiplication by 3600 converts hours into seconds.
    t_evening = 19.0*3600;    % 19:00

    % Define the width of the morning Gaussian peak.
    % A larger width gives a broader and smoother demand peak.
    % 1.2 hours is converted into seconds.
    w_morning = 1.2*3600;

    % Define the width of the evening Gaussian peak.
    % A larger width gives a broader and smoother demand peak.
    % 1.8 hours is converted into seconds.
    w_evening = 1.8*3600;

    % Calculate the complete 24-hour heat demand profile.
    %
    % The profile is made from:
    %   base demand
    %   + morning Gaussian peak
    %   + evening Gaussian peak
    %
    % The Gaussian expression is:
    %
    %   exp(-0.5*((t - peak_time)/peak_width)^2)
    %
    % This creates smooth peaks centered at t_morning and t_evening.
    PdGrid = Pd_base ...
        + Pd_morning_peak * exp(-0.5*((tGrid - t_morning)/w_morning).^2) ...
        + Pd_evening_peak * exp(-0.5*((tGrid - t_evening)/w_evening).^2);

    % Ensure that the demand profile never becomes negative.
    % This is a safety step, although the current formula already produces
    % positive values because it is made from positive terms.
    PdGrid = max(PdGrid, 0);

    % Force the final point of the 24-hour profile to match the first point.
    % This makes the profile periodic and avoids a discontinuity when the
    % simulation wraps from 24 hours back to 0 hours.
    PdGrid(end) = PdGrid(1);
end

% Convert the current simulation time into a time within one 24-hour day.
% mod(t,86400) makes the profile repeat every 24 hours.
% For example:
%   t = 90000 s gives tau = 3600 s,
%   meaning the second day at 01:00.
tau = mod(t, 86400);

% Interpolate the heat demand from the stored profile.
%
% 'previous' means zero-order hold interpolation:
% the function uses the most recent previous grid value.
%
% This is suitable for a sampled demand signal held constant between
% sampling instants.
%
% 'extrap' allows the function to return a value even if the query point is
% slightly outside the grid range. In normal use, tau is always between
% 0 and 86400 because of the mod operation.
Pd = interp1(tGrid, PdGrid, tau, 'previous', 'extrap');

% End of the Pd_24h function.
end