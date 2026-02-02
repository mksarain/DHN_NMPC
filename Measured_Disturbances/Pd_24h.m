function Pd = Pd_24h(t)
%Pd_24h  24-hour heat demand profile (measured disturbance) in Watts.
%
% Notes:
% - Uses array-based profile (ZOH / piecewise constant) -> MPC-friendly.
% - Repeats every 24h via mod(t,86400).

persistent tGrid PdGrid Ts

if isempty(tGrid)
    % ---- sampling for the disturbance profile ----
    Ts = 60;                      % [s] profile resolution (same as MPC Ts )
    tGrid = 0:Ts:86400;           % [s] 24h grid (inclusive endpoint)

    % ---- Simple residential-like daily heat demand (edit these) ----
    Pd_base = 40e3;               % 40 kW base load
    Pd_morning_peak = 20e3;      % +20 kW morning peak
    Pd_evening_peak = 30e3;      % +30 kW evening peak

    % Peak timing (seconds from midnight)
    t_morning = 7.5*3600;         % 07:30
    t_evening = 19.0*3600;        % 19:00

    % Peak widths (standard deviation-like)
    w_morning = 1.2*3600;         % ~1.2 h
    w_evening = 1.8*3600;         % ~1.8 h

    % Smooth peaks (Gaussian bumps), then ZOH is applied at output
    PdGrid = Pd_base ...
           + Pd_morning_peak * exp(-0.5*((tGrid - t_morning)/w_morning).^2) ...
           + Pd_evening_peak * exp(-0.5*((tGrid - t_evening)/w_evening).^2);

    % Keep nonnegative (heat consumption cannot be negative)
    PdGrid = max(PdGrid, 0);

    % Optional: enforce last sample equals first for clean wrap
    PdGrid(end) = PdGrid(1);
end

% ---- Wrap time to 0..86400 (repeat daily) ----
tau = mod(t, 86400);

% ---- MPC-friendly lookup: piecewise constant (zero-order hold) ----
Pd = interp1(tGrid, PdGrid, tau, 'previous', 'extrap');

end
