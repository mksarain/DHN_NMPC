function Pd = Pd_24h(t)
%Pd_24h 24-hour heat demand profile in Watts.
% Pd represents thermal power extracted from the DHN at node 3.

persistent tGrid PdGrid Ts

if isempty(tGrid)

    Ts = 60;                  % [s]
    tGrid = 0:Ts:86400;       % [s]

    % Realistic single-radiator / small-load demand
    Pd_base = 0.5e3;          % 0.5 kW base load
    Pd_morning_peak = 1.0e3;  % +1.0 kW
    Pd_evening_peak = 1.5e3;  % +1.5 kW

    t_morning = 7.5*3600;     % 07:30
    t_evening = 19.0*3600;    % 19:00

    w_morning = 1.2*3600;
    w_evening = 1.8*3600;

    PdGrid = Pd_base ...
        + Pd_morning_peak * exp(-0.5*((tGrid - t_morning)/w_morning).^2) ...
        + Pd_evening_peak * exp(-0.5*((tGrid - t_evening)/w_evening).^2);

    PdGrid = max(PdGrid, 0);
    PdGrid(end) = PdGrid(1);
end

tau = mod(t, 86400);
Pd = interp1(tGrid, PdGrid, tau, 'previous', 'extrap');

end