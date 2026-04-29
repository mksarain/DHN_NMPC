clear; clc; close all;

%% -------------------- Tuning --------------------
Ts = 60;          % [s]
Np = 10;          % prediction horizon
Nc = 5;           % control horizon

x0    = [80; 80; 70; 65];                   % Initial temperatures (degC)
mv0 = [0.05; 0.05; 0.05; 0.05; 1.0e3];      % [q12;q23;q34;q41;P]
md0 = Pd_24h(0);                            % realistic initial demand
yref0 = x0;                                 % Temperature reference

%% -------------------- Build nlmpc object --------------------
nx  = 4;     % states
ny  = 4;     % outputs
nmv = 5;     % manipulated variables: q12,q23,q34,q41,P
nmd = 1;     % measured disturbance: Pd

% Total inputs = 6
% MV indices = 1:5, MD index = 6
nlobj = nlmpc(nx, ny, 'MV', 1:nmv, 'MD', nmv + (1:nmd));

nlobj.Ts = Ts;
nlobj.PredictionHorizon = Np;
nlobj.ControlHorizon    = Nc;

% Continuous-time model
nlobj.Model.IsContinuousTime = true;

% Model callbacks
nlobj.Model.StateFcn  = @DHN_mpcStateFcn;
nlobj.Model.OutputFcn = @DHN_mpcOutputFcn;

%% -------------------- Custom cost + constraints --------------------
nlobj.Optimization.ReplaceStandardCost = true;
nlobj.Optimization.CustomCostFcn       = @DHN_costFcn;
nlobj.Optimization.CustomEqConFcn      = @DHN_eqConFcn;

% ---------- State bounds ----------
Tmax_safe = [90; 90; 75; 75];

for i = 1:4
    nlobj.States(i).Min = 0;
    nlobj.States(i).Max = Tmax_safe(i);
    nlobj.States(i).ScaleFactor = Tmax_safe(i);
end

for i = 1:4
    nlobj.MV(i).ScaleFactor = 0.1;   % expected useful flow scale
end

% ---------- MV bounds for flows ----------
qmax = [0.1; 0.1; 0.1; 0.1];   % kg/s
qmin = 1e-3;

for i = 1:4
    nlobj.MV(i).Min = qmin;
    nlobj.MV(i).Max = qmax(i);
    nlobj.MV(i).ScaleFactor = 0.05;
end

% ---------- MV bounds for producer power ----------

nlobj.MV(5).Min = 0;
nlobj.MV(5).Max = 5e3;
nlobj.MV(5).ScaleFactor = 5e3;

%% -------------------- Solver options --------------------
nlobj.Optimization.SolverOptions = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','none', ...
    'MaxIterations', 50, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-6);

nlobj.Optimization.UseSuboptimalSolution = true;

%% -------------------- Optional validation --------------------
% try
%     validateFcns(nlobj, x0(:), mv0(:), md0);
% catch ME
%     warning('validateFcns failed: %s', ME.message);
% end

%% -------------------- Export to base workspace for Simulink --------------------
assignin('base','nlobj',  nlobj);
assignin('base','x0',     x0(:));
assignin('base','mv0',    mv0(:));
assignin('base','md0',    md0);
assignin('base','yref0',  yref0(:));

disp('Created nlobj in base workspace. Use it in the Simulink Nonlinear MPC Controller block.');

%% -------------------- Run simulation --------------------
mdl = 'DHN_Simulink';
simOut = sim(mdl);
logsout = simOut.logsout;

%% -------------------- Save results for plotting --------------------
save('DHN_sim_results.mat', 'simOut', 'logsout');

disp('Simulation completed and results saved to DHN_sim_results.mat');

%% =====================================================================
%% Local model callbacks
%% =====================================================================

function Tdot = DHN_mpcStateFcn(T, u)
% T : 4x1
% u : 6x1 -> [q12 q23 q34 q41 P Pd]

    T = T(:);
    u = u(:);

    q  = u(1:4);
    P  = u(5);
    Pd = u(6);

    Tdot = DHN_ct_plant(T, q, P, Pd);
end

function y = DHN_mpcOutputFcn(T, u) %#ok<INUSD>
% Return output as column vector (4x1)
    y = T(:);
end