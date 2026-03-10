clear; clc;

%% -------------------- Tuning --------------------
Ts = 60;          % [s]
Np = 10;          % prediction horizon
Nc = 5;           % control horizon

x0    = [85; 85; 65; 65];       % T0 (degC)
mv0   = [1; 1; 1; 1; 40e3];     % [q12;q23;q34;q41;P]
md0   = 40e3;                   % Pd (W)
yref0 = x0;                     % temperature reference (degC)

%% -------------------- Build nlmpc object --------------------
nx  = 4;     % states
ny  = 4;     % outputs
nmv = 5;     % q's + P
nmd = 1;     % Pd
nu  = nmv + nmd;   %#ok<NASGU>

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
nlobj.Optimization.CustomIneqConFcn    = @DHN_ineqConFcn;
nlobj.Optimization.CustomEqConFcn      = @DHN_eqConFcn;

%% -------------------- Solver options --------------------
nlobj.Optimization.SolverOptions = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','none', ...
    'MaxIterations', 50, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-6);

nlobj.Optimization.UseSuboptimalSolution = true;

%% -------------------- Validate (kept commented) --------------------
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

%% -------------------- Open model + run simulation --------------------
mdl = 'DHN_Simulink';

open_system(mdl);
set_param(mdl,'SimulationCommand','start');