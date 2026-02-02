clear; clc;

%% -------------------- Tuning --------------------
Ts = 60;          % [s]
Np = 10;          % prediction horizon
Nc = 5;           % control horizon

x0    = [60; 55; 40; 45];              % T0 (degC)
mv0   = [5; 5; 5; 5; 5e4];             % [q12;q23;q34;q41;P]
md0   = 1e4;                           % Pd (W)
yref0 = x0;                            % temperature reference (degC)

%% -------------------- Build nlmpc object --------------------
nx  = 4;     % states
ny  = 4;     % outputs (all T)
nmv = 5;     % q's + P
nmd = 1;     % Pd

%  (MV channels + MD channels)
nlobj = nlmpc(nx, ny, 'MV', 1:nmv, 'MD', nmv + (1:nmd));   % MV=1..5, MD=6

nlobj.Ts = Ts;
nlobj.PredictionHorizon = Np;
nlobj.ControlHorizon    = Nc;

% Continuous-time model
nlobj.Model.IsContinuousTime = true;

% State function: Tdot = f(T, mv, md)
nlobj.Model.StateFcn  = @DHN_mpcStateFcn;
nlobj.Model.OutputFcn = @DHN_mpcOutputFcn;

%% -------------------- Custom cost + constraints --------------------
nlobj.Optimization.ReplaceStandardCost = true;
nlobj.Optimization.CustomCostFcn       = @DHN_costFcn;
nlobj.Optimization.CustomIneqConFcn    = @DHN_ineqConFcn;

%% -------------------- Solver options (keep simple) --------------------
nlobj.Optimization.SolverOptions = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','none', ...
    'MaxIterations', 50, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-6);

nlobj.Optimization.UseSuboptimalSolution = true;

%% -------------------- Validate --------------------
% x0v   = x0(:);
% mv0v  = mv0(:).';        % 1x5
% md0v  = md0;             % 1x1
% yrefv = yref0(:).';      % 1x4
% 
% try
%     validateFcns(nlobj, x0v, mv0v, yrefv, md0v);
% catch
%     try
%         validateFcns(nlobj, x0v, mv0v, md0v);
%     catch ME
%         warning('validateFcns failed: %s', ME.message);
%     end
% end

%% -------------------- Export to base workspace for Simulink --------------------
assignin('base','nlobj',  nlobj);
assignin('base','x0',     x0(:));
assignin('base','mv0',    mv0(:));
assignin('base','md0',    md0);
assignin('base','yref0',  yref0(:));

disp('Created nlobj in base workspace. Use it in the Simulink Nonlinear MPC Controller block.');

%% =====================================================================
%% Local model callbacks (required by nlobj)
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
    y = T(:).';   % 1x4
end


%% -------------------- Open model + run simulation --------------------
mdl = 'DHN_Simulink';

% opens Simulink window 
open_system(mdl);   

% Start simulation
set_param(mdl,'SimulationCommand','start');

