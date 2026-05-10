clear;
clc;
close all;

%% -------------------- Tuning --------------------

% Define the NMPC sampling time.
% The controller computes a new control action every Ts seconds.
Ts = 60;          % [s]

% Define the prediction horizon.
% The NMPC predicts the system response over Np future sampling steps.
% With Ts = 60 s and Np = 10, the prediction horizon is 10 minutes.
Np = 10;          % prediction horizon

% Define the control horizon.
% The optimizer is allowed to independently adjust the control moves over
% the first Nc future steps. After that, the moves are usually held constant.
% With Ts = 60 s and Nc = 5, the independent control horizon is 5 minutes.
Nc = 5;           % control horizon

% Define the initial state vector.
% These are the initial temperatures of the four DHN nodes.
%
% x0 = [T1; T2; T3; T4]
%
% T1 = producer-side heat exchanger temperature
% T2 = supply pipe temperature
% T3 = consumer/load node temperature
% T4 = return pipe temperature
x0 = [80; 80; 70; 65];                   % Initial temperatures (degC)

% Define the initial manipulated-variable vector.
%
% mv0 = [q12; q23; q34; q41; P]
%
% q12 = mass flow from node 1 to node 2
% q23 = mass flow from node 2 to node 3
% q34 = mass flow from node 3 to node 4
% q41 = mass flow from node 4 to node 1
% P   = producer heat input
%
% The initial flows are all equal because the network is a single
% hydraulically connected loop.
mv0 = [0.05; 0.05; 0.05; 0.05; 1.0e3];      % [q12;q23;q34;q41;P]

% Define the initial measured disturbance.
% Pd_24h(0) returns the heat demand at time t = 0 seconds.
% Pd represents the thermal power extracted at node 3.
md0 = Pd_24h(0);                            % realistic initial demand

% Define the initial output reference.
% Since the outputs are the four temperatures, the reference is also a
% four-element temperature vector.
%
% Here, the initial reference is set equal to the initial state.
yref0 = x0;                                 % Temperature reference

%% -------------------- Build nlmpc object --------------------

% Define the number of states.
% The four states are:
%
% x = [T1; T2; T3; T4]
nx = 4;     % states

% Define the number of outputs.
% The output function returns all four temperatures.
%
% y = [T1; T2; T3; T4]
ny = 4;     % outputs

% Define the number of manipulated variables.
%
% Manipulated variables are variables predicted by the NMPC optimizer.
%
% MV = [q12; q23; q34; q41; P]
nmv = 5;    % manipulated variables: q12,q23,q34,q41,P

% Define the number of measured disturbances.
%
% Pd is treated as a measured disturbance because it is imposed on the
% system and not directly controlled by the NMPC.
nmd = 1;    % measured disturbance: Pd

% Total input vector used by the state function:
%
% u = [q12; q23; q34; q41; P; Pd]
%
% Total number of inputs = nmv + nmd = 6.
%
% MV indices = 1:5
% MD index   = 6
%
% Create the nonlinear MPC object.
% The syntax specifies:
%   nx = number of states
%   ny = number of outputs
%   MV = indices of manipulated variables in the input vector
%   MD = indices of measured disturbances in the input vector
nlobj = nlmpc(nx, ny, 'MV', 1:nmv, 'MD', nmv + (1:nmd));

% Assign the controller sampling time.
nlobj.Ts = Ts;

% Assign the prediction horizon.
nlobj.PredictionHorizon = Np;

% Assign the control horizon.
nlobj.ControlHorizon = Nc;

% Specify that the plant model is continuous-time.
% Then the controller automatically discretizes the continuous-time prediction 
% model using the implicit trapezoidal rule with the sample time nlobj.Ts. 
% MathWorks states that nonlinear MPC is a discrete-time controller, so 
% continuous-time state functions are automatically discretized internally for prediction.
nlobj.Model.IsContinuousTime = true;

% Assign the state function used by the NMPC.
% This function receives the current state and input vector and returns Tdot.
nlobj.Model.StateFcn = @DHN_mpcStateFcn;

% Assign the output function used by the NMPC.
% This function maps the state vector to the output vector.
nlobj.Model.OutputFcn = @DHN_mpcOutputFcn;

%% -------------------- Custom cost + constraints --------------------

% Replace MATLAB's standard quadratic tracking cost with a custom cost.
% The custom cost is defined in DHN_costFcn.
nlobj.Optimization.ReplaceStandardCost = true;

% Assign the custom NMPC cost function.
% This function includes temperature tracking, flow penalty, power-use
% penalty, and move suppression.
nlobj.Optimization.CustomCostFcn = @DHN_costFcn;

% Assign the custom equality constraint function.
% This function enforces:
%
% q12 = q23 = q34 = q41
%
% so that all pipe sections have the same loop mass flow.
nlobj.Optimization.CustomEqConFcn = @DHN_eqConFcn;

%% -------------------- State bounds --------------------

% Define safe maximum temperatures for each node.
%
% Tmax_safe = [T1_max; T2_max; T3_max; T4_max]
%
% These limits prevent the optimizer from allowing unrealistic or unsafe
% temperatures during prediction.
Tmax_safe = [90; 90; 75; 75];

% Apply lower bounds, upper bounds, and scale factors to all four states.
for i = 1:4

    % Set minimum allowed temperature for state i.
    % A value of 0 degC is used as a broad physical lower bound.
    nlobj.States(i).Min = 0;

    % Set maximum allowed temperature for state i.
    nlobj.States(i).Max = Tmax_safe(i);

    % Set the scale factor for state i.
    % Scale factors improve numerical conditioning of the optimization.
    nlobj.States(i).ScaleFactor = Tmax_safe(i);
end

% Assign initial scale factors for the first four manipulated variables.
% These are the mass-flow inputs.
%
% A value of 0.1 kg/s corresponds to the expected upper useful flow range.
for i = 1:4

    % Set the flow scale factor.
    nlobj.MV(i).ScaleFactor = 0.1;   % expected useful flow scale
end

%% -------------------- MV bounds for flows --------------------

% Define maximum mass-flow values for the four pipe sections.
%
% qmax = [q12_max; q23_max; q34_max; q41_max]
%
% The same maximum value is used for all four flows.
qmax = [0.1; 0.1; 0.1; 0.1];   % kg/s

% Define minimum mass-flow value.
% A small positive lower bound avoids exactly zero flow, which can cause
% poor controllability and numerical issues in thermal transport models.
qmin = 1e-3;

% Apply flow bounds and flow scale factors to the first four MVs.
for i = 1:4

    % Set minimum allowed mass flow for pipe section i.
    nlobj.MV(i).Min = qmin;

    % Set maximum allowed mass flow for pipe section i.
    nlobj.MV(i).Max = qmax(i);

    % Set the scale factor used by the optimizer for flow variable i.
    % A value of 0.05 kg/s represents a typical useful operating flow.
    nlobj.MV(i).ScaleFactor = 0.05;
end

%% -------------------- MV bounds for producer power --------------------

% Set minimum producer power.
% The producer cannot supply negative heat in this model.
nlobj.MV(5).Min = 0;

% Set maximum producer power.
% This value should be consistent with the power scale used in DHN_costFcn.
nlobj.MV(5).Max = 5e3;

% Set producer power scale factor.
% This improves numerical conditioning because P is much larger in magnitude
% than the flow variables.
nlobj.MV(5).ScaleFactor = 5e3;

%% -------------------- Solver options --------------------

% Define optimization solver settings for fmincon.
%
% The nonlinear MPC toolbox uses an optimization solver internally.
% Here, fmincon with the SQP algorithm is selected.

nlobj.Optimization.SolverOptions = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...              % Use Sequential Quadratic Programming algorithm.
    'Display', 'none', ...               % Suppress solver iteration output in the Command Window.
    'MaxIterations', 50, ...             % Limit the number of solver iterations at each NMPC step.
    'OptimalityTolerance', 1e-3, ...     % Stop when first-order optimality is below this tolerance.
    'StepTolerance', 1e-6);              % Stop when optimization step size becomes very small.

% Allow the NMPC to use a suboptimal solution if the solver does not fully
% converge within the iteration limit.
nlobj.Optimization.UseSuboptimalSolution = true;

%% -------------------- Optional validation --------------------

% Validate the NMPC functions before running the Simulink simulation.
%
% validateFcns checks whether:
%   1) the state function runs correctly,
%   2) the output function runs correctly,
%   3) the custom cost function runs correctly,
%   4) the custom constraint function runs correctly,
%   5) dimensions of states, inputs, and outputs are consistent.
try

    % Validate all NMPC callback functions using the initial condition,
    % initial manipulated variables, and initial measured disturbance.
    validateFcns(nlobj, x0(:), mv0(:), md0);

catch ME

    % If validation fails, show a warning but do not stop the script.
    % This allows the user to inspect the warning and still continue if needed.
    warning('validateFcns failed: %s', ME.message);
end

%% -------------------- Export to base workspace for Simulink --------------------

% Export the NMPC object to the MATLAB base workspace.
% The Simulink Nonlinear MPC Controller block can access this variable.
assignin('base','nlobj', nlobj);

% Export the initial state vector to the base workspace.
assignin('base','x0', x0(:));

% Export the initial manipulated-variable vector to the base workspace.
assignin('base','mv0', mv0(:));

% Export the initial measured disturbance to the base workspace.
assignin('base','md0', md0);

% Export the initial reference vector to the base workspace.
assignin('base','yref0', yref0(:));

% Display a message confirming that the NMPC object has been created.
disp('Created nlobj in base workspace. Use it in the Simulink Nonlinear MPC Controller block.');

%% -------------------- Run simulation --------------------

% Define the Simulink model name.
% The model file should be named DHN_Simulink.slx.
mdl = 'DHN_Simulink';

% Run the Simulink model.
% The simulation output is stored in simOut.
simOut = sim(mdl);

% Extract logged simulation signals from the simulation output.
% The Simulink model should have signal logging enabled and should store
% relevant signals in logsout.
logsout = simOut.logsout;

%% -------------------- Save results for plotting --------------------

% Save the simulation output and logged signals to a MAT-file.
% This file can later be loaded by a separate plotting script.
save('DHN_sim_results.mat', 'simOut', 'logsout');

% Display a message confirming that the simulation has finished.
disp('Simulation completed and results saved to DHN_sim_results.mat');

%% =====================================================================
%% Local model callbacks
%% =====================================================================

function Tdot = DHN_mpcStateFcn(T, u)
%DHN_mpcStateFcn
%
% This local callback function is used by the nonlinear MPC object.
%
% It maps the NMPC input vector into the arguments required by the
% continuous-time DHN plant model.
%
% Inputs:
%   T = current state vector:
%
%       T = [T1; T2; T3; T4]
%
%   u = complete NMPC input vector:
%
%       u = [q12; q23; q34; q41; P; Pd]
%
% Outputs:
%   Tdot = temperature derivative vector:
%
%          Tdot = [dT1/dt; dT2/dt; dT3/dt; dT4/dt]

    % Force the state vector to be a column vector.
    % This prevents dimension errors if MATLAB passes T as a row vector.
    T = T(:);

    % Force the input vector to be a column vector.
    % This prevents dimension errors if MATLAB passes u as a row vector.
    u = u(:);

    % Extract the four mass-flow inputs from the complete input vector.
    %
    % q = [q12; q23; q34; q41]
    q = u(1:4);

    % Extract producer heat input.
    % P is the fifth input.
    P = u(5);

    % Extract heat demand / disturbance.
    % Pd is the sixth input.
    Pd = u(6);

    % Call the continuous-time DHN plant model.
    % The plant model calculates the temperature derivatives.
    Tdot = DHN_ct_plant(T, q, P, Pd);

% End of DHN_mpcStateFcn.
end

function y = DHN_mpcOutputFcn(T, u) %#ok<INUSD>
%DHN_mpcOutputFcn
%
% This local callback function defines the NMPC output equation.
%
% In this model, all four states are directly used as outputs.
%
% Inputs:
%   T = current state vector:
%
%       T = [T1; T2; T3; T4]
%
%   u = complete input vector.
%       It is not used in this output function.
%
% Output:
%   y = output vector:
%
%       y = [T1; T2; T3; T4]

    % Return the temperature state vector as a column vector.
    % This ensures that the output has the required 4-by-1 format.
    y = T(:);

% End of DHN_mpcOutputFcn.
end