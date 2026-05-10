function [Tdot] = DHN_ct_plant(T, q, P, Pd)
%DHN_ct_plant
%
% This function defines the continuous-time thermal model of the 4-node
% district heating network.
%
% The model calculates the temperature derivatives of the four DHN nodes.
%
% Network topology:
%
%   Node 1  ->  Node 2  ->  Node 3  ->  Node 4  ->  Node 1
%
% Physical meaning of nodes:
%
%   Node 1 = producer-side heat exchanger node
%   Node 2 = supply pipe node
%   Node 3 = consumer/load node
%   Node 4 = return pipe node
%
% Inputs:
%   T  = 4-by-1 temperature vector in degC:
%
%        T = [T1; T2; T3; T4]
%
%   q  = 4-by-1 mass-flow vector in kg/s:
%
%        q = [q12; q23; q34; q41]
%
%        q12 = mass flow from node 1 to node 2
%        q23 = mass flow from node 2 to node 3
%        q34 = mass flow from node 3 to node 4
%        q41 = mass flow from node 4 to node 1
%
%   P  = producer heat input at node 1 in Watts.
%
%   Pd = heat demand / heat extraction at node 3 in Watts.
%
% Output:
%   Tdot = 4-by-1 vector of temperature derivatives in degC/s:
%
%          Tdot = [dT1/dt; dT2/dt; dT3/dt; dT4/dt]
%
% Model form:
%
%   M*Tdot = L*T + B*P + d
%
% Therefore:
%
%   Tdot = M \ (L*T + B*P + d)
%
% where:
%
%   M = diagonal mass matrix
%   L = flow mixing matrix minus heat-loss matrix
%   B = producer heat input distribution vector
%   d = disturbance and ambient heat-loss contribution

    % -------------------------------------------------------------
    % Ambient temperatures
    % -------------------------------------------------------------

    % Define constant ambient temperatures for each node.
    %
    % Ta(1) = ambient around producer node
    % Ta(2) = ambient around supply pipe
    % Ta(3) = ambient around consumer/load node
    % Ta(4) = ambient around return pipe
    %
    % Units: degC
    Ta = [10; 10; 20; 10];

    % -------------------------------------------------------------
    % Input formatting
    % -------------------------------------------------------------

    % Force T to be a column vector.
    % This avoids dimension errors if T is passed as a row vector.
    T = T(:);

    % Force q to be a column vector.
    % This avoids dimension errors if q is passed as a row vector.
    q = q(:);

    % -------------------------------------------------------------
    % Load model parameters
    % -------------------------------------------------------------

    % Declare par as a persistent variable.
    % Persistent variables keep their value between function calls.
    %
    % This avoids recalculating physical parameters at every simulation step.
    persistent par

    % Check whether the parameter structure has already been created.
    if isempty(par)

        % If par is empty, this is the first function call.
        % Therefore, calculate the default DHN parameters.
        par = DHN_defaultParams();
    end

    % Extract the mass matrix from the parameter structure.
    % M is a 4-by-4 diagonal matrix containing nodal water masses.
    M = par.M;          % 4x4, kg

    % Extract the producer input distribution vector.
    % B maps producer power P into node 1.
    B = par.B;          % 4x1

    % Extract the demand disturbance distribution vector.
    % Bd maps heat demand Pd into node 3.
    Bd = par.Bd;        % 4x1

    % Extract heat-loss coefficients divided by cp.
    % kappa has units of kg/s and represents UA/cp for each node.
    kappa = par.kappa;  % 4x1, kg/s

    % -------------------------------------------------------------
    % Extract individual flows
    % -------------------------------------------------------------

    % Extract mass flow from node 1 to node 2.
    q12 = q(1);

    % Extract mass flow from node 2 to node 3.
    q23 = q(2);

    % Extract mass flow from node 3 to node 4.
    q34 = q(3);

    % Extract mass flow from node 4 to node 1.
    q41 = q(4);

    % -------------------------------------------------------------
    % Flow mixing matrix
    % -------------------------------------------------------------

    % Construct the flow mixing matrix for the ring topology:
    %
    %   1 -> 2 -> 3 -> 4 -> 1
    %
    % Each row represents the net mass-flow effect on one node.
    %
    % For example, row 1:
    %
    %   -q12*T1 + q41*T4
    %
    % means node 1 loses heat through outgoing flow q12 and gains heat
    % through incoming flow q41 from node 4.
    %
    % Similarly:
    %
    %   row 2 = q12*T1 - q23*T2
    %   row 3 = q23*T2 - q34*T3
    %   row 4 = q34*T3 - q41*T4
    Lflow = [ -q12,   0,    0,   q41;
               q12, -q23,   0,    0;
                0,   q23, -q34,   0;
                0,    0,   q34, -q41 ];

    % -------------------------------------------------------------
    % Heat loss and total system matrix
    % -------------------------------------------------------------

    % Subtract the heat-loss terms from the flow matrix.
    %
    % diag(kappa) creates a diagonal matrix:
    %
    %   [kappa1   0       0       0
    %      0    kappa2    0       0
    %      0      0     kappa3    0
    %      0      0       0     kappa4]
    %
    % The term -diag(kappa)*T represents heat loss from each node to ambient.
    L = Lflow - diag(kappa);

    % -------------------------------------------------------------
    % Disturbance and ambient contribution
    % -------------------------------------------------------------

    % Calculate the disturbance and ambient-temperature contribution.
    %
    % The term -Bd*Pd removes heat from node 3 because Pd is consumer demand.
    %
    % The term kappa.*Ta adds the ambient-temperature part of the heat-loss
    % expression.
    %
    % Heat loss is originally:
    %
    %   -kappa_i * (T_i - Ta_i)
    %
    % which expands to:
    %
    %   -kappa_i*T_i + kappa_i*Ta_i
    %
    % The -kappa_i*T_i part is already included in L.
    % The +kappa_i*Ta_i part is included here in d.
    d = -Bd*Pd + kappa.*Ta;

    % -------------------------------------------------------------
    % Continuous-time state equation
    % -------------------------------------------------------------

    % Calculate the temperature derivatives using:
    %
    %   Tdot = M \ (L*T + B*P + d)
    %
    % This is equivalent to:
    %
    %   Tdot = inv(M) * (L*T + B*P + d)
    %
    % but the backslash operator is numerically better than using inv(M).
    Tdot = M \ (L*T + B*P + d);

% End of main continuous-time DHN plant function.
end


function par = DHN_defaultParams()
%DHN_defaultParams
%
% This helper function defines and returns the physical parameters of the
% 4-node district heating network.
%
% The parameters include:
%
%   1) Water specific heat capacity.
%   2) Water densities at reference temperatures.
%   3) Node volumes.
%   4) Node masses.
%   5) Heat-loss coefficients.
%   6) Input distribution vectors.
%   7) Mass matrix.
%
% Output:
%   par = structure containing all model parameters.

    % -------------------------------------------------------------
    % Constants
    % -------------------------------------------------------------

    % Specific heat capacity of water.
    %
    % Units:
    %   J/(kg*K)
    %
    % This value converts thermal power terms into temperature-rate effects.
    cp = 4180;   % J/(kg*K)

    % -------------------------------------------------------------
    % Densities
    % -------------------------------------------------------------

    % Water densities corresponding to updated reference temperatures:
    %
    %   Node 1 reference temperature = 80 degC
    %   Node 2 reference temperature = 80 degC
    %   Node 3 reference temperature = 70 degC
    %   Node 4 reference temperature = 60 degC
    %
    % Density decreases as water temperature increases.
    %
    % Units:
    %   kg/m^3
    rho = [971.79; 971.79; 977.76; 980.55];   % kg/m^3

    % -------------------------------------------------------------
    % Volumes
    % -------------------------------------------------------------

    % -------------------------------------------------------------
    % V1: Producer volume
    % -------------------------------------------------------------

    % V1 represents the water volume inside the producer-side heat exchanger.
    %
    % Heat exchanger:
    %   SWEP B25T brazed plate heat exchanger
    %
    % Given:
    %   Number of plates, NoP = 140
    %   Channel volume = 0.111 dm^3 per channel per side
    %
    % Total channels:
    %   140 plates create 139 flow channels.
    %
    % Since the channels alternate between two sides:
    %   one side has approximately 70 channels,
    %   the other side has approximately 69 channels.
    %
    % Chosen producer-side channel count:
    %   69 channels
    %
    % Producer-side volume:
    %   V1 = 69 * 0.111 dm^3
    %      = 7.659 dm^3
    %      = 7.659 L
    %      = 7.659e-3 m^3
    %
    % Rounded value used:
    V1 = 7.66e-3;   % m^3

    % -------------------------------------------------------------
    % V2 and V4: Pipe volumes
    % -------------------------------------------------------------

    % Internal pipe diameter.
    %
    % Units:
    %   m
    D_pipe_in = 0.012;   % m

    % Pipe length for one pipe section.
    %
    % This length is used for both the supply pipe node and return pipe node.
    %
    % Units:
    %   m
    L_pipe = 21.0;       % m

    % Calculate supply pipe volume using:
    %
    %   V = pi * D^2 / 4 * L
    %
    % where:
    %   D = internal pipe diameter
    %   L = pipe length
    %
    % Units:
    %   m^3
    V2 = (pi * D_pipe_in^2 / 4) * L_pipe;   % m^3

    % Return pipe volume is assumed equal to supply pipe volume.
    V4 = V2;                                % m^3

    % -------------------------------------------------------------
    % V3: Consumer/load node volume
    % -------------------------------------------------------------

    % V3 represents the water volume in the consumer radiator.
    %
    % Radiator:
    %   Stelrad Classic Compact K2 600x1000
    %
    % Given water content:
    %   6.60 L
    %
    % Conversion:
    %   6.60 L = 6.60e-3 m^3
    V3 = 6.60e-3;   % m^3

    % Collect all node volumes into one vector.
    %
    % V(1) = producer volume
    % V(2) = supply pipe volume
    % V(3) = consumer/load volume
    % V(4) = return pipe volume
    V = [V1; V2; V3; V4];

    % -------------------------------------------------------------
    % Masses
    % -------------------------------------------------------------

    % Calculate nodal water masses using:
    %
    %   m_i = rho_i * V_i
    %
    % where:
    %   rho_i = density of water at node i
    %   V_i   = water volume of node i
    %
    % Units:
    %   kg
    m = rho .* V;

    % -------------------------------------------------------------
    % Heat-loss coefficients
    % -------------------------------------------------------------

    % Heat-loss coefficients are represented using:
    %
    %   kappa_i = UA_i / cp
    %
    % where:
    %   UA_i = heat transfer coefficient times area for node i
    %   cp   = specific heat capacity of water
    %
    % Units:
    %   UA_i    = W/K
    %   cp      = J/(kg*K)
    %   kappa_i = kg/s
    %
    % The heat-loss term in the temperature model is:
    %
    %   -kappa_i * (T_i - Ta_i)

    % -------------------------------------------------------------
    % K1: Producer heat-loss coefficient
    % -------------------------------------------------------------

    % The producer is modelled as an insulated rectangular body.
    %
    % Assumptions:
    %   Heat exchanger = SWEP B25T
    %   Number of plates = 140
    %   A = 0.526 m
    %   B = 0.119 m
    %   F = 0.319 m
    %   insulation thickness = 0.020 m
    %   insulation thermal conductivity = 0.035 W/(m*K)
    %   outside heat-transfer coefficient = 8 W/(m^2*K)

    % Heat exchanger dimension A.
    HX_A = 0.526;   % m

    % Heat exchanger dimension B.
    HX_B = 0.119;   % m

    % Heat exchanger dimension F.
    HX_F = 0.319;   % m

    % Calculate approximate external surface area of the rectangular body:
    %
    %   A_surf = 2*(A*B + A*F + B*F)
    %
    % Units:
    %   m^2
    A_surf = 2 * (HX_A*HX_B + HX_A*HX_F + HX_B*HX_F);   % m^2

    % Producer insulation thickness.
    t_ins_prod = 0.020;   % m

    % Thermal conductivity of producer insulation.
    lambda_ins_prod = 0.035;   % W/(m*K)

    % Outside convective heat-transfer coefficient for producer body.
    h_out_prod = 8.0;     % W/(m^2*K)

    % Calculate overall heat-transfer coefficient for the insulated producer.
    %
    % Thermal resistance per unit area is:
    %
    %   R = t_ins/lambda_ins + 1/h_out
    %
    % Therefore:
    %
    %   U = 1/R
    %
    % Units:
    %   W/(m^2*K)
    U1 = 1 / (t_ins_prod/lambda_ins_prod + 1/h_out_prod); % W/(m^2*K)

    % Calculate total producer UA value:
    %
    %   UA1 = U1 * A_surf
    %
    % Units:
    %   W/K
    UA1 = U1 * A_surf;                                       % W/K

    % -------------------------------------------------------------
    % K2 and K4: Pipe heat-loss coefficients
    % -------------------------------------------------------------

    % Pipe heat loss is calculated using cylindrical insulation resistance.
    %
    % Assumptions:
    %   insulation thickness = 0.020 m
    %   insulation thermal conductivity = 0.035 W/(m*K)
    %   outside convective coefficient = 10 W/(m^2*K)
    %   pipe outer diameter = 0.016 m
    %   pipe length = 21 m

    % Pipe outer diameter.
    d_o = 0.016;   % m

    % Pipe outer radius before insulation.
    r1 = d_o / 2;  % m

    % Pipe insulation thickness.
    t_ins_pipe = 0.020;   % m

    % Outer radius after insulation.
    r2 = r1 + t_ins_pipe; % m

    % Outside convective heat-transfer coefficient for insulated pipe.
    h_o = 10.0;    % W/(m^2*K)

    % Thermal conductivity of pipe insulation.
    lambda_ins = 0.035;   % W/(m*K)

    % Calculate insulation thermal resistance per meter:
    %
    %   Rins' = ln(r2/r1) / (2*pi*lambda_ins)
    %
    % Units:
    %   K/W per meter
    Rins_p = log(r2/r1) / (2*pi*lambda_ins); % K/W per m

    % Calculate outside convective thermal resistance per meter:
    %
    %   Rout' = 1 / (h_o * 2*pi*r2)
    %
    % Units:
    %   K/W per meter
    Rout_p = 1 / (h_o * 2*pi*r2);            % K/W per m

    % Calculate total thermal resistance per meter:
    %
    %   Rtot' = Rins' + Rout'
    %
    % Units:
    %   K/W per meter
    Rtot_p = Rins_p + Rout_p;                % K/W per m

    % Calculate total UA value for one 21 m pipe section:
    %
    %   UA_pipe = L_pipe / Rtot'
    %
    % Units:
    %   W/K
    UA_pipe = L_pipe / Rtot_p;               % W/K for one 21 m pipe

    % Assign pipe UA value to supply pipe node.
    UA2 = UA_pipe;

    % Assign pipe UA value to return pipe node.
    UA4 = UA_pipe;

    % -------------------------------------------------------------
    % K3: Consumer/load node heat-transfer coefficient
    % -------------------------------------------------------------
    % radiator heat transfer from the rated output:
    %
    %   Rated output = 1732 W at 75/65/20 degC
    %
    % Mean temperature difference:
    %
    %   DeltaT_m = (75 + 65)/2 - 20 = 50 K
    %
    % Then:
    %
    %   UA3 = Q3_rated / DeltaT_m
    % % ---- K3 Consumer (Stelrad radiator) ----
    % % Stelrad Classic Compact K2, 600x1000
    % % Rated output = 1732 W at 75/65/20
    % % Mean temperature difference:
    % % DeltaT_m = (75 + 65)/2 - 20 = 50 K
    % % UA = Q / DeltaT_m
    %
    % % Q3_rated = 1732;   % W
    % % DeltaT_m = 50;     % K
    % % UA3 = Q3_rated / DeltaT_m;   % W/K

    % In this corrected model, Pd already represents the prescribed consumer
    % heat extraction at node 3.
    %
    % Therefore, the radiator heat-transfer term:
    %
    %   UA3*(T3 - Ta3)
    %
    % is not included.
    %
    % If UA3 were included together with Pd, the consumer heat removal would
    % be counted twice. because node 3 is already a radiator.
    %
    % Node 3 heat removal is now represented only by:
    %
    %   -Bd * Pd
    %
    % Therefore, UA3 is set to zero.
    UA3 = 0;   % W/K

    % -------------------------------------------------------------
    % Collect UA values and calculate kappa
    % -------------------------------------------------------------

    % Collect heat-loss UA values for all four nodes.
    %
    % UA(1) = producer UA
    % UA(2) = supply pipe UA
    % UA(3) = consumer/load node UA
    % UA(4) = return pipe UA
    UA = [UA1; UA2; UA3; UA4];

    % Convert UA values into kappa values using:
    %
    %   kappa = UA / cp
    %
    % Units:
    %   kg/s
    kappa = UA / cp;   % kg/s

    % -------------------------------------------------------------
    % Input distribution vectors
    % -------------------------------------------------------------

    % Define producer input distribution vector.
    %
    % Producer power P enters only at node 1.
    %
    % The factor 1/cp converts power in Watts into the corresponding
    % temperature-balance term.
    B = (1/cp) * [1; 0; 0; 0];

    % Define demand disturbance distribution vector.
    %
    % Demand Pd is extracted only at node 3.
    %
    % The sign of Pd is handled in the main plant function through:
    %
    %   -Bd * Pd
    Bd = (1/cp) * [0; 0; 1; 0];

    % -------------------------------------------------------------
    % Pack parameters into structure
    % -------------------------------------------------------------

    % Store specific heat capacity.
    par.cp = cp;

    % Store density vector.
    par.rho = rho;

    % Store volume vector.
    par.V = V;

    % Store mass vector.
    par.m = m;

    % Store diagonal mass matrix.
    %
    % M = diag(m) gives:
    %
    %   [m1  0   0   0
    %    0   m2  0   0
    %    0   0   m3  0
    %    0   0   0   m4]
    par.M = diag(m);

    % Store UA heat-loss vector.
    par.UA = UA;

    % Store kappa heat-loss coefficient vector.
    par.kappa = kappa;

    % Store producer input vector.
    par.B = B;

    % Store demand disturbance vector.
    par.Bd = Bd;

% End of DHN_defaultParams helper function.
end