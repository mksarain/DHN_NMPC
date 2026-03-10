function [Tdot] = DHN_ct_plant(T, q, P, Pd)
%DHN_CT_PLANT 4-node DHN continuous-time thermal model
%
% Inputs :
%   T  (4x1) degC : [T1; T2; T3; T4]
%   q  (4x1) kg/s : [q12; q23; q34; q41]
%   P  (1x1) W    : producer heat input (node 1)
%   Pd (1x1) W    : disturbance power (node 3)
%
% Outputs:
%   Tdot (4x1) degC/s : temperature derivatives

    % --- Constant ambient temperatures (degC) ---
    Ta = [10; 10; 20; 10];

    % Force column vectors
    T = T(:);
    q = q(:);

    % Params (computed once)
    persistent par
    if isempty(par)
        par = DHN_defaultParams();
    end

    M     = par.M;          % 4x4, kg
    B     = par.B;          % 4x1, kg/(J/K) = 1/cp form
    Bd    = par.Bd;         % 4x1
    kappa = par.kappa;      % 4x1, kg/s
    % UA  = par.UA;         % 4x1, W/K (not needed directly here)
    % m   = par.m;          % 4x1, kg   (not needed directly here)
    % cp  = par.cp;         % scalar    (not needed directly here)

    % Flows (ring 1->2->3->4->1)
    q12 = q(1); q23 = q(2); q34 = q(3); q41 = q(4);

    Lflow = [ -q12,   0,    0,   q41;
               q12, -q23,   0,    0;
                0,   q23, -q34,   0;
                0,    0,   q34, -q41 ];

    % Heat loss term and total L
    L = Lflow - diag(kappa);

    % d(t) = -Bd*Pd + kappa.*Ta
    d = -Bd*Pd + kappa.*Ta;

    % ODE
    Tdot = M \ (L*T + B*P + d);

end


function par = DHN_defaultParams()
    % -------------------------------------------------------------
    % Constants
    % -------------------------------------------------------------
    cp = 4180;   % J/(kg*K)

    % Densities (as specified)
    rho_85C = 968.39;   % kg/m^3  -> for T1, T2
    rho_65C = 980.45;   % kg/m^3  -> for T3, T4

    % Node densities [1..4]
    rho = [rho_85C; rho_85C; rho_65C; rho_65C];

    % -------------------------------------------------------------
    % Volumes (m^3)
    % -------------------------------------------------------------

    % ---- V1 Producer volume: SWEP B25T ----
    % Given:
    % NoP = 140, channel volume = 0.111 dm^3 per channel per side
    % Total channels = 139, alternating -> one side 70, other side 69
    % Chosen producer-side volume = V_out = 69 * 0.111 dm^3 = 7.659 dm^3
    % Convert dm^3 (L) to m^3
    V1 = 7.66e-3;   % m^3

    % ---- V2 and V4 pipes ----
    % V = pi*D_in^2/4 * L
    D_pipe_in = 0.012;   % m
    L_pipe    = 3.0;     % m
    V2 = (pi * D_pipe_in^2 / 4) * L_pipe;   % m^3
    V4 = V2;                                    % m^3

    % ---- V3 Consumer volume: Stelrad Classic Compact K2 600x1000 ----
    % Water content = 6.60 L for 1000 mm length
    V3 = 6.60e-3;   % m^3

    V = [V1; V2; V3; V4];

    % Masses (kg)
    m = rho .* V;

    % -------------------------------------------------------------
    % Heat-loss coefficients
    % kappa_i = UA_i / cp   [kg/s]
    % -------------------------------------------------------------

    % ---- K1 Producer (SWEP B25T, insulated rectangular body) ----
    % Assumptions:
    %   NoP = 140
    %   A = 0.526 m, B = 0.119 m, F = 0.319 m
    %   t_ins = 0.020 m
    %   lambda_ins = 0.035 W/(m*K)
    %   h_out = 8 W/(m^2*K)

    HX_A = 0.526;   % m
    HX_B = 0.119;   % m
    HX_F = 0.319;   % m

    A_surf = 2 * (HX_A*HX_B + HX_A*HX_F + HX_B*HX_F);   % m^2

    t_ins_prod      = 0.020;   % m
    lambda_ins_prod = 0.035;   % W/(m*K)
    h_out_prod      = 8.0;     % W/(m^2*K)

    U1  = 1 / (t_ins_prod/lambda_ins_prod + 1/h_out_prod); % W/(m^2*K)
    UA1 = U1 * A_surf;                                       % W/K

    % ---- K2 and K4 Pipes ----
    % Assumptions:
    %   insulation thickness = 0.020 m
    %   lambda_ins = 0.035 W/(m*K)
    %   h_o = 10 W/(m^2*K)
    %   outer diameter d_o = 0.016 m
    %   length = 3 m

    d_o        = 0.016;   % m
    r1         = d_o / 2; % m
    t_ins_pipe = 0.020;   % m
    r2         = r1 + t_ins_pipe; % m

    h_o        = 10.0;    % W/(m^2*K)
    lambda_ins = 0.035;   % W/(m*K)

    Rins_p = log(r2/r1) / (2*pi*lambda_ins); % K/W per m
    Rout_p = 1 / (h_o * 2*pi*r2);            % K/W per m
    Rtot_p = Rins_p + Rout_p;                % K/W per m

    UA_pipe = L_pipe / Rtot_p;               % W/K for one 3 m pipe
    UA2 = UA_pipe;
    UA4 = UA_pipe;

    % ---- K3 Consumer (Stelrad radiator) ----
    % Stelrad Classic Compact K2, 600x1000
    % Rated output = 1732 W at 75/65/20
    % Mean temperature difference:
    % DeltaT_m = (75 + 65)/2 - 20 = 50 K
    % UA = Q / DeltaT_m

    Q3_rated = 1732;   % W
    DeltaT_m = 50;     % K
    UA3 = Q3_rated / DeltaT_m;   % W/K

    % Collect UA and kappa
    UA = [UA1; UA2; UA3; UA4];
    kappa = UA / cp;   % kg/s

    % -------------------------------------------------------------
    % Inputs
    % -------------------------------------------------------------
    B  = (1/cp) * [1; 0; 0; 0];
    Bd = (1/cp) * [0; 0; 1; 0];

    % -------------------------------------------------------------
    % Pack parameters
    % -------------------------------------------------------------
    par.cp    = cp;
    par.rho   = rho;
    par.V     = V;
    par.m     = m;
    par.M     = diag(m);
    par.UA    = UA;
    par.kappa = kappa;
    par.B     = B;
    par.Bd    = Bd;
end