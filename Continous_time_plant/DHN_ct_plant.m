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
%   Tdot (4x1)    : temperature derivatives

    % --- Constant ambient temperatures (degC) ---
    Ta = [10; 10; 20; 10];

    % Force column vectors
    T  = T(:);
    q  = q(:);

    % Params (computed once)
    persistent par
    if isempty(par)
        par = DHN_defaultParams();
    end

    M     = par.M;          % 4x4
    B     = par.B;          % 4x1
    Bd    = par.Bd;         % 4x1
    kappa = par.kappa;      % 4x1
    UA    = par.UA;         % 4x1 (W/C)
    m     = par.m;          % 4x1 (kg)
    cp    = par.cp;         % scalar

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

    % ODE: Tdot
    Tdot = M \ (L*T + B*P + d);

end


function par = DHN_defaultParams()
    % Constants
    cp = 4180;                 % J/(kg*C)
    rho_60C = 983.13;
    rho_30C = 995.71;

    % Node densities [1..4]
    rho = [rho_60C; rho_60C; rho_60C; rho_30C];

    % Convection & insulation
    h_o = 10;                  % W/(m2*C)
    t_ins = 0.10;              % m
    lambda_ins = 0.35;         % W/(m*C)

    % ---- Volumes ----
    % V1 tank
    D_in = 60; H_tank = 10; fill_fraction = 0.9;
    H_water = fill_fraction * H_tank;
    V1 = (pi*D_in^2/4) * H_water;

    % V2 and V4 pipes
    D_pipe_in = 0.7; L_pipe = 20;
    V2 = (pi*D_pipe_in^2/4) * L_pipe;
    V4 = V2;

    % V3 tube bundle 
    N_tubes = 20;
    d_in = 0.01;         % m  
    L_tube = 0.5;        % m
    V3 = N_tubes * (pi*d_in^2/4) * L_tube;

    V = [V1; V2; V3; V4];
    m = rho .* V;

    % ---- Heat-loss (kappa = UA/cp) ----
    % Node 1 tank: planar insulation + outside convection
    A_side = pi*D_in*H_water;
    A_bottom = pi*D_in^2/4;
    A1 = A_side + A_bottom;

    Rins1 = t_ins / (lambda_ins * A1);
    Rout1 = 1 / (h_o * A1);
    UA1 = 1 / (Rins1 + Rout1);

    % Node 2/4 pipes: log insulation + outside convection
    r1 = D_pipe_in/2;
    r2 = r1 + t_ins;
    Rins_p = log(r2/r1) / (2*pi*lambda_ins);
    Rout_p = 1 / (h_o * 2*pi*r2);
    UA_pipe = L_pipe / (Rins_p + Rout_p);
    UA2 = UA_pipe;
    UA4 = UA_pipe;

    % Node 3 radiator to room: UA = (h_conv + h_rad)*A
    h_conv = 10; h_rad = 5.5; h_total = h_conv + h_rad;
    hx_L = 1.0; hx_W = 0.5; hx_H = 0.15;
    A3 = 2*(hx_L*hx_W + hx_L*hx_H + hx_W*hx_H);
    UA3 = h_total * A3;

    UA = [UA1; UA2; UA3; UA4];
    kappa = UA / cp;

    % Inputs
    B  = (1/cp) * [1; 0; 0; 0];
    Bd = (1/cp) * [0; 0; 1; 0];

    % Pack
    par.cp    = cp;
    par.V     = V;
    par.m     = m;
    par.M     = diag(m);
    par.UA    = UA;
    par.kappa = kappa;
    par.B     = B;
    par.Bd    = Bd;
end
