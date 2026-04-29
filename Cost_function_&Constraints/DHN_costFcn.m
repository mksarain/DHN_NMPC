function J = DHN_costFcn(X, U, e, data) %#ok<INUSD>
%DHN_costFcn
%
% Updated NMPC cost for corrected small-scale DHN model.
%
% Objective:
%   1) Track realistic nodal temperature references.
%   2) Penalize excessive mass flow.
%   3) Penalize producer power use.
%   4) Penalize sharp MV changes.
%
% Important modelling choice:
%   Pd is already the heat extracted from node 3.
%   Therefore, producer power P is NOT forced to track Pd + fixed margin.
%   Instead, P is chosen by the NMPC to maintain temperatures while avoiding
%   unnecessary power use.

nx_expected = 4;
nu_expected = 6;   % [q12 q23 q34 q41 P Pd]

Xm = local_makeXMatrix(X, nx_expected);
Um = local_makeUMatrix(U, nu_expected);

if isempty(Um)
    J = 0;
    return
end

% nlmpc usually gives:
% X : (N+1)-by-nx
% U : N-by-(nmv+nmd)
if size(Xm,1) == size(Um,1) + 1
    Xp = Xm(2:end,:);     % predicted states aligned with U rows
else
    Nmin = min(size(Xm,1), size(Um,1));
    Xp = Xm(1:Nmin,:);
    Um = Um(1:Nmin,:);
end

N = size(Um,1);

%% -------------------------------------------------------------
% Design/scaling values
% --------------------------------------------------------------

% Default temperature reference.
% This is used only if the Simulink NMPC reference signal is unavailable.
Tref_default = [80 80 70 65];     % degC

% Temperature scaling values.
% These should be close to safe/maximum operating values.
Tscale = [90 90 75 75];           % degC

% Flow scale.
% For the corrected small-scale system:
% qmax = 0.1 kg/s
% qscale = 0.05 kg/s represents a normal useful operating flow.
qscale = [0.05 0.05 0.05 0.05];   % kg/s

% Producer power scale.
% Keep this consistent with the main script:
% nlobj.MV(5).Max = 5e3;
Pmax = 5e3;                       % W

%% -------------------------------------------------------------
% Weights
% --------------------------------------------------------------

% Temperature tracking weights
aT = [100 120 100 5];

% Flow-use penalty.
% Since q is normalized by 0.05 kg/s, this weight should not be too large.
aq = [0.02 0.02 0.02 0.02];

% Producer power-use penalty
aPuse = 0.3;

% Move suppression for [q12 q23 q34 q41 P]
adu = [0.05 0.05 0.05 0.05 0.2];

%% -------------------------------------------------------------
% Reference handling
% --------------------------------------------------------------

% Use Simulink reference if available; otherwise use Tref_default.
TrefMat = local_makeRefMatrix(data, N, Tref_default);

%% -------------------------------------------------------------
% Stage cost
% --------------------------------------------------------------

J = 0;

for k = 1:N

    Tk = Xp(k,1:4);
    qk = Um(k,1:4);
    Pk = Um(k,5);

    Trefk = TrefMat(k,1:4);

    % Temperature tracking error
    eT = (Tk - Trefk) ./ Tscale;

    % Flow-use error
    eq = qk ./ qscale;

    % Producer power-use error
    ePuse = Pk / Pmax;

    J = J ...
        + sum(aT .* (eT.^2)) ...
        + sum(aq .* (eq.^2)) ...
        + aPuse * (ePuse^2);
end

%% -------------------------------------------------------------
% Move suppression penalty
% --------------------------------------------------------------

if N > 1

    dU = diff(Um(:,1:5), 1, 1);     % changes in [q12 q23 q34 q41 P]

    scale = [qscale Pmax];

    dUn = dU ./ repmat(scale, size(dU,1), 1);

    J = J + sum(sum((dUn.^2) .* repmat(adu, size(dUn,1), 1)));
end

end

%% =====================================================================
% Helper functions
% =====================================================================

function Xm = local_makeXMatrix(X, nx_expected)

Xm = squeeze(X);

if isempty(Xm)
    Xm = zeros(0, nx_expected);
    return
end

if isvector(Xm)

    if numel(Xm) == nx_expected
        Xm = reshape(Xm, 1, nx_expected);
    else
        error('DHN_costFcn:BadXShape', ...
            'X has %d elements; expected %d.', numel(Xm), nx_expected);
    end

else

    if size(Xm,2) == nx_expected
        % already N-by-4

    elseif size(Xm,1) == nx_expected
        Xm = Xm.';

    else
        error('DHN_costFcn:BadXShape', ...
            'X has size %dx%d; expected N-by-%d or %d-by-N.', ...
            size(Xm,1), size(Xm,2), nx_expected, nx_expected);
    end
end

end

function Um = local_makeUMatrix(U, nu_expected)

Um = squeeze(U);

if isempty(Um)
    Um = zeros(0, nu_expected);
    return
end

if isvector(Um)

    if numel(Um) == nu_expected
        Um = reshape(Um, 1, nu_expected);
    else
        error('DHN_costFcn:BadUShape', ...
            'U has %d elements; expected %d.', numel(Um), nu_expected);
    end

else

    if size(Um,2) == nu_expected
        % already N-by-6

    elseif size(Um,1) == nu_expected
        Um = Um.';

    else
        error('DHN_costFcn:BadUShape', ...
            'U has size %dx%d; expected N-by-%d or %d-by-N.', ...
            size(Um,1), size(Um,2), nu_expected, nu_expected);
    end
end

end

function TrefMat = local_makeRefMatrix(data, N, Tref_default)

% Default reference
TrefMat = repmat(Tref_default, N, 1);

% Try to use reference signal from Simulink NMPC block.
try
    if isstruct(data) && isfield(data,'References') && ~isempty(data.References)

        R = squeeze(data.References);

        if isvector(R)

            if numel(R) == 4
                TrefMat = repmat(R(:).', N, 1);
            end

        else

            if size(R,2) == 4

                nr = min(N, size(R,1));
                TrefMat(1:nr,:) = R(1:nr,:);

                if nr < N
                    TrefMat(nr+1:end,:) = repmat(R(nr,:), N-nr, 1);
                end

            elseif size(R,1) == 4

                R = R.';
                nr = min(N, size(R,1));
                TrefMat(1:nr,:) = R(1:nr,:);

                if nr < N
                    TrefMat(nr+1:end,:) = repmat(R(nr,:), N-nr, 1);
                end
            end
        end
    end
catch
    % If reference extraction fails, continue with default reference.
end

end