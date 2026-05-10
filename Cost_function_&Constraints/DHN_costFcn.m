function J = DHN_costFcn(X, U, e, data) %#ok<INUSD>
%DHN_costFcn
%
% This function defines the custom objective/cost function for the NMPC.
%
% Inputs:
%   X    = predicted state trajectory from the NMPC solver.
%          For this DHN model, the states are:
%          X = [T1 T2 T3 T4]
%
%   U    = predicted input trajectory from the NMPC solver.
%          For this DHN model, the inputs are:
%          U = [q12 q23 q34 q41 P Pd]
%
%   e    = slack variable used by MATLAB NMPC for soft constraints.
%          It is not used in this cost function.
%
%   data = additional NMPC data structure.
%          It may contain the reference signal from the Simulink NMPC block.
%
% Output:
%   J    = scalar total cost value minimized by the NMPC optimizer.
%
% Updated NMPC cost for corrected small-scale DHN model.
%
% Objective:
%   1) Track realistic nodal temperature references.
%   2) Penalize excessive mass flow.
%   3) Penalize producer power use.
%   4) Penalize sharp manipulated-variable changes.
%
% Important modelling choice:
%   Pd is already the heat extracted from node 3.
%   Therefore, producer power P is NOT forced to track Pd + fixed margin.
%   Instead, P is chosen by the NMPC to maintain temperatures while avoiding
%   unnecessary power use.

% Expected number of thermal states: [T1 T2 T3 T4].
nx_expected = 4;

% Expected number of input columns: [q12 q23 q34 q41 P Pd].
% The first five are manipulated variables or control-related variables.
% Pd is the measured disturbance.
nu_expected = 6;

% Convert X into a clean matrix form with size N-by-4 or (N+1)-by-4.
% This helper handles row vectors, column vectors, and squeezed arrays.
Xm = local_makeXMatrix(X, nx_expected);

% Convert U into a clean matrix form with size N-by-6.
% This helper handles row vectors, column vectors, and squeezed arrays.
Um = local_makeUMatrix(U, nu_expected);

% If the input trajectory is empty, there is no prediction horizon to cost.
if isempty(Um)

    % Return zero cost because no control moves are available.
    J = 0;

    % Exit the function immediately.
    return
end

% MATLAB nlmpc usually provides:
% X : (N+1)-by-nx matrix because it includes the current state and N future states.
% U : N-by-(nmv+nmd) matrix because there are N future input moves.
%
% Therefore, if X has one more row than U, the first row of X corresponds
% to the current state and should be skipped for stage-cost alignment.
if size(Xm,1) == size(Um,1) + 1

    % Use predicted future states from row 2 onward.
    % This aligns Xp(k,:) with Um(k,:).
    Xp = Xm(2:end,:);

else

    % If X and U are not in the standard (N+1)-by-nx and N-by-nu format,
    % use the minimum available number of rows to avoid indexing errors.
    Nmin = min(size(Xm,1), size(Um,1));

    % Keep only the common number of state rows.
    Xp = Xm(1:Nmin,:);

    % Keep only the common number of input rows.
    Um = Um(1:Nmin,:);
end

% Number of prediction steps available in the input trajectory.
N = size(Um,1);

%% -------------------------------------------------------------
% Design/scaling values
% --------------------------------------------------------------

% Default temperature reference.
% This is used only if the Simulink NMPC reference signal is unavailable.
% The values correspond to realistic nominal temperatures of the four nodes.
Tref_default = [80 80 70 65];     % degC

% Temperature scaling values.
% These are used to normalize the temperature tracking error.
% They should be close to safe or maximum operating temperatures.
Tscale = [90 90 75 75];           % degC

% Flow scaling values.
% These normalize the mass-flow terms in the cost function.
% For the corrected small-scale system, qmax is assumed to be 0.1 kg/s.
% A qscale value of 0.05 kg/s represents a normal useful operating flow.
qscale = [0.05 0.05 0.05 0.05];   % kg/s

% Producer power scaling value.
% This should be consistent with the MV upper bound in the main NMPC script.
% Example:
% nlobj.MV(5).Max = 5e3;
Pmax = 5e3;                       % W

%% -------------------------------------------------------------
% Weights
% --------------------------------------------------------------

% Temperature tracking weights for [T1 T2 T3 T4].
% Higher values make the controller prioritize temperature tracking more.
% T2 is weighted slightly higher here because supply-side temperature
% tracking is usually important in a DHN.
% T4 is weighted lightly because return temperature may be less critical.
aT = [100 120 100 5];

% Flow-use penalty weights for [q12 q23 q34 q41].
% These discourage unnecessarily high circulation flow.
% Since q is already normalized by 0.05 kg/s, the weight is kept small.
aq = [0.02 0.02 0.02 0.02];

% Producer power-use penalty.
% This discourages unnecessary heat production.
% A larger value would make the controller more energy-saving but may
% reduce temperature tracking performance.
aPuse = 0.3;

% Move suppression weights for [q12 q23 q34 q41 P].
% These penalize rapid changes in flow and producer power.
% The producer power move penalty is higher than the flow move penalty.
adu = [0.05 0.05 0.05 0.05 0.2];

%% -------------------------------------------------------------
% Reference handling
% --------------------------------------------------------------

% Build the reference matrix for the full prediction horizon.
% If data.References is available from Simulink, it will be used.
% Otherwise, the constant default reference Tref_default is used.
TrefMat = local_makeRefMatrix(data, N, Tref_default);

%% -------------------------------------------------------------
% Stage cost
% --------------------------------------------------------------

% Initialize total cost before accumulating stage costs.
J = 0;

% Loop over each prediction step.
for k = 1:N

    % Extract predicted nodal temperatures at prediction step k.
    % Tk = [T1 T2 T3 T4].
    Tk = Xp(k,1:4);

    % Extract predicted mass flows at prediction step k.
    % qk = [q12 q23 q34 q41].
    qk = Um(k,1:4);

    % Extract predicted producer power at prediction step k.
    % Pk is the fifth input column.
    Pk = Um(k,5);

    % Extract the temperature reference for prediction step k.
    % Trefk = [T1_ref T2_ref T3_ref T4_ref].
    Trefk = TrefMat(k,1:4);

    % Calculate normalized temperature tracking error.
    % Normalization prevents large numerical temperature values from
    % dominating the optimization only because of their units.
    eT = (Tk - Trefk) ./ Tscale;

    % Calculate normalized flow-use term.
    % This penalizes the absolute use of circulation flow.
    eq = qk ./ qscale;

    % Calculate normalized producer power-use term.
    % This penalizes the absolute use of heating power.
    ePuse = Pk / Pmax;

    % Add the stage cost at prediction step k.
    % The stage cost contains:
    %   1) weighted squared temperature tracking error,
    %   2) weighted squared flow-use penalty,
    %   3) weighted squared producer power-use penalty.
    J = J ...
        + sum(aT .* (eT.^2)) ...
        + sum(aq .* (eq.^2)) ...
        + aPuse * (ePuse^2);
end

%% -------------------------------------------------------------
% Move suppression penalty
% --------------------------------------------------------------

% Move suppression can only be calculated if there is more than one
% prediction step, because diff() needs at least two rows.
if N > 1

    % Calculate the step-to-step changes in manipulated variables.
    % Only the first five columns are used:
    % [q12 q23 q34 q41 P].
    % Pd is excluded because it is a measured disturbance, not a manipulated
    % variable chosen by the controller.
    dU = diff(Um(:,1:5), 1, 1);

    % Define scaling values for the manipulated-variable changes.
    % The first four scales correspond to mass-flow changes.
    % The fifth scale corresponds to producer-power changes.
    scale = [qscale Pmax];

    % Normalize the manipulated-variable changes.
    % repmat repeats the scale row so that its size matches dU.
    dUn = dU ./ repmat(scale, size(dU,1), 1);

    % Add the move suppression cost to the total cost.
    % Each squared normalized move is weighted by adu.
    J = J + sum(sum((dUn.^2) .* repmat(adu, size(dUn,1), 1)));
end

% End of the main cost function.
end

%% =====================================================================
% Helper functions
% =====================================================================

function Xm = local_makeXMatrix(X, nx_expected)
%local_makeXMatrix
%
% This helper function converts the state prediction array X into a
% consistent matrix format.
%
% Expected final format:
%   Xm = N-by-4 or (N+1)-by-4
%
% where the columns are:
%   [T1 T2 T3 T4]

% Remove singleton dimensions from X.
% This is useful because MATLAB nlmpc may pass arrays with extra dimensions.
Xm = squeeze(X);

% Check whether X is empty.
if isempty(Xm)

    % Return an empty matrix with the expected number of columns.
    Xm = zeros(0, nx_expected);

    % Exit the helper function.
    return
end

% Check whether the squeezed X is a vector.
if isvector(Xm)

    % If the vector has exactly four elements, it represents one state row.
    if numel(Xm) == nx_expected

        % Convert the vector into a 1-by-4 row matrix.
        Xm = reshape(Xm, 1, nx_expected);

    else

        % Throw an error if the vector does not contain exactly four states.
        error('DHN_costFcn:BadXShape', ...
            'X has %d elements; expected %d.', numel(Xm), nx_expected);
    end

else

    % If X already has four columns, it is already in N-by-4 format.
    if size(Xm,2) == nx_expected

        % No action is needed because Xm already has the correct orientation.

    % If X has four rows, it is likely in 4-by-N format.
    % Transpose it to obtain N-by-4 format.
    elseif size(Xm,1) == nx_expected

        % Transpose the matrix so that each row is one prediction step.
        Xm = Xm.';

    else

        % Throw an error if X is neither N-by-4 nor 4-by-N.
        error('DHN_costFcn:BadXShape', ...
            'X has size %dx%d; expected N-by-%d or %d-by-N.', ...
            size(Xm,1), size(Xm,2), nx_expected, nx_expected);
    end
end

% End of local_makeXMatrix helper function.
end

function Um = local_makeUMatrix(U, nu_expected)
%local_makeUMatrix
%
% This helper function converts the input prediction array U into a
% consistent matrix format.
%
% Expected final format:
%   Um = N-by-6
%
% where the columns are:
%   [q12 q23 q34 q41 P Pd]

% Remove singleton dimensions from U.
% This is useful because MATLAB nlmpc may pass arrays with extra dimensions.
Um = squeeze(U);

% Check whether U is empty.
if isempty(Um)

    % Return an empty matrix with the expected number of columns.
    Um = zeros(0, nu_expected);

    % Exit the helper function.
    return
end

% Check whether the squeezed U is a vector.
if isvector(Um)

    % If the vector has exactly six elements, it represents one input row.
    if numel(Um) == nu_expected

        % Convert the vector into a 1-by-6 row matrix.
        Um = reshape(Um, 1, nu_expected);

    else

        % Throw an error if the vector does not contain exactly six inputs.
        error('DHN_costFcn:BadUShape', ...
            'U has %d elements; expected %d.', numel(Um), nu_expected);
    end

else

    % If U already has six columns, it is already in N-by-6 format.
    if size(Um,2) == nu_expected

        % No action is needed because Um already has the correct orientation.

    % If U has six rows, it is likely in 6-by-N format.
    % Transpose it to obtain N-by-6 format.
    elseif size(Um,1) == nu_expected

        % Transpose the matrix so that each row is one prediction step.
        Um = Um.';

    else

        % Throw an error if U is neither N-by-6 nor 6-by-N.
        error('DHN_costFcn:BadUShape', ...
            'U has size %dx%d; expected N-by-%d or %d-by-N.', ...
            size(Um,1), size(Um,2), nu_expected, nu_expected);
    end
end

% End of local_makeUMatrix helper function.
end

function TrefMat = local_makeRefMatrix(data, N, Tref_default)
%local_makeRefMatrix
%
% This helper function builds the reference temperature matrix used by the
% cost function.
%
% Expected final format:
%   TrefMat = N-by-4
%
% where each row contains:
%   [T1_ref T2_ref T3_ref T4_ref]
%
% If Simulink provides a reference through data.References, that reference
% is used. Otherwise, the default reference is repeated over the horizon.

% Create a default reference matrix by repeating Tref_default for all N steps.
TrefMat = repmat(Tref_default, N, 1);

% Use try-catch so that the cost function does not fail if the reference
% signal is missing or has an unexpected format.
try

    % Check whether data is a structure, contains the field 'References',
    % and that the reference field is not empty.
    if isstruct(data) && isfield(data,'References') && ~isempty(data.References)

        % Remove singleton dimensions from the reference signal.
        R = squeeze(data.References);

        % Check whether the reference signal is a vector.
        if isvector(R)

            % If the vector has four elements, it represents one temperature
            % reference vector [T1_ref T2_ref T3_ref T4_ref].
            if numel(R) == 4

                % Repeat the single reference vector over the full horizon.
                TrefMat = repmat(R(:).', N, 1);
            end

        else

            % If R has four columns, it is already in N-by-4 style format.
            if size(R,2) == 4

                % Use the smaller of:
                %   N              = required prediction length,
                %   size(R,1)      = available reference rows.
                nr = min(N, size(R,1));

                % Copy the available reference rows into TrefMat.
                TrefMat(1:nr,:) = R(1:nr,:);

                % If fewer reference rows are available than prediction steps,
                % repeat the last available reference row for the remaining steps.
                if nr < N

                    % Fill the remaining prediction steps using the last
                    % available reference row.
                    TrefMat(nr+1:end,:) = repmat(R(nr,:), N-nr, 1);
                end

            % If R has four rows, it is likely in 4-by-N format.
            % It must be transposed to become N-by-4.
            elseif size(R,1) == 4

                % Transpose R so that each row corresponds to one prediction step.
                R = R.';

                % Use the smaller of the required and available number of rows.
                nr = min(N, size(R,1));

                % Copy the available reference rows into TrefMat.
                TrefMat(1:nr,:) = R(1:nr,:);

                % If fewer reference rows are available than prediction steps,
                % repeat the last available reference row for the remaining steps.
                if nr < N

                    % Fill the remaining prediction steps using the last
                    % available reference row.
                    TrefMat(nr+1:end,:) = repmat(R(nr,:), N-nr, 1);
                end
            end
        end
    end

catch

    % If reference extraction fails for any reason, the function continues
    % using Tref_default. This avoids stopping the NMPC simulation because
    % of reference formatting issues.
end

% End of local_makeRefMatrix helper function.
end