function ceq = DHN_eqConFcn(X, U, data) %#ok<INUSD>
%DHN_eqConFcn
%
% This function defines the equality constraints for the NMPC controller.
%
% Equality constraints are written in the form:
%
%   ceq = 0
%
% In this DHN model, the equality constraint is used to enforce the same
% mass flow rate through all four pipe sections of the single hydraulic loop.
%
% Enforced hydraulic-loop condition:
%
%   q12 = q23 = q34 = q41
%
% This means the mass flow leaving node 1 and entering node 2 must be equal
% to the mass flow from node 2 to node 3, from node 3 to node 4, and from
% node 4 back to node 1.
%
% Inputs:
%   X    = predicted state trajectory.
%          It is not used directly in this equality constraint function.
%
%   U    = predicted input trajectory.
%          U contains all control and disturbance inputs:
%
%          U(:,1) = q12   mass flow from node 1 to node 2
%          U(:,2) = q23   mass flow from node 2 to node 3
%          U(:,3) = q34   mass flow from node 3 to node 4
%          U(:,4) = q41   mass flow from node 4 to node 1
%          U(:,5) = P     producer heat input
%          U(:,6) = Pd    heat demand / heat extracted at node 3
%
%   data = additional NMPC data structure.
%          It is not used directly in this equality constraint function.
%
% Output:
%   ceq  = equality constraint vector.
%          The NMPC optimizer enforces all elements of ceq to be zero.

    % Expected number of input columns:
    % [q12 q23 q34 q41 P Pd].
    nu_expected = 6;

    % Convert U into a clean matrix form with size Nu-by-6.
    % This helper handles row vectors, column vectors, and squeezed arrays.
    Um = local_makeUMatrix(U, nu_expected);

    % Extract the predicted mass flow from node 1 to node 2.
    q12 = Um(:,1);

    % Extract the predicted mass flow from node 2 to node 3.
    q23 = Um(:,2);

    % Extract the predicted mass flow from node 3 to node 4.
    q34 = Um(:,3);

    % Extract the predicted mass flow from node 4 to node 1.
    q41 = Um(:,4);

    % Define equality constraints for the hydraulic loop.
    %
    % First constraint:
    %   q12 - q23 = 0
    % which forces:
    %   q12 = q23
    %
    % Second constraint:
    %   q23 - q34 = 0
    % which forces:
    %   q23 = q34
    %
    % Third constraint:
    %   q34 - q41 = 0
    % which forces:
    %   q34 = q41
    %
    % Together, these three constraints enforce:
    %   q12 = q23 = q34 = q41
    ceq = [ q12 - q23;
            q23 - q34;
            q34 - q41 ];

% End of the main equality constraint function.
end

function Um = local_makeUMatrix(U, nu_expected)
%local_makeUMatrix
%
% This helper function converts the input prediction array U into a
% consistent matrix format.
%
% Expected final format:
%
%   Um = Nu-by-6
%
% where each row corresponds to one prediction step and the columns are:
%
%   [q12 q23 q34 q41 P Pd]
%
% This helper is useful because MATLAB may pass U in slightly different
% shapes depending on the prediction horizon, array dimensions, or whether
% the horizon contains only one row.

    % Remove singleton dimensions from U.
    % This is useful when U is passed with extra dimensions such as
    % 1-by-6-by-1 or similar forms.
    Um = squeeze(U);

    % Check whether U is empty.
    if isempty(Um)

        % Return an empty matrix with the expected number of columns.
        % This keeps the output format consistent even when no input rows exist.
        Um = zeros(0, nu_expected);

        % Exit the helper function immediately.
        return
    end

    % Check whether the squeezed U is a vector.
    % This can happen when there is only one prediction step.
    if isvector(Um)

        % If the vector has exactly six elements, it represents one complete
        % input row:
        % [q12 q23 q34 q41 P Pd].
        if numel(Um) == nu_expected

            % Convert the vector into a 1-by-6 row matrix.
            % This ensures that later indexing using Um(:,1), Um(:,2), etc.
            % works consistently.
            Um = reshape(Um, 1, nu_expected);

        else

            % Throw an error if the vector does not contain exactly six inputs.
            % This prevents silent mistakes caused by incorrect input dimensions.
            error('DHN_eqConFcn:BadUShape', ...
                'U has %d elements; expected %d.', numel(Um), nu_expected);
        end

    else

        % If U already has six columns, then it is already in the desired
        % Nu-by-6 format.
        if size(Um,2) == nu_expected

            % No action is needed because Um already has the correct orientation.

        % If U has six rows, it is likely in 6-by-Nu format.
        % Transpose it to obtain the desired Nu-by-6 format.
        elseif size(Um,1) == nu_expected

            % Transpose the matrix so that each row corresponds to one
            % prediction step and each column corresponds to one input.
            Um = Um.';

        else

            % Throw an error if U is neither Nu-by-6 nor 6-by-Nu.
            % This helps identify formatting problems during NMPC simulation.
            error('DHN_eqConFcn:BadUShape', ...
                'U has size %dx%d; expected Nu-by-%d or %d-by-Nu.', ...
                size(Um,1), size(Um,2), nu_expected, nu_expected);
        end
    end

% End of local_makeUMatrix helper function.
end