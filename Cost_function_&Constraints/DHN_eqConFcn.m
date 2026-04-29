function ceq = DHN_eqConFcn(X, U, data) %#ok<INUSD>
%DHN_eqConFcn Equality constraints ceq = 0
%
% Enforce:
%   q12 = q23 = q34 = q41
%
% U contains ALL inputs:
%   U(:,1) = q12
%   U(:,2) = q23
%   U(:,3) = q34
%   U(:,4) = q41
%   U(:,5) = P
%   U(:,6) = Pd

    nu_expected = 6;
    Um = local_makeUMatrix(U, nu_expected);

    q12 = Um(:,1);
    q23 = Um(:,2);
    q34 = Um(:,3);
    q41 = Um(:,4);

    ceq = [ q12 - q23;
            q23 - q34;
            q34 - q41 ];
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
            error('DHN_eqConFcn:BadUShape', ...
                'U has %d elements; expected %d.', numel(Um), nu_expected);
        end
    else
        if size(Um,2) == nu_expected
            % already Nu-by-6
        elseif size(Um,1) == nu_expected
            Um = Um.';
        else
            error('DHN_eqConFcn:BadUShape', ...
                'U has size %dx%d; expected Nu-by-%d or %d-by-Nu.', ...
                size(Um,1), size(Um,2), nu_expected, nu_expected);
        end
    end
end