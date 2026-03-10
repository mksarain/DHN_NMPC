function J = DHN_costFcn(X, U, e, data) %#ok<INUSD>
%DHN_costFcn Simple quadratic penalty to make P follow Pd
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

    if isempty(Um)
        J = 0;
        return
    end

    Pcol  = Um(:,5);
    Pdcol = Um(:,6);

    J = sum((Pcol - Pdcol).^2);
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
            % already Nu-by-6
        elseif size(Um,1) == nu_expected
            Um = Um.';
        else
            error('DHN_costFcn:BadUShape', ...
                'U has size %dx%d; expected Nu-by-%d or %d-by-Nu.', ...
                size(Um,1), size(Um,2), nu_expected, nu_expected);
        end
    end
end