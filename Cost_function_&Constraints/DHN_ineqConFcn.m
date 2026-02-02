function c = DHN_ineqConFcn(X, U, e, data) %#ok<INUSD>
%DHN_ineqConFcn Robust inequality constraints c<=0

    q_eps = 1e-3;           % kg/s
    ENFORCE_P_GE_PD = true; % set false if infeasible

    nmv = 5;

    % ---- make U = Nu-by-nmv
    Um = U;
    if size(Um,2) == nmv
        % ok
    elseif size(Um,1) == nmv
        Um = Um.'; % transpose to Nu-by-nmv
    end

    Nu = size(Um,1);

    % q >= q_eps  ->  q_eps - q <= 0
    qMat = Um(:,1:4);
    c_q  = q_eps - qMat(:);

    if ENFORCE_P_GE_PD
        PdVec = local_getPdVec(data, Nu);
        Pcol  = Um(:,5);
        c_p   = PdVec - Pcol;   % Pd - P <= 0
        c     = [c_q; c_p(:)];
    else
        c = c_q;
    end
end

function PdVec = local_getPdVec(data, N)
    Pd = 0;

    if isstruct(data)
        if isfield(data,'MD'), Pd = data.MD;
        elseif isfield(data,'MeasuredDisturbances'), Pd = data.MeasuredDisturbances;
        elseif isfield(data,'md'), Pd = data.md;
        elseif isfield(data,'UserData') && isstruct(data.UserData) && isfield(data.UserData,'Pd')
            Pd = data.UserData.Pd;
        end
    end

    Pd = Pd(:);
    if isempty(Pd), Pd = 0; end

    if numel(Pd) == 1
        PdVec = repmat(Pd, N, 1);
    else
        if numel(Pd) < N
            PdVec = [Pd; Pd(end)*ones(N-numel(Pd),1)];
        else
            PdVec = Pd(1:N);
        end
    end
end
