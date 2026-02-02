function J = DHN_costFcn(X, U, e, data) %#ok<INUSD>


    % ---------- weights  ----------
    wT   = [1 1 1 1];
    wU   = [1e-3 1e-3 1e-3 1e-3 1e-8];
    wDU  = [1e-2 1e-2 1e-2 1e-2 1e-7];
    wDef = 1e-6;

    ny_expected  = 4;
    nmv_expected = 5;

    % ---- Make e a 2-D matrix with ny columns 
    E = squeeze(e);
    if ~ismatrix(E), E = reshape(E, size(E,1), []); end
    if size(E,2) ~= ny_expected && size(E,1) == ny_expected
        E = E.';  % make it Ne-by-4
    end

    % ---- Make U a 2-D matrix with nmv columns
    Um = squeeze(U);
    if ~ismatrix(Um), Um = reshape(Um, size(Um,1), []); end
    if size(Um,2) ~= nmv_expected && size(Um,1) == nmv_expected
        Um = Um.'; % make it Nu-by-5
    end

    % ---------- tracking term ----------
    if isempty(E)
        J_track = 0;
    else
        ny = size(E,2);
        wT_use = fitWeights(wT, ny);
        J_track = sum(sum(E.^2 .* reshape(wT_use,1,ny)));
    end

    % ---------- MV magnitude + move terms ----------
    if isempty(Um)
        J_mv = 0;
        J_dmv = 0;
        J_def = 0;
        J = J_track;
        return
    end

    nU = size(Um,2);
    wU_use  = fitWeights(wU,  nU);
    wDU_use = fitWeights(wDU, nU);

    % Column-wise energy (no implicit expansion problems)
    J_mv  = sum( sum(Um.^2,1)  .* (wU_use(:).') );

    dU = diff(Um,1,1);
    if isempty(dU)
        J_dmv = 0;
    else
        J_dmv = sum( sum(dU.^2,1) .* (wDU_use(:).') );
    end

    % ---------- Soft encourage P >= Pd (only if column 5 exists) ----------
    if size(Um,2) >= 5
        Nu = size(Um,1);
        PdVec = local_getPdVec(data, Nu);
        Pcol  = Um(:,5);
        deficit = max(0, PdVec - Pcol);
        J_def = wDef * sum(deficit.^2);
    else
        J_def = 0;
    end

    J = J_track + J_mv + J_dmv + J_def;
end

function w = fitWeights(wIn, n)
    wIn = wIn(:);
    if isempty(wIn)
        w = ones(n,1);
        return
    end
    if numel(wIn) >= n
        w = wIn(1:n);
    else
        w = [wIn; wIn(end)*ones(n-numel(wIn),1)];
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
        PdVec = Pd(1:min(end,N));
        if numel(PdVec) < N
            PdVec(end+1:N,1) = PdVec(end);
        end
    end
end
