function [problem_data] = construct_lm_data(measurements_lm, measurements_pg)
% Build landmark data matrix representation 
% Inputs:    measurements_lm - SE-sync input data structure (landmark)
%            measurements_pg - SE-sync input data structure (pose graph)
tic;
% Get number of poses and landmarks 
Np = max(measurements_pg.edges(:,1));        % Number of Poses
Nl = max(measurements_lm.edges(:,2))-Np;     % Number of Landmarks

% Initialize
M1 = sparse(3*Np,3*Np);
M2 = sparse(Np*3,Np-1);
M3 = sparse(Np-1,Np-1);

%% LANDMARK TRANSLATION DATA
% Permute edges so that landmarks are grouped together. 
[~,perm] = sort(measurements_lm.edges(:,2));
edges = measurements_lm.edges(perm,:);
% Permute input data the same way
m_ij = [measurements_lm.t{perm}];
weight_sqr = sqrt([measurements_lm.tau{perm}]);
% Construct incidence matrix
nEdges = size(edges,1);
rows = [edges(:,1);edges(:,2)];
set = [(1:nEdges),(1:nEdges)]';
vals = [-ones(nEdges,1); ones(nEdges,1)];
B_lm = sparse(rows, set, vals, Np + Nl, nEdges);
% Construct weighted incidence matrix for landmark part
V_lm = B_lm.*weight_sqr;
% "Edge Leaving" part of incidence matrix
B_lm_1 = sparse(rows(1:nEdges), set(1:nEdges), vals(1:nEdges), Np, nEdges);
V_lm_1 = -B_lm_1.*weight_sqr;
% Loop over landmarks and construct matrices
curEdge = 1;
nextEdge = 1;
E = cell(Nl,1);
for l = 1 : Nl
    % Get the set of edges/cols in incidence matrix corresponding to this
    % landmark
    while nextEdge <= nEdges && edges(nextEdge,2) == edges(curEdge,2)
        nextEdge = nextEdge + 1;
    end
    set = curEdge : nextEdge-1;
    curEdge = nextEdge;
    % get degree
    deg = length(set);
    % Construct E matrix (schur compelment) component
    E{l} = sparse(eye(deg) - ones(deg)/deg);
    % Measurement data (weighted)
    m_set_wght = m_ij(:,set).*weight_sqr(set);
    % Pose part of landmark incidence matrix and reduced incidence
    % matrix (weighted)
    V_p_r_l = V_lm(2:Np,set);
    % Compute kronecker product as sparesly to avoid costly implementation
    % Original Code:
    %     tmp = mat2cell(m_ij(:,set),3,ones(1,deg));
    %     Y{l} = sparse(blkdiag(tmp{:}));
    %     V_p_l = V_lm(1:Np, set);
    %     Vpl_Y = kron(V_p_l,speye(3))*Y{l};
    rows = ones(3,1)*((edges(set,1)-1)*3)' + [1;2;3]*ones(1,deg);
    cols = ones(3,1)*(1:deg);
    Vpl_Y = -sparse(rows(:),cols(:),m_set_wght(:),Np*3,deg);
    % Eq 17a
    M1 = M1 + Vpl_Y*E{l}*Vpl_Y';
    % Eq 17b
    M2 = M2 - Vpl_Y*E{l}*V_p_r_l';
    % Eq 17c 
    M3 = M3 + V_p_r_l*E{l}*V_p_r_l';
end

%% POSE GRAPH TRANSLATION DATA
% Get data from input struct
edges = measurements_pg.edges;
t_ij = measurements_pg.t;
weight = [measurements_pg.tau{:}]';
% Construct incidence matrix
nEdges = size(edges,1);
rows = [edges(:,1);edges(:,2)];
set = [(1:nEdges),(1:nEdges)]';
vals = [-ones(nEdges,1); ones(nEdges,1)];
B_pg = sparse(rows, set, vals, Np, nEdges);
% Weighted incidence matrix and reduced version
V_p_t = B_pg.*sqrt(weight)';
V_p_t_r = V_p_t(2:Np,:);
% "Edge Leaving" part of incidence matrix
B_pg_1 = sparse(rows(1:nEdges), set(1:nEdges), vals(1:nEdges), Np, nEdges);
V_p_t_1 = -B_pg_1.*sqrt(weight)';
% Generate measurement matrix
Y_t = sparse(blkdiag(t_ij{:}));
% Add data to matrices
Vpt1Y = kron(V_p_t_1, speye(3))*Y_t;
M1 = M1 + Vpt1Y*Vpt1Y';
M2 = M2 + Vpt1Y*V_p_t_r';
M3 = M3 + V_p_t_r*V_p_t_r';
%% COMBINE INTO CLASS.
% Create implicit data object. This avoids creating a dense Q matrix, and
% allows us to exploit sparsity when the number of poses is large.
problem_data.Q_bt = lm_data(M1,M2,M3);
fprintf('Done constructing landmark-marginalized data matrices in %g seconds\n', toc);

%% Compute Matrices for translation recovery
tic;
% Combine edges
edges = [measurements_pg.edges;measurements_lm.edges];
nEdgesAll = size(edges,1);
% Combine data and weights
tm = [[measurements_pg.t{:}],[measurements_lm.t{:}]];
weights_sqr = sqrt([[measurements_pg.tau{:}],[measurements_lm.tau{:}]]);
tm_wght = tm.*weights_sqr;
% Build V_p_bar*Y_s sparsely
rows = ones(3,1)*((edges(:,1)-1)*3)' + [1;2;3]*ones(1,nEdgesAll);
cols = ones(3,1)*(1:nEdgesAll);
Vpbar_Y = sparse(rows(:),cols(:),tm_wght(:),Np*3,nEdgesAll);
% Make weighted incidence matrix
V_s = [[V_p_t; sparse(Nl,nEdges)],V_lm];
VpBar_Y_Vst = Vpbar_Y*V_s';
% Construct Laplacian and store
problem_data.L_s = V_s*V_s';
problem_data.VpBar_Y_Vst = VpBar_Y_Vst;
fprintf('Done constructing auxiliary data matrices: \t%g\n', toc);

end


