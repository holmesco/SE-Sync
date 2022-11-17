function [M1,M2,M3,V_p_t] = construct_pg_data(measurements_pg, Np)
% Construct data matrices for the pose graph part of the problem
% Inputs    measurements_pg - SE-sync input data structure (pose graph)
%           Np - number of poses
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
M1 = Vpt1Y*Vpt1Y';
M2 = Vpt1Y*V_p_t_r';
M3 = V_p_t_r*V_p_t_r';

end