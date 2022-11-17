function [L_s, VpBar_Y_Vst] = construct_aux_matrices(measurements_pg, ...
                          measurements_lm,V_p_t, V_lm, Np, Nl)
% Construct auxiliary matrices that are used for recovery of translation
% measurements.
% Inputs:       measurements_pg:  SE-sync measurements structure with only
%               pose graph measurements.
%               measurements_lm:  SE-sync measurements structure with only
%               pose-to-landmarks measurements.
%               V_p_t:  Weighted incidence matrix for pose graph
%               V_lm:   Weighted incidence matrix for pose-to-landmark
%               measurements
%               Np:  Number of poses
%               Nl:  Number of landmarks

% Time it
tic;
% Combine edges

edges = [measurements_pg.edges;measurements_lm.edges];
nEdgesPg = length(measurements_pg.edges);
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
V_s = [[V_p_t; sparse(Nl,nEdgesPg)],V_lm];
VpBar_Y_Vst = Vpbar_Y*V_s';
% Construct Laplacian and store
L_s = V_s*V_s';


end