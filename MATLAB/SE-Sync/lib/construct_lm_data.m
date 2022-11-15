function [landmark_data] = construct_lm_data(meas, Np, Nl)
% Build landmark data matrix representation 
% Inputs:    meas - measurement structure (same as SE-Sync
%            Np   - Number of poses
%            Nl   - Number of landmarks
% Initialize
M1 = sparse(3*Np,3*Np);
M2 = sparse(Np*3,Np-1);
M3 = sparse(Np-1,Np-1);
% Permute edges so that landmarks are grouped together. 
[~,perm] = sort(meas.edges(:,2));
edges = meas.edges(perm,:);
% Permute input data the same way
t_ij = [meas.t{perm}];
weight = [meas.tau{perm}]';
% Construct incidence matrix
nEdges = size(edges,1);
rows = [edges(:,1);edges(:,2)];
cols = [(1:nEdges),(1:nEdges)]';
vals = [-ones(nEdges,1); ones(nEdges,1)];
B = sparse(rows, cols, vals, Np + Nl, nEdges);
% Construct weighted incidence matrix for landmark part
V_lm = B.*sqrt(weight)';
% Loop over landmarks and construct matrices
curEdge = 1;
nextEdge = 1;
E = cell(Nl,1);
for l = 1 : Nl
    % Get the set of columns in incidence matrix corresponding to this
    % landmark
    while nextEdge <= nEdges && edges(nextEdge,2) == edges(curEdge,2)
        nextEdge = nextEdge + 1;
    end
    cols = curEdge : nextEdge-1;
    curEdge = nextEdge;
    % get degree
    deg = length(cols);
    % Construct E matrix (schur compelment) component
    E{l} = sparse(eye(deg) - ones(deg)/deg);
    % Construct measurement matrix
    tmp = mat2cell(t_ij(:,cols),3,ones(1,deg));
    Y_l = sparse(blkdiag(tmp{:}));
    % Pose part of landmark incidence matrix and reduced incidence
    % matrix (weighted)
    V_p_l = V_lm(1:Np, cols);
    V_p_r_l = V_lm(2:Np,cols);
    % Eq 17a - Should be able to make this part much faster: don't use
    % kron to place the matrices
    M1 = M1 + kron(V_p_l,eye(3))*Y_l*E{l}*Y_l'*kron(V_p_l,eye(3))';
    % Eq 17b
    M2 = M2 - kron(V_p_l,eye(3))*Y_l*E{l}*V_p_r_l';
    % Eq 17c 
    M3 = M3 + V_p_r_l*E{l}*V_p_r_l';
end

% Create implicit Q object. This avoids creating a dense Q matrix, and
% allows us to exploit sparsity when the number of poses is large.
landmark_data = lm_data(M1,M2,M3);

end


