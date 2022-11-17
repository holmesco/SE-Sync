function [M1,M2,M3,V_lm] = construct_lm_data(measurements_lm,Np,Nl)
% Build landmark data matrix representation 
% Inputs:    measurements_lm - SE-sync input data structure (landmark)


% Initialize
M1 = sparse(3*Np,3*Np);
M2 = sparse(Np*3,Np-1);
M3 = sparse(Np-1,Np-1);
% Permute edges so that landmarks are grouped together. 
[~,perm] = sort(measurements_lm.edges(:,2));
edges = measurements_lm.edges(perm,:);
% Permute input data the same way
m_ij = [measurements_lm.t{:}];
m_ij = m_ij(:,perm);
weight = [measurements_lm.tau{:}];
weight_sqr = sqrt(weight(:,perm));
% Construct incidence matrix
nEdges = size(edges,1);
rows = [edges(:,1);edges(:,2)];
set = [(1:nEdges),(1:nEdges)]';
vals = [-ones(nEdges,1); ones(nEdges,1)];
B_lm = sparse(rows, set, vals, Np + Nl, nEdges);
% Construct weighted incidence matrix for landmark part
V_lm = B_lm.*weight_sqr;
% Loop over landmarks and construct matrices
curEdge = 1;
nextEdge = 1;
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
    E_l = eye(deg) - ones(deg)/deg;
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
    M1 = M1 + Vpl_Y*E_l*Vpl_Y';
    % Eq 17b
    M2 = M2 - Vpl_Y*E_l*V_p_r_l';
    % Eq 17c 
    M3 = M3 + V_p_r_l*E_l*V_p_r_l';
end


end


