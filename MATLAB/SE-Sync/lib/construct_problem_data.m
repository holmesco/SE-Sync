function problem_data = construct_problem_data(measurements)
%function problem_data = construct_problem_data(measurements)
%
% This helper function accepts a MATLAB struct containing the raw 
% measurements specifying a special Euclidean synchronization problem, and 
% returns another struct containing the data matrices required by the 
% SE-Sync algorithm.  Formally:
%
% INPUT: A MATLAB struct 'measurements' containing the following fields
% (see eq. (11) in the long-form version of the paper for details):
% edges:  An (mx2)-dimension encoding the edges in the measurement network;
%     edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform from pose i to pose j.  NB:  This indexing scheme 
%     requires that the states x_i are numbered sequentially as 
%     x_1, ... x_n.
% R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
% t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
% kappa:  An m-dimensional cell array whose kth element gives the precision
%     of the rotational part of the kth measurement. 
% tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% 
%
% OUTPUT:  A MATLAB struct 'problem_data' containing various data matrices
% that the SE-Sync algorithm requires:
%
% d:  The dimensional parameter of the special Euclidean group over which 
%     the estimation takes place (typically d = 2 or 3).
% n:  The number of states to be estimated.
% m:  The number of available measurements.
% LWtau:  The Laplacian for the translational weight graph W^tau.
% ConLap:  The connection Laplacian for the rotational measurements (see
%     eq. (15) in the paper.
% A:  The oriented incidence matrix for the underlying directed graph of
%     measurements (see eq. (7) in the paper).
% Ared:  The reduced incidence matrix obtained by removing the final 
%     row of A.
% L:  The lower-triangular Cholesky factor of the reduced translational
%     weight graph Laplacian.
% p:  Minimum degree ordering permutation vector for above Cholesky
%     Factorization.
% Omega:  The diagonal matrix of translational matrix precisions (see eq.
%     (23) in the paper).
% T:  The sparse matrix of translational observations definedin equation 
%     (24) in the paper.
% V:  The matrix of translational observations defined in equation 
%     (16) in the paper

%Given the 'measurements' struct returned by 'readG2oDataset3D', this
%function constructs and returns the data matrices defining the pose-graph
%relaxation problem

% Copyright (C) 2016 by David M. Rosen

% If landmark flag is available, split the measurement into pose graph and
% landmakr measurements
if isfield(measurements,'lmFlag') && ~isempty(measurements.lmFlag)
    tic;
    lmFlag = measurements.lmFlag;
    % Split the data 
    fnames = fieldnames(measurements);
    for iName = 1:length(fnames)
        fname = fnames{iName};
        measurements_pg.(fname) = measurements.(fname)(~lmFlag,:);
        measurements_lm.(fname) = measurements.(fname)(lmFlag,:);
    end
    % Get number of poses and landmarks 
    Np = max(measurements_pg.edges(:,1));        % Number of Poses
    Nl = max(measurements_lm.edges(:,2))-Np;     % Number of Landmarks
    
    % Pose graph data matrices
    [M1pg, M2pg, M3pg, V_p_t] = construct_pg_data(measurements_pg, Np);
    % Landmark graph data matrices
    [M1lm, M2lm, M3lm, V_lm] = construct_lm_data(measurements_lm, Np, Nl);
    
    % Create implicit data object. This avoids creating a dense Q matrix, and
    % allows us to exploit sparsity when the number of poses is large.
    Q_bt = lm_data(M1pg+M1lm,M2pg+M2lm,M3pg+M3lm);
    fprintf('Done constructing landmark-marginalized data matrices in %g seconds\n', toc);
    
    % Construct connection Laplacian for the rotational measurements
    tic();
    problem_data.ConLap = construct_connection_Laplacian(measurements_pg);
    t = toc();
    fprintf('Constructed rotational connection Laplacian in %g seconds\n', t);
    
    % Auxiliary matrix generation
    [L_s, VpBar_Y_Vst] = construct_aux_matrices(measurements_pg, measurements_lm,...
                          V_p_t, V_lm, Np, Nl);
    fprintf('Done constructing auxiliary data matrices: \t%g\n', toc);
    
    % Collect data in problem_data struct
    problem_data.d = length(measurements_pg.t{1});
    problem_data.n = max(max(measurements_pg.edges));
    problem_data.m = size(measurements_pg.edges, 1);
    problem_data.Q_bt = Q_bt;
    problem_data.L_s = L_s;
    problem_data.VpBar_Y_Vst = VpBar_Y_Vst;

else
%     measurements_pg = measurements;
    % Set additional variables
    problem_data.d = length(measurements.t{1});
    problem_data.n = max(max(measurements.edges));
    problem_data.m = size(measurements.edges, 1);

    % Construct connection Laplacian for the rotational measurements
    tic();
    problem_data.ConLap = construct_connection_Laplacian(measurements);
    t = toc();
    fprintf('Constructed rotational connection Laplacian in %g seconds\n', t);


    % Construct the oriented incidence matrix for the underlying directed graph
    % of measurements
    tic();
    problem_data.A = construct_incidence_matrix(measurements);
    t = toc();
    fprintf('Constructed oriented incidence matrix in %g seconds\n', t);

    % Construct the reduced oriented incidence matrix
    problem_data.Ared = problem_data.A(1:problem_data.n-1, :);

    tic();
    [T, Omega] = construct_translational_matrices(measurements);
    V = construct_V_matrix(measurements);
    t = toc();
    fprintf('Constructed translational observation and measurement precision matrices in %g seconds\n', t);

    problem_data.T = T;
    problem_data.Omega = Omega;
    problem_data.V = V;


    % Construct the Laplacian for the translational weight graph
    tic();
    LWtau = problem_data.A * problem_data.Omega * problem_data.A';
    t = toc();
    fprintf('Constructed Laplacian for the translational weight graph in %g seconds\n', t);
    problem_data.LWtau = LWtau;

    % Construct the Cholesky factor for the reduced translational weight graph
    % Laplacian 
    % CTH modified to minimize fill in (return permutation matrix), because
    % running into memory issues for large problems
    % NOTE: (for later) use the permuted version of the Cholesky decomposition 
    % tic();
    % [problem_data.L,flag,problem_data.p] = chol(LWtau(1:end-1, 1:end-1), 'lower','vector');
    % [problem_data.L] = chol(LWtau(1:end-1, 1:end-1), 'lower');
    % t = toc();
    % fprintf('Computed lower-triangular factor of reduced translational weight graph Laplacian in %g seconds\n', t);

    % Cache a couple of various useful products
    fprintf('Caching additional product matrices ... \n');

    problem_data.sqrt_Omega = spdiags(sqrt(spdiags(Omega)), 0, problem_data.m, problem_data.m);
    problem_data.sqrt_Omega_AredT = problem_data.sqrt_Omega * problem_data.Ared';
    problem_data.sqrt_Omega_T = problem_data.sqrt_Omega * problem_data.T;

end

