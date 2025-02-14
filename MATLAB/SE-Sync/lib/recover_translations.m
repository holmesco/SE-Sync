function t = recover_translations(R, problem_data)
% function t = recover_translations(R, problem_data)
%
% Given a dense d x dN block matrix R whose (dxd)-block elements are
% rotations comprising an optimal solution to the reduced rotation-only
% maximum-likelihood estimation, this function returns a dense (d x n)
% matrix t whose columns are the oorresponding optimal translations

% Copyright (C) 2016 by David M. Rosen


% We have that vec(t) = -(LWtau^pinv V \otimes I_3) vec(R)  (1)
%
% Using the fact that vec(AB) = (B' \otimes I) vec(A),
% we see that (1) is equivalent to
%
% t = - R * V' * LWtau^pinv
% Check if we are running landmark based SLAM
if isfield(problem_data,'Q_bt')
    % Load stored RHS
    b = problem_data.VpBar_Y_Vst'*R';
%     t = lsqminnorm(problem_data.L_s,b)';
    % Invert the laplacian matrix to get all translations.
    t = (problem_data.L_s \ b)';
else
    t = - (problem_data.LWtau \ (problem_data.V * R'))';
end


end

