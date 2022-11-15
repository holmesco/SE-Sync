function [problem_data] = landmark_wrapper(measurements)
% This helper function adds landmark support to SE-sync algorithm without
% naively representing landmarks with poses.

% If landmark flag is available, split the measurement into pose graph and
% landmakr measurements
if isfield(measurements,'lmFlag') && ~isempty(measurements.lmFlag)
    lmFlag = measurements.lmFlag;
    % Split the data 
    fnames = fieldnames(measurements);
    for iName = 1:length(fnames)
        fname = fnames{iName};
        measurements_pg.(fname) = measurements.(fname)(~lmFlag,:);
        measurements_lm.(fname) = measurements.(fname)(lmFlag,:);
    end
    % Get problem data for pose graph measurements
    problem_data = construct_problem_data(measurements_pg);
    % Add problem data for landmark measurements
    lmDataTime = tic;
    Np = max(measurements_pg.edges(:,1));        % Number of Poses
    Nl = max(measurements_lm.edges(:,2))-Np;     % Number of Landmarks
    problem_data.lmMat = construct_lm_data(measurements_lm,Np, Nl);
    fprintf('Computed Landmark Information:\t%g secs\n',toc(lmDataTime))
else
    % If no landmark data then construct normally 
    problem_data = construct_problem_data(measurements);
end
    
end

