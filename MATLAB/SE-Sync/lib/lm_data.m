classdef lm_data
    % This class defines an implicit version of the data matrix Q
    properties
        M1      % Constructed Matrix definition
        M2      %    Q_bt = M1 + M2*pinv(M3)*M2'
        M3      % 
        R       % Thin QR decomposition of M3 matrix
        L       % Cholesky factor of M3 projection matrix
        sz      % Size of the problem
    end
    methods
        % get size
        function [n,m] = size(self)
            [n,m] = size(self.Q_r);
        end
        % Contructor, stores matrices and performs decomposition of M3
        function self = lm_data(M1,M2,M3)
            % Make composite matrices sparse
            self.M1 = M1;
            self.M2 = M2;
            self.M3 = M3;
            % Problem size
            self.sz = size(M1);
            % Cholesky factorization of M3 matrix
            self.L = chol(self.M3);
        end
        % Overload matrix multiply
        function [Y] = mtimes(self, X)
            % Check size
            [n,nVec] = size(X);
            if n ~= self.sz
                error('Incorrect dims for matrix multiplication')
            end
            Y_1 = self.M1 * X;
            % For the rest of the elements, implement columnwize
            y_2 = cell(1,nVec);
            for i = 1:nVec
                x1 = self.M2'*X(:,i);
                x2 = self.L' \ x1;
                x3 = self.L \ x2;
                y_2{i} = self.M2*x3;
            end
            % Assemble result
            Y_2 = cat(2,y_2{:});
            Y = Y_1 - Y_2;
        end
        % Make function handle (for getting eigenvalues). Allow auxiliary
        % matrix A which is added to function
        function [f] = getFunc(self,A)
            f = @(x) self.mtimes(x) + A*x; 
        end

        % Explicitly construct Data Matrix - tolerance required on
        % pseudoinverse. There are numerical issues when the number of 
        % landmarks is large.
        function Q = double(self)
            Q_bt = self.M1 + - (self.M2 * (self.M3 \ self.M2')) ;
            % Combine Matrices
            Q = Q_bt;
        end
    end
end