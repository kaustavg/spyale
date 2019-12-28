function varargout = spyale(varargin)
%SPYALE Convert to and from CSR or CSC (Yale) sparse matrix formats.
%   The Yale sparse matrix format represents a m-by-n sparse matrix M
%   having nnz non-zero elements as three vectors:
%   A (nnz-by-1) containing the values of non-zero elements.
%   iA (m+1-by-1) containing number of non-zero elements upto but not
%       including ith row (for csr) or ith column (for csc).
%   j (nnz-by-1) containing column index for each non-zero element.
%   
%   S = spyale(A,iA,j) returns the MATLAB sparse representation given the 
%   input Yale-formatted arrays. Size of S is (length(iA)-1)-by-max(j) for 
%   CSR and max(j)-by-(length(iA)-1) for CSC.
%   
%   S = spyale(A,iA,j,m,n) specifies the size of S as m-by-n.
%   
%   [A,iA,j,m,n] = spyale(S) returns the Yale-formatted arrays for a given 
%   MATLAB sparse matrix.
%   
%   spyale(...,FORMAT) specifies whether to use CSR (default) or CSC.
%   FORMAT can be:
%   'csr'   - Use compressed sparse row format for conversions.
%   'csc'   - Use compressed sparse column format for conversions.
%   Default is 'csr'.

% Last modified by: Kaustav Gopinathan 12/28/2019

% Parse format input
if ischar(varargin{end}), fmt = varargin{end}; else fmt = 'csr'; end

% Are we converting MATLAB to Yale or Yale to MATLAB?
if issparse(varargin{1}) % MATLAB to Yale
    % Parse inputs
    S = varargin{1};
    
    % Get subscript indices
    if strcmpi(fmt,'csr')
        [j,i,A] = find(S'); % Transpose because matlab finds rows then cols
    else
        [j,i,A] = find(S); % Transpose for CSC
    end
    
    % Run-length encode to get iA
    iA = [0;cumsum(accumarray(i,1))];

    % Return Yale format arrays
    [m,n] = size(S);
    out = {A,iA,j-1,m,n};
    varargout = out(1:nargout);
    
else % Yale to MATLAB
    % Parse inputs
    [A,iA,j] = varargin{1:3};
    A = A(:); iA = iA(:); j = j(:);
    if nargin >= 5
        [m,n] = varargin{4:5};
    else
        if strcmpi(fmt,'csr')
            m = length(iA)-1; n = max(j)+1;
        else
            m = max(j)+1; n = length(iA)-1;
        end
    end
    
    % Get subscript indices using run-length decode
    [r,~,v] = find(diff(iA));
    iv = cumsum([1;v]);
	acc = zeros(iv(end)-1,1);
	acc(iv(1:end-1)) = 1;
	i = r(cumsum(acc));

    % Return sparse matrix
    if strcmpi(fmt,'csr')
        S = sparse(i,j+1,A,m,n); 
    else
        S = sparse(j+1,i,A,m,n); 
    end
    varargout = {S};
end
end

