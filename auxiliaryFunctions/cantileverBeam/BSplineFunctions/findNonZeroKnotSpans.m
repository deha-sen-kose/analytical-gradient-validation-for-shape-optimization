function nonZeroKnotSpans = findNonZeroKnotSpans(knotVector)
%% FUNCTION findNonZeroKnotSpans
%
%   This function finds the knot span indices of the non-zero knot spans in 
%   the given knot vector. The indexing starts from 1 in MATLAB. Thus all
%   the indexing in this framework also start from 1.
%
%   Author(s)       :  Deha Şen Köse, deha.koese@tum.de
%
%% Input(s):
%
%   knotVector      : Knot vector as a non-descending array.
%
%% Output(s):
%
%   nonZeroKnotSpans: The non-zeros knot span indices of the given knot 
%                     vector. Can be a vector or a matrix depending on the 
%                     given configuration.
%
%% End of function definition - Code

% Check for input.
if ~issorted(knotVector)

    error("Given knot vector does not satisfy the non-descending order character!")

end

% Initialize output.
nonZeroKnotSpans = [];

% Loop through the knot vector.
for ii=1:length(knotVector)-1

    % Non-zero knot span criteria.
    if knotVector(ii) ~= knotVector(ii+1)
        nonZeroKnotSpans = [nonZeroKnotSpans; ii];
    end

end


end
%% End of code