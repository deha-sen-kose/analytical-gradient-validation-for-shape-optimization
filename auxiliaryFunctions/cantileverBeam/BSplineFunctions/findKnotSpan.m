function knotSpan = findKnotSpan(knotVector, curveParameter)
%% FUNCTION findKnotSpan
%
%   This function finds the knot span index of the given parameter in a
%   given knot vector. The indexing starts from 1 in MATLAB. Thus all the
%   indexing in this framework also start from 1.
%
%   Author(s)       :  Deha Şen Köse, dehasenkose@gmail.com
%
%% Input(s):
%
%   knotVector      : Knot vector as a non-descending array.
%
%   curveParameter  : The curve parameter where the knot span index is
%                     required to be evaluated.
%
%% Output(s):
%
%   knotSpan        : The knot span index of the given curve parameter in
%                     the given knot vector.
%
%% End of function definition - Code

% Input validation.
if curveParameter < knotVector(1) || curveParameter > knotVector(end)

    error("Given curve parameter (%d) is outside of the knot vector!", curveParameter)

elseif ~issorted(knotVector)

    error("Given knot vector does not satisfy the non-descending order character!")

end

% Loop in the knot vector.
for i=1:length(knotVector)-1

    % Check if the knot vector is bigger or smaller than the knot spans.
    if curveParameter >= knotVector(i) && curveParameter < knotVector(i+1) 
        knotSpan = i;
    
    % Check for the equality especially for the last element.
    elseif curveParameter == knotVector(i+1)
        knotSpan = i;
    end

end

% % Check for input data if no index is found.
% if knotSpan == null
% 
%     error("Given curve parameter %d is outside of the knot vector %d !", curveParameter, knotVector)
% 
% end

end
%% End of code