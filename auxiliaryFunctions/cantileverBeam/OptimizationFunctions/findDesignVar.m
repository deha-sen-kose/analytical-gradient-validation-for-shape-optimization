function [designVarNum, desingVariablesFlat] = findDesignVar...
         (designVariablesStr, controlPointsSurface)
%% FUNCTION findDesignVar
%   This function finds and returns the design variables' numerical values 
%   from the control point coordinates. 
%
%   Author(s): Deha Şen Köse, dehasenkose@gmail.com
%
%% Input(s):
%
%   designVariablesStr      : An array containing the design variables as
%                             strings, e.g., ["x1", "y1", "x6", "y6"]. The
%                             desing variables have to be coupled in this
%                             notation! This means that x- and
%                             y-coordinates cannot be put into this array
%                             alone. Also the numbering of the degrees of
%                             freedoom or the design variable numbering,
%                             refer to the following diagram (stars are the
%                             control points.
%
%                             *5 ---------*10 --------- *15 --------- *20
%                             *4 --------- *9 --------- *14 --------- *19
%                             *3 --------- *8 --------- *13 --------- *18
%                             *2 --------- *7 --------- *12 --------- *17
%                             *1 --------- *6 --------- *11 --------- *16 
%
%   controlPointsSurface    : (nx2) array of the coordinates of the control
%                             points. n is the number of total coordinate
%                             points. The problem is simplified in 2D. The
%                             rows of this matrix should contain the both
%                             the x- and y-coordinates of the control
%                             points of the surface in the following
%                             sequence:
%                             
%                             *5 ---------*10 --------- *15 --------- *20
%                             *4 --------- *9 --------- *14 --------- *19
%                             *3 --------- *8 --------- *13 --------- *18
%                             *2 --------- *7 --------- *12 --------- *17
%                             *1 --------- *6 --------- *11 --------- *16
%
%% Output(s):
%
%   designVarNum            : (nx2) array of x- and y-coordinates of the
%                             desing variables. n is the number of desing
%                             variables divided by 2.
%
%   desingVariablesFlat     : Flattened version of the designVarNum. MATLAB
%                             expects vectors in symbolic functions and not
%                             matrices mostly. Thus the flattening helps.
%
%% End of documentation - Code

% Get the numbering of the desing variables.
strExp = regexp(designVariablesStr,'\d*','Match');

% Create a matrix for the ids extracted from the strings.
ids = zeros(length(designVariablesStr),1);

% Loop over the numberings to get the ids.
for ii= 1:length(strExp)
  if ~isempty(strExp{ii})
      ids(ii,1)=str2double(strExp{ii}(end));
  else
      ids(ii,1)=NaN;
  end
end

% Check for any errors when the designVariablesStr is provided wrong.
if any(isnan(ids))

    error("Error in finding the design variable ids!")

end

% Return the sorted unique ids of the design variables.
idsUnique = unique(ids);

% Check if the provided design variables are if in the fashion explained in
% the function description above.
if length(idsUnique) ~= length(ids)/2

    error("Design variables in x- and y-direction are not coupled! " + ...
    "If the problem statement requires this, please change this implementation!!!")

end

% Extract the coordinates (numerical values) of the design variables from
% the control point matrix.
designVarNum = zeros(length(idsUnique),2);
for ii=1:height(designVarNum)
    
    designVarNum(ii, :) = controlPointsSurface(idsUnique(ii),:);

end

% Flatten the numerical values of the design variables found above for the
% afforementioned reasons.
desingVariablesFlat = zeros(height(designVarNum)*2,1);
for ii=1:height(designVarNum)

    desingVariablesFlat(2*ii-1, 1) = designVarNum(ii,1);
    desingVariablesFlat(2*ii, 1) = designVarNum(ii,2);

end


end
%% End of code.
