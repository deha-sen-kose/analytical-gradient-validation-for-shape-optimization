function [masterStiffnessMtx, sensitivities] = assembleMasterStiffMtxandSensitivities(...
          knotVectorXi, knotVectorEta, polOrderXi, polOrderEta, designVariablesSym, ...
          designVariablesNum, controlPointsSurface, materialProp, symbolicFlag)
%%  Function assembleMasterStiffMtxandSensitivities
%
%   This function assebles the SYMBOLIC master FEM equations and
%   sensitivities. This is NOT used in the thesis because of efficiency
%   reasons.
%
%   Author(s)               :  Deha Şen Köse, deha.koese@tum.de
%
%%  Input(s):
%
%       knotVectorXi        : The knot vector in xi-direction defined in
%                             the parametric space. m=n+p+1 and the
%                             non-descending character must be satisfied.
%
%       knotVectorEta       : The knot vector in eta-direction defined in
%                             the parametric space. m=n+p+1 and the 
%                             non-descending character must be satisfied.
%
%   designVariablesSym      : Symbolic design variable vector. 
%
%   designVariablesNum      : Numerical values of the symbolic design
%                             variables.
%
%       polOrderXi          : Polynomial order (int) in xi-direction.
%
%       polOrderEta         : Polynomial order (int) in Eta-direction.
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
%       materialProp        : Matlab structure that contains the following
%                             material properties. .E = Young's Modulus.
%                             .nu = Poisson's Ratio. .t = Thickness of the
%                             plate.
%
%   symbolicFlag            : Flag for the symbolic calculations. Keep it
%                             true as more efficient code for numerical
%                             computations are available in this framework.
%
%%  Output(s):
%
%   masterStiffnessMtx      : The sybmolically computed master stiffness
%                             matrix.
%
%   sensitivities           : The element stiffness sensitivity obtained
%                             from symbolic derivations.
%
%%  End of function definition - Code.

% Check for the input flag.
if symbolicFlag
    % Initialize output.
    masterStiffnessMtx = sym(zeros(height(controlPointsSurface)*2));
    sensitivities = sym(cell(height(controlPointsSurface)*2));
    
    for ii=1:length(sensitivities)

        sensitivities{ii} = zeros(height(controlPointsSurface)*2);
    end
else
    % Initialize output.
    masterStiffnessMtx = zeros(height(controlPointsSurface)*2);
    sensitivities = cell(height(controlPointsSurface)*2,1);
    for ii=1:height(sensitivities)
        sensitivities{ii} = zeros(height(controlPointsSurface)*2);
    end
end

% Get the number of elements in xi and eta-directions.
xiNumElements = length(findNonZeroKnotSpans(knotVectorXi));  
etaNumElements = length(findNonZeroKnotSpans(knotVectorEta));

% Create the element index vectors for both xi- and eta-directions.
elementIndexVectorXi = zeros(xiNumElements*etaNumElements,1);
counter = 1;
for ii=1:xiNumElements
    for jj=1:etaNumElements

        elementIndexVectorXi(counter) = ii; 
        counter = counter +1;
    end
end


elementIndexVectorEta = zeros(xiNumElements*etaNumElements,1);
counter = 1;
for ii=1:xiNumElements
    for jj=1:etaNumElements

        elementIndexVectorEta(counter) = jj; 
        counter = counter +1;
    end
end

% Concatenate the element index vectors.
elementIndexVector = cat(2, elementIndexVectorXi, elementIndexVectorEta);

% Loop through the elements.
for elemIndex=1:height(elementIndexVector)

% Compute the symbolica expressions of the element stiffness matrices.
[elStiffnessMtx, elStiffnessSensitivity, elDOFTable] = computeSurfaceElementStiffnessMtx(...
         controlPointsSurface, knotVectorXi, knotVectorEta, polOrderXi, ...
         polOrderEta, elementIndexVector(elemIndex,:), designVariablesSym, designVariablesNum, materialProp, symbolicFlag);

% Perform the assembly.
masterStiffnessMtx(elDOFTable, elDOFTable) = masterStiffnessMtx(elDOFTable, elDOFTable) + elStiffnessMtx;

% Please use other functionality if no symbolic calculations are used.
if ~symbolicFlag

    for ii=1:1:length(sensitivities)

        sensitivities{ii}(elDOFTable, elDOFTable) = sensitivities{ii}(elDOFTable, elDOFTable) + elStiffnessSensitivity{ii};

    end


end



end

% Use the symbolic toolbox and get the derivatives of the stiffness
% matrices with respect to the design variables.
if symbolicFlag
    desingVariablesFlat = zeros(height(designVariablesNum)*2,1);
    for ii=1:height(designVariablesNum)
    
        desingVariablesFlat(2*ii-1, 1) = designVariablesNum(ii,1);
        desingVariablesFlat(2*ii, 1) = designVariablesNum(ii,2);
    
    end
    
    desingVariablesFlat = transpose(desingVariablesFlat);
    
    sensitivities = cell(length(designVariablesSym),1);
    for ii=1:length(designVariablesSym)
        
        sensitivities{ii} = double(subs(diff(masterStiffnessMtx, designVariablesSym(ii)),designVariablesSym,desingVariablesFlat));
    
    end
    
    masterStiffnessMtx = double(subs(masterStiffnessMtx,designVariablesSym,desingVariablesFlat));
end

end
% End of code.