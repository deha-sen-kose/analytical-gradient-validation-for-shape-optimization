function [conIneqValue, conEqValue, conIneqGrad, conEqGrad] = ...
          volumeConstraint(volumeFunHandle, knotVectorXi, knotVectorEta,...
          polOrderXi, polOrderEta, controlPointsSurfaceOrig, propMaterial, ...
          designVariables, designVariableDOF, initialVolume, percentage, ...
          numCPy, symtryC, symtrynum, boundStyle, p1, p2)
%% Function volumeConstraint
% 
%   This function converts a volume calculation function handle to a
%   constraint form as required in Matlab optimizers. Returns only one
%   constraint, which can be chosen as;
%   "UpperBound" in the form of: 
%   current volume  <= initial volume * percentage, where percentage > 1. 
%   "LowerBound" in the form of:
%   current volume  >= initial volume * percentage, where percentage < 1.
%   "Equality" in the form of:
%   current volume  = initial volume * percentage.
%   "LowerAndUpper in the form of:
%   current volume  <= initial volume * p1, where p1 > 1 and
%   current volume  >= initial volume * p2, where p2 < 1.
%
%   Author(s)       : Deha Şen Köse, deha.koese@tum.de
%
%% Input(s):
%
%   volumeFunHandle : Function handle for the volume calculation. This
%                     function requires the following input:
%                     
%                     knotVectorXi :
%                     The knot vector in xi-direction defined in the 
%                     parametric space. m=n+p+1 and the non-descending 
%                     character must be satisfied.
%
%                     knotVectorEta:
%                     The knot vector in eta-direction defined in the 
%                     parametric space. m=n+p+1 and the non-descending 
%                     character must be satisfied.
%
%                     polOrderXi   :
%                     Polynomial order (int) in xi-direction.
%
%                     polOrderEta   :
%                     Polynomial order (int) in xi-direction.
%
%                     controlPointsSurface:
%                     (nx2) array of the coordinates of the control points.
%                     n is the number of total coordinate points. The 
%                     problem is simplified in 2D. The rows of this matrix
%                     should contain the both the x- and y-coordinates of 
%                     the control points of the surface in the following
%                     sequence:
%                             *5 ---------*10 --------- *15 --------- *20
%                             *4 --------- *9 --------- *14 --------- *19
%                             *3 --------- *8 --------- *13 --------- *18
%                             *2 --------- *7 --------- *12 --------- *17
%                             *1 --------- *6 --------- *11 --------- *16
%
%                     propMaterial:
%                     Matlab structure that contains the following material
%                     properties. .E = Young's Modulus. .nu = Poisson's 
%                     Ratio. .t = Thickness of the plate.
%
%   designVariables : The numeric values of the design variables in the
%                     curret optimization iteration. Function can also be 
%                     used without the optimization.
%
%  designVariableDOF: The corresponding DOFs of the design variables.
%
%   initialVolume   : The volume of the initial design. The volume
%                     constraints can be defined through the initial
%                     volume.
%
%   percentage      : Percentage of the initial volume that is required to
%                     build the volume constraint. 
%
%   numCPy          : The number of control points in y-direction.
%
%   symtryC         : A boolean for the symmetrical control point linking. 
%                     If true, a symmetry is assumed w.r.t. the mid x-line 
%                     using the symtrynum input.
%
%   symtrynum       : A vector containing the -X control point (meaning 
%                     that they are below the assumed symmetry line for the
%                     design variables) numbers that are assumed to be 
%                     optimized symmetrical to the corresponding above 
%                     control points. Can be left [] if symtryC is false.
%
%   p1, p2          : percentages for the "LowerAndUpper" volumetric 
%                     constraint.
%
%% Output(s):
%
%   conIneqValue    : The inequality constraint value (ineq. volume value).
%
%   conEqValue      : The equality constraint value (eq. volume value).
%                     This is an empty array in this code because there is
%                     no equality constraints at the moment.
%
%   conIneqGrad     : Inequality constraint gradient.
%
%   conEqGrad       : Equality constraint gradient.
%                     This is an empty array in this code because there is
%                     no equality constraints at the moment.
%
%   boundStyle      : A string from the following "UpperBound",
%                     "LowerBound", "Equality" and "LowerAndUpper".
%                     These terms are explained in the beginning of the
%                     documentation.
%
%% End of function definition - Code

% Get the number of design variables in the optimization process.
numDesignVar = length(designVariables);
controlPointsSurface = controlPointsSurfaceOrig;
% Loop through the design variables to update the surface control points 
% matrix respecting the control polygons. 
for ii=1:numDesignVar

    % Get the DOF of the current desing variable.
    dof = designVariableDOF(ii);
    % Find the row of that design variable in the surface control point
    % matrix.
    r = ceil(dof/2);

    % Find the column of that desing variable in the surface control point
    % matrix.
    if mod(dof,2) == 0
        c = 2;
    else
        c = 1;
    end
    
    % Update the control point matrix with the current desing variables.
    controlPointsSurface(r,c) = designVariables(ii);

end

% If requested, apply the symmetrical parameter linking.
if symtryC == true

    controlPointsSurface = symmetryFun(controlPointsSurfaceOrig,...
        controlPointsSurface, numCPy,symtrynum);
end

% Compute the volume and the volume gradient using the function handle.
volume = volumeFunHandle(knotVectorXi, ...
          knotVectorEta, polOrderXi, polOrderEta, controlPointsSurface,...
          propMaterial);

% Using the constraint type and percentages, build the constraints.
if boundStyle == "UpperBound"
    conIneqValue =  volume - initialVolume*percentage;
    conEqValue = [];
    conIneqGrad = [];
    conEqGrad = [];
elseif boundStyle == "LowerBound"
    conIneqValue =  - volume + initialVolume*percentage;
    conEqValue = [];
    conIneqGrad = [];
    conEqGrad = [];
elseif boundStyle == "Equality"
    conIneqValue =  [];
    conEqValue = volume - initialVolume*percentage;
    conIneqGrad = [];
    conEqGrad = [];
elseif boundStyle =="LowerAndUpper"
    conIneqValue =  [volume - initialVolume*p1;...
                    - volume + initialVolume*p2];
    conEqValue = [];
    conIneqGrad = [];
    conEqGrad = [];
else
    error("Given Bound String Does not Match the Avalaible Cases...." + ...
        "Please Check the Defined Bound Strings in The Documentation!")
end


end
% End of code.
