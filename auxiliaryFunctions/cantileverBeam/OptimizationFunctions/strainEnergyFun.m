function [strainEnergy, gradient, hessian, A] = strainEnergyFun(...
          knotVectorXi, knotVectorEta, polOrderXi, polOrderEta, ...
          controlPointsSurfaceOrig, materialProp, distributedLoadVector, ...
          homDOFs, designVariablesDOFs, quadraticProgramming, ...
          designVariables, forceType, pointLoadVector, symtryC, symtrynum)
%% FUNCTION strainEnergyFun
%
%   This function runs the codes for element stiffness matrices and force 
%   vectors as well as the sensitivities, and returns the strain energy 
%   together with the gradient of the strain energy since it is easier to
%   to use optimizers in Matlab.
%
%   Author(s)           : Deha Şen Köse, dehasenkose@gmail.com
%
%% Reference(s):
%
%   Kiendl, J., Schmidt, R., Wuchner, R., & Bletzinger, K.-U. (2014). Isogeometric shape
%	optimization of shells using semi-analytical sensitivity analysis and sensitivity
%	weighting. Computer Methods in Applied Mechanics and Engineering, 274, 148–
%	167. https://doi.org/https://doi.org/10.1016/j.cma.2014.02.001
%
%% Input(s):
%
%       knotVectorXi        : The knot vector in xi-direction defined in
%                             the parametric space. m=n+p+1 and the
%                             non-descending character must be satisfied.
%
%       knotVectorEta       : The knot vector in eta-direction defined in
%                             the parametric space. m=n+p+1 and the 
%                             non-descending character must be satisfied.
%
%       polOrderXi          : Polynomial order (int) in xi-direction.
%
%       polOrderEta         : Polynomial order (int) in Eta-direction.
%
%   controlPointsSurfaceOrig: (nx2) array of the coordinates of the control
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
%     distributedLoadVector : Distributed load in x- and y-directions. Has
%                             to be of size (1,2).
%
%       homDOFs             : Numbering of the homogeneous degrees of 
%                             freedoom. In this framework, the degrees of 
%                             freedoom are in x- and y-displacement pairs. 
%                             Also, the global numbering of DOFs or rather 
%                             control points are first incremented in y- 
%                             and then x-direction. In this code, always 
%                             the first elements are considered to be on
%                             the cantilever edge. e.g., here the
%                             constrainted CPs should be 1,2,3,4,5.
%                     
%                             *5 ---------*10 --------- *15 --------- *20
%                             *4 --------- *9 --------- *14 --------- *19
%                             *3 --------- *8 --------- *13 --------- *18
%                             *2 --------- *7 --------- *12 --------- *17
%                             *1 --------- *6 --------- *11 --------- *16
%
%   designVariablesDOFs     : Numbering of degrees of freedoom of the
%                             design variables.
%
%   quadraticProgramming    : A boolean for calculation of the hessian.
%
%   designVariables         : The design variable values. They are
%                             incorporated into the controlPointsSurface in
%                             the optimization iterations through 
%                             designVariablesDOFs. 
%
%   forceType               : The type of the force acting on the plate.
%                             Currently there is only one option, that is
%                             "RightEdgePointLoad". If this option is not
%                             given as a string, the program uses the given
%                             distributed force vector and calculates the
%                             forces accordingly. For the force calculation
%                             please check the corresponding code files.
%
%   pointLoadVector         : The point load vector (1xn). Where n is the
%                             number of control points on the right edge
%                             times 2. Can be left [] if the load is
%                             distributed load.
%
%   symtryC                 : A boolean for the symmetrical control point
%                             linking. If true, a symmetry is assumed
%                             w.r.t. the mid x-line using the symtrynum
%                             input.
%
%   symtrynum               : A vector containing the -X control point 
%                             (meaning that they are below the assumed 
%                             symmetry line for the design variables) 
%                             numbers that are assumed to be optimized 
%                             symmetrical to the corresponding above
%                             control points. Can be left [] if symtryC is
%                             false.
%
%% Output(s):
%
%   strainEnergy            : Scalar strain energy of the surface under 
%                             considered loading.
%
%   gradient                : Is the gradient of the strain energy with
%                             respect to the design variables. They will be
%                             in the same order as the given desing 
%                             variable degree of freedoom vector. 
%
%   hessian                 : Hessian of the strain energy function w.r.t.
%                             the design variable coordinates. Is in the
%                             same order with the given designVariablesDOFs
%                             which should be ascending. They are not 
%                             validated. BFGS approximation is used in 
%                             optimization.
%
%   A                       : The diagonal filtering matrix.
%
%% End of function definition - Code

% Get the number of design variables in the optimization process.
numDesignVar = length(designVariables);

% Get the number of control points in y-direction.
numCPy = length(knotVectorEta) - polOrderEta - 1;
numCP = height(controlPointsSurfaceOrig);
controlPointsSurface = controlPointsSurfaceOrig;
% Loop through the design variables to update the surface control points 
% matrix respecting the control polygons. 
for ii=1:numDesignVar

    % Get the DOF of the current desing variable.
    dof = designVariablesDOFs(ii);
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

% If the symmetric linking is assumed, then the control points are oriented
% accordingly.
if symtryC == true
    controlPointsSurface = symmetryFun(controlPointsSurfaceOrig, ...
        controlPointsSurface, numCPy, symtrynum);
end


% Compute the master stiffness matrix, force vector and sensitivities of
% those those with respect to ALL control point coordinates
[masterStiffnessMtx, forceVector, stiffnessSensitivities, ...
    stiffnessHessian, forceSensitivities, A] = ...
    assembleMasterSystemAndSensitivitiesNum(...
    knotVectorXi, knotVectorEta, polOrderXi, polOrderEta, ...
    controlPointsSurface, materialProp, distributedLoadVector, ...
    quadraticProgramming);

% Apply homogeneous dirichlet boundary conditions and reduce the system
[reducedStiffnessMtx, reducedForceVector] =applyHomDirichletBCs(homDOFs,...
    masterStiffnessMtx, forceVector);

% Change the force vector if point load is assumed on the rhs.
if forceType == "RightEdgePointLoad"
    if length(pointLoadVector) ~= 2*numCPy
        error("Wrong Number of forces in the force vector!")
    else
        reducedForceVector = [zeros(1,(2*numCP-2*numCPy-length(homDOFs))),...
            pointLoadVector]';
    end
    forceSensitivities = cell(2*numCP,1);
    for ii=1:length(forceSensitivities)
        forceSensitivities{ii} = zeros(2*numCP,1);
    end
elseif forceType ~= "DistributedLoad"
    error("Given loading condition is not supported!")
end

% Solve the IGA Ku = F linear system of equations
[displacementField] = solveIGASystem(reducedStiffnessMtx,...
    reducedForceVector);

% Reconstruct the displacement field by adding the homogeneous degrees of
% freedoom.
u = zeros(length(masterStiffnessMtx),1);
counter = 1;
for ii=1:length(u)
    k = find(homDOFs==ii, 1);
    if isempty(k)
        u(ii) = displacementField(counter);
        counter = counter+1;
    end
end
displacementField = u;

% Return the gradient vector of the strain energy with respect to the
% control point coordinates.
[gradient, hessian] = computeSensitivities(displacementField, ... 
    designVariablesDOFs, stiffnessSensitivities, forceSensitivities, ...
    stiffnessHessian);

% Compute the strain energy.
strainEnergy = 0.5.*transpose(displacementField)*masterStiffnessMtx*...
    displacementField;

% Making strain energy go to infinite when the shape malforms
if strainEnergy < 0
    strainEnergy = Inf;
end


end
% End of code