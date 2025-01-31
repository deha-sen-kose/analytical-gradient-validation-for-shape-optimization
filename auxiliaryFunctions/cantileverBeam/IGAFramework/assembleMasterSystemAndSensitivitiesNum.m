function [masterStiffnessMtx, forceVector, stiffnessSensitivities,...
    stiffnessHessian, forceSensitivities, A] = ...
    assembleMasterSystemAndSensitivitiesNum(knotVectorXi, knotVectorEta,...
    polOrderXi, polOrderEta, controlPointsSurface, materialProp, ...
    distributedLoadVector, quadraticProgramming)
%%  FUNCTION assembleMasterStiffMtxandSensitivitiesNum
%
%   The function calls the element stiffness matrix, force vector and 
%   sensitivities function, namely computeSurfaceElementStiffnessMtxNum. 
%   It uses the element freedom tables to perform the assembly. It must be 
%   noted that the returned master stiffness matrix is the general 
%   stiffness matrix and it is not reduced.
%
%   Author(s)               : Deha Şen Köse, dehasenkose@gmail.com
%
%%  Reference(s):
%
%   [1] Kiendl, J., Schmidt, R., Wüchner, R., & Bletzinger, K. (2014). 
%   Isogeometric shape optimization of shells using semi-analytical 
%   sensitivity analysis and sensitivity weighting. Computer Methods in 
%   Applied Mechanics and Engineering, 274, 148-167. 
%   https://doi.org/10.1016/j.cma.2014.02.001
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
%       polOrderXi          : Polynomial order (int) in xi-direction.
%
%       polOrderEta         : Polynomial order (int) in eta-direction.
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
%   distributedLoadVector   : Distributed load in x- and y-directions. Has
%                             to be of size (1,2).
%
%   quadraticProgramming    : A boolean for calculation of the hessian. If
%                             true, the Hessian of the stiffness matrix is
%                             computed. It must be noted that these
%                             matrices are not used. BFGS
%                             approximation is used in the optimization. In
%                             addition, these hessians are not validated,
%                             so please use cautiously. 
%
%%  Output(s):
%
%   
%   masterStiffnessMtx      : The master stiffness matrix of the given
%                             plate (2*n_CP, 2*n_CP). It is not reduced
%
%   forceVector             : Force vector of the finite element system.
%                             Currently only the distributed forces are
%                             considered.
%
%   
%   stiffnessSensitivities  : The sensitivities (gradients) of the master 
%                             stiffness matrix with respect to the control
%                             point coordinates x and y. The sequence of
%                             the output is as follows: [dK/dx1, dK/dy1,
%                             dK/dx2, dK/dy2, ..., dK/dxN, dK/dyN].
%
%   stiffnessHessian        : Hessian of the master stiffness matrix. It is
%                             a cell matrix where each cell contains a
%                             double derivative of the master stiffness
%                             matrix. BFGS approximation of the hessians are
%                             used for the SQP algorithm. Thus exists the
%                             quadraticProgramming boolean as an input to
%                             control these computations.
%
%   forceSensitivities      : The sensitivities (gradients) of the global
%                             force vector with respect to the control
%                             point coordinates x and y. The sequence of
%                             the output is as follows: [dF/dx1, dF/dy1,
%                             dF/dx2, dF/dy2, ..., dF/dxN, dF/dyN].
%
%   A                       : The diagonal filtering matrix from the
%                             references above.
%

%% End of documentation - Code

% Get the number of elements in both xi and eta directions.
xiNumElements = length(findNonZeroKnotSpans(knotVectorXi));
etaNumElements = length(findNonZeroKnotSpans(knotVectorEta));

% Create the element indices for xi- and eta-directions.
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

% Concatenate the index vectors to obtain the final index matrix.
elementIndexVector = cat(2, elementIndexVectorXi, elementIndexVectorEta);

% Initialize the master stiffness matrix and the gradient matrices.
masterStiffnessMtx = zeros(height(controlPointsSurface)*2);
forceVector = zeros(height(controlPointsSurface)*2, 1);
stiffnessSensitivities = cell(height(controlPointsSurface)*2, 1);
stiffnessHessian = cell(height(controlPointsSurface)*2);
for ii=1:length(stiffnessSensitivities)

    stiffnessSensitivities{ii} = zeros(height(controlPointsSurface)*2);

end

forceSensitivities = cell(height(controlPointsSurface)*2, 1);
for ii=1:length(forceSensitivities)

    forceSensitivities{ii} = zeros(height(controlPointsSurface)*2,1);

end

for ii=1:width(stiffnessHessian)
    for jj=1:height(stiffnessHessian)

        stiffnessHessian{ii, jj} = zeros(height(controlPointsSurface)*2);

    end
end

% Initialize sensitivity weighting matrix [1]. 
A = zeros(height(controlPointsSurface));

% Loop over the elements.
for elemIndex=1:height(elementIndexVector)

    % Call the function for each element and compute the stiffness matrices
    % , sensitivities and element freedoom tables.
    [elStiffnessMtx, elStiffnessSensitivity, elConsistentForceVecVoigt, elForceSensitivity, elStiffnessHessian, a, elCPTable, elDOFTable] = ...
         computeSurfaceElStffnssMtxFrcVcNum(...
         controlPointsSurface, knotVectorXi, knotVectorEta, polOrderXi, ...
         polOrderEta, elementIndexVector(elemIndex,:), materialProp, distributedLoadVector, quadraticProgramming);

    % Assemble the stiffness matrices. 
    masterStiffnessMtx(elDOFTable, elDOFTable) = masterStiffnessMtx(elDOFTable, elDOFTable) + elStiffnessMtx;

    % Assemble the force vectors.
    forceVector(elDOFTable) = forceVector(elDOFTable) + elConsistentForceVecVoigt;

    % Assemble the sensitivity filtering [1].
    A(elCPTable, elCPTable) = A(elCPTable, elCPTable) + a;  

    % Loop over the design variables.
    for ii=1:1:length(stiffnessSensitivities)
        
        % Assemble the stiffness matrix sensitivities.
        stiffnessSensitivities{ii}(elDOFTable, elDOFTable) = stiffnessSensitivities{ii}(elDOFTable, elDOFTable) + elStiffnessSensitivity{ii};

        % Assemble the force vector sensitivities.
        forceSensitivities{ii}(elDOFTable,1) = forceSensitivities{ii}(elDOFTable,1) + elForceSensitivity{ii};


    end

    % Check if the hessians are wanted to be computed.
    if quadraticProgramming
        for ii=1:width(stiffnessHessian)
            for jj=1:length(stiffnessHessian)

                % Assemble the stiffness hessians.
                stiffnessHessian{ii, jj}(elDOFTable, elDOFTable) = stiffnessHessian{ii, jj}(elDOFTable, elDOFTable) + ...
                    elStiffnessHessian{ii, jj};


            end
        end
    end

end

end
%  End of code

