function surfaceVolume = computeSurfaceVolume(knotVectorXi, knotVectorEta,...
    polOrderXi, polOrderEta, controlPointsSurface, propMaterial)
%%  FUNCTION computeSurfaceVolume
%
%   This function computes the surface volume of the plate.
%
%   Author(s)               : Deha Şen Köse, dehasenkose@gmail.com
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
%       propMaterial        : Matlab structure that contains the following
%                             material properties. .E = Young's Modulus.
%                             .nu = Poisson's Ratio. .t = Thickness of the
%                             plate.
%
%% Output(s):
%
%   surfaceVolume           : Volume of the plate, found by multplying with
%                             thickess.
%
%% End of function definition - Code

% Get the number of elements in both xi and eta directions.
xiNumElements = length(findNonZeroKnotSpans(knotVectorXi));
etaNumElements = length(findNonZeroKnotSpans(knotVectorEta));

% Create the element indices
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

% Concatenate the index vectors for xi and eta and obtain the overall
% vector.
elementIndexVector = cat(2, elementIndexVectorXi, elementIndexVectorEta);


% Get gauss points. The polynomial order values considered in the
% gaussWeightsandPoints function are the integrand polynomial orders.
[weightsU, quadraturePointsU] = gaussWeightsandPoints(2*(polOrderXi+polOrderEta));
[weightsV, quadraturePointsV] = gaussWeightsandPoints(2*(polOrderXi+polOrderEta));


% Initialize the first derivative matrix.
%elementBasisFunctionVector = zeros(length(xiActive)*length(etaActive), 1);
elementBasisFunctionFirstDerivs = zeros(2, (polOrderXi+1)*(polOrderEta+1));
%elementBasisFunctionVector = zeros(1, (polOrderXi+1)*(polOrderEta+1));

% Initialize constraint gradient vector.
% constraintGradient = zeros(height(controlPointsSurface)*2,1);
% elConstraintGradient = zeros(height(controlPointsSurface)*2,1);
% constraintHessian = zeros(height(controlPointsSurface)*2);
% elConstraintHessian = zeros(height(controlPointsSurface)*2);

% Initialize the surface area variable
totalArea = 0;
% Loop through all the elements. The integration takes place on the
% parametric space, thus no element overlapping is present.
for elemIndex=1:height(elementIndexVector)
    
    % Initialize the element area.
    elementArea = 0;
    % Get the element index.
    elementIndex = elementIndexVector(elemIndex,:);
    % Get the number of control points in eta-direction (1D). Used later in
    % creating the element control point table.
    numCPY = length(knotVectorEta) - polOrderEta - 1;

    % Finding the active control points for all the elements in the plate.
    [xiActive, etaActive] = findActiveCPs(knotVectorXi, knotVectorEta, polOrderXi, polOrderEta);

    % Finding the non-zero knot spans (elements in IGA) in knot vectors.
    % The indices of the element start and end is also found.
    elemBoundXi = findNonZeroKnotSpans(knotVectorXi);
    elemLowerBoundXi = elemBoundXi(elementIndex(1));
    elemUpperBoundXi = elemBoundXi(elementIndex(1)) + 1;

    elemBoundEta = findNonZeroKnotSpans(knotVectorEta);
    elemLowerBoundEta = elemBoundEta(elementIndex(2));
    elemUpperBoundEta = elemBoundEta(elementIndex(2)) + 1;

    % Compute the Jacobian between the parametric space and the Gaussian
    % space.
    J_xi_u = [(knotVectorXi(elemUpperBoundXi)-knotVectorXi(elemLowerBoundXi))/2, 0;
        0, (knotVectorEta(elemUpperBoundEta)-knotVectorEta(elemLowerBoundEta))/2];

    % Active control points for the specific element is indexed.
    xiActive = xiActive(elementIndex(1),:);
    etaActive = etaActive(elementIndex(2),:);

    % Element Control Point table is computed.
    elCPTable = zeros((polOrderEta+1)*(polOrderXi+1), 1);
    counter = 1;
    for ii=1:length(xiActive)
        for jj=1:length(etaActive)
            % get the CP numbers. Note: the CP numbering goes as: x1,y1 -
            % x1,y2 - x1,y3 - x1,y4 ... x2,y1 ... xn,ym
            elCPTable(counter) = (xiActive(ii)-1)*numCPY + etaActive(jj);
            counter = counter + 1;
        end
    end


    % Getting the active control point coordinates using the element control
    % point table.
    elControlPoints = zeros(length(xiActive)*length(etaActive),2);
    for ii=1:length(elCPTable)

        elControlPoints(ii,:) = controlPointsSurface(elCPTable(ii),:);

    end


    % Loop through the 2D Gauss points.
    for jj=1:length(weightsU)
        for mm=1:length(weightsV)

            % As the quadrature points are in the Gaussian space, and the
            % basis functions are computed in the parametric space, the
            % quadrature points are mapped from the Gaussian space to the
            % parametric space.
            xi_num =  quadraturePointsU(jj)*(knotVectorXi(elemUpperBoundXi) - knotVectorXi(elemLowerBoundXi))/2 ...
                + (knotVectorXi(elemUpperBoundXi) + knotVectorXi(elemLowerBoundXi))/2;

            eta_num = quadraturePointsV(mm)*(knotVectorEta(elemUpperBoundEta) - knotVectorEta(elemLowerBoundEta))/2 ...
                + (knotVectorEta(elemUpperBoundEta) + knotVectorEta(elemLowerBoundEta))/2;


            % Fill in the shape function derivative matrix.
            counter = 1;
            for kk=1:length(xiActive)

                [xiBasisFunc, xiBasisFuncDerivative] = computeBSplineBasisFunctionAndDerivatives(xiActive(kk), polOrderXi, xi_num, knotVectorXi);

                for zz=1:length(etaActive)

                    [etaBasisFunc, etaBasisFuncDerivative] = computeBSplineBasisFunctionAndDerivatives(etaActive(zz), polOrderEta, eta_num, knotVectorEta);
                    % the numbering goes as R_i,j = [R_1,1 , R_1,2, R_1,3, R_2,1 .....
                    % R_m,n]. So first the x-directional (or xi-) CPs are filled.
                    %elementBasisFunctionVector(counter) = xiBasisFunc*etaBasisFunc;
                    elementBasisFunctionFirstDerivs(1, counter) = xiBasisFuncDerivative*etaBasisFunc;
                    elementBasisFunctionFirstDerivs(2, counter) = etaBasisFuncDerivative*xiBasisFunc;
                    counter = counter +1;

                end
            end

            % Compute the jacobian of the element and its inverse.
            J_x_xi = elementBasisFunctionFirstDerivs*elControlPoints;

            % Correction for making the element stiffness matrices SPD. This
            % way the solver works faster.
            if abs(J_x_xi(2,1)) <= 1e-6
                J_x_xi(2,1) = 0;
            end
            if abs(J_x_xi(1,2)) <= 1e-6
                J_x_xi(1,2) = 0;
            end

            % Compute the numerical integration in parametric domain.
            elementArea = elementArea + norm(cross([J_x_xi(1,1),0,0],[0,J_x_xi(2,2),0]))*det(J_xi_u)*weightsU(jj)*weightsV(mm);


            % Compute the derivative of the Jacobian to calculate the
            % constraint gradient.
            % dX_ds = cell(height(controlPointsSurface)*2,1);
            % for ii=1:height(controlPointsSurface)*2
            %     dX_ds{ii} = zeros(length(elCPTable),2);
            % end
            % 
            % % Fill in the gradient matrices for the active control points in the
            % % element.
            % for ii=1:height(elCPTable)
            %     dX_ds{2*elCPTable(ii)}(ii, 2) = 1;
            %     dX_ds{2*elCPTable(ii)-1}(ii,1) = 1;
            % end
            % 
            % for kk=1:length(constraintGradient)
            % 
            % 
            %     dJ_x_xi_dsk = elementBasisFunctionFirstDerivs*dX_ds{kk};
            % 
            %     elConstraintGradient(kk) = (dJ_x_xi_dsk(1,1)*J_x_xi(2,2) + ...
            %         J_x_xi(1,1)*dJ_x_xi_dsk(2,2))*det(J_xi_u)*weightsU(jj)*weightsV(mm);
            % 
            %     for zz=1:length(constraintGradient)
            % 
            %         dJ_x_xi_dsz = elementBasisFunctionFirstDerivs*dX_ds{zz};
            % 
            %         elConstraintHessian(kk, zz) = elConstraintHessian(kk, zz) + (dJ_x_xi_dsk(1,1)*dJ_x_xi_dsz(2,2) + ...
            %         dJ_x_xi_dsz(1,1)*dJ_x_xi_dsk(2,2))*det(J_xi_u)*weightsU(jj)*weightsV(mm); 
            % 
            %     end
            % 
            % 
            % end


        end
    end



    % Add the element area into the total area.
    totalArea = totalArea + elementArea;

    % Add all the contributions from element derivatives
    % constraintGradient = constraintGradient + elConstraintGradient;
    % constraintHessian = constraintHessian + elConstraintHessian;
end

% compute the surface volume.
surfaceVolume = totalArea*propMaterial.t;

% Finalize the gradient value
% constraintGradient = constraintGradient.*propMaterial.t;
% constraintHessian = constraintHessian.*propMaterial.t;

end
%% End of code