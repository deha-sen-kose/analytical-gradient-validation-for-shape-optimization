function [elStiffnessMtx, elStiffnessSensitivity, ...
         elConsistentForceVecVoigt, elForceSensitivityVoigt,...
         elStiffnessHessian, a, elCPTable, elDOFTable] = ...
         computeSurfaceElStffnssMtxFrcVcNum(controlPointsSurface,...
         knotVectorXi, knotVectorEta, polOrderXi, polOrderEta,...
         elementIndex, materialProp, distributedLoadVector,...
         quadraticProgramming)
%%  FUNCTION computeSurfaceElementStiffnessMtxNum
%
%   This function computes the element stiffness matrix, force vector of a 
%   plate in membrane action, as well as its element freedoom table that can 
%   be used in the assembly process and the sensitivities of it w.r.t. the 
%   control point coordinates. The thickness is assumed to be much smaller 
%   than the other two dimensions of the plate, thus the formulation is a 
%   plane-stress problem (a parallel force in either or both two dimensions 
%   are considered).
%
%   Author(s)               : Deha Şen Köse, deha.koese@tum.de
% 
%% Reference(s): 
%                 ISOGEOMETRIC STRUCTURAL ANALYSIS AND DESIGN, 
%                 Dr.-Ing. Roland Wüchner, Michael Breitenberger M.Sc.
%                 (hons), Anna Bauer M.Sc. Summer Term 2019.
%
%                 Structural sensitivity analysis and optimization 1.
%                 Kyung K. Choi, Nam Ho Kim. Springer-Verlag New York 2005.
%                 https://doi.org/10.1007/b138709 
%
%                 Kiendl, J., Schmidt, R., Wüchner, R., & Bletzinger, K. 
%                 (2014). Isogeometric shape optimization of shells using 
%                 semi-analytical sensitivity analysis and sensitivity 
%                 weighting. Computer Methods in Applied Mechanics and 
%                 Engineering, 274, 148-167. 
%                 https://doi.org/10.1016/j.cma.2014.02.001
% 
%% Input(s):
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
%       elementIndex        : The index of the current element. It is a 
%                             (1x2) array of integers. e.g.,
%                             *----------------------------------------*
%                             |                    |                   |
%                             |       (1,2)        |         (2,2)     |
%                             |                    |                   |
%                             |--------------------|-------------------|
%                             |                    |                   |
%                             |       (1,1)        |         (2,1)     |
%                             |                    |                   |
%                             *----------------------------------------*
%
%       materialProp        : Matlab structure that contains the following
%                             material properties. .E = Young's Modulus.
%                             .nu = Poisson's Ratio. .t = Thickness of the
%                             plate.
%
%     distributedLoadVector : Distributed load in x- and y-directions. Has
%                             to be of size (1,2).
%
%     quadraticProgramming  : A boolean for calculation of the hessian. If
%                             true the hessian will be computed. These
%                             hessians are not used in the thesis and are
%                             not validated!
%
%% Output(s):
%
%       elStiffnessMtx      : (2*nCP_el x 2*nCP_el) element stiffness
%                             matrix. nCP_el is the number of total active
%                             control points on that element. The active
%                             control points are found in this script.
%       
%   elStiffnessSensitivity  : Cell array of (2*nCPx1). Each cell contins the
%                             sensitivity of the element stiffness matrix
%                             (2*nCP_el x 2*nCP_el) with respect to the 
%                             design variables (control points x- and 
%                             y-coordinates). The sequence is pointwise and
%                             first x and second y comes.
%                             e.g., transpose({dKe/dx1, dKe/dy1, dKe/dx2, 
%                             dKe/dy2, ..., dKe/dxnCP, dKe/dynCP})
%
% elConsistentForceVecVoigt : Element consistent force vector for 
%                             distributed load in matrix-vector 
%                             multiplication notation.
%
%   elForceSensitivityVoigt : Sensitivities of the element force vector
%                             w.r.t. the control point coordinates in 
%                             matrix-vector multiplication notation.
%
%     elStiffnessHessian    : Hessian of the element stiffness matrix. It 
%                             is a cell matrix where each cell contains a
%                             double derivative of the element stiffness
%                             matrix.
%
%
%            a              : The filtering matrix components that are
%                             integrated elementwise.
%
%         elCPTable         : The element control point table. A vector
%                             containig the control point indices that are
%                             active in an element.
%
%         elDOFTable        : The element freedoom table to be used in the
%                             assembly process.
%
%% End of documentation - Code

% Get the number of control points in eta-direction (1D). Used later in
% creating the element control point table.
numCPY = length(knotVectorEta) - polOrderEta - 1;

% Finding the active control points for all the elements in the plate.
[xiActive, etaActive] = findActiveCPs(knotVectorXi, knotVectorEta,...
    polOrderXi, polOrderEta);

% Finding the non-zero knot spans (elements in IGA) in knot vectors.
% The indices of the element start and end is also found. 
elemBoundXi = findNonZeroKnotSpans(knotVectorXi);
elemLowerBoundXi = elemBoundXi(elementIndex(1));
elemUpperBoundXi = elemBoundXi(elementIndex(1)) + 1;

elemBoundEta = findNonZeroKnotSpans(knotVectorEta);
elemLowerBoundEta = elemBoundEta(elementIndex(2));
elemUpperBoundEta = elemBoundEta(elementIndex(2)) + 1;

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

% Compute the element freedoom table.
elDOFTable = zeros(length(elCPTable)*2,1);
for ii=1:length(elCPTable)

    elDOFTable(2*ii-1) = 2*elCPTable(ii)-1;
    elDOFTable(2*ii) = 2*elCPTable(ii);

end


% Getting the active control point coordinates using the element control
% point table.
elControlPoints = zeros(length(xiActive)*length(etaActive),2);
for ii=1:length(elCPTable)

    elControlPoints(ii,:) = controlPointsSurface(elCPTable(ii),:);

end


% Get the loads on the element
% elementLoadVector = loadVector(elDOFTable);


% Compute the jacobian of transformation between parametric space (xi,eta)
% and Gaussian space (xi,eta).
J_xi_u = ...
 [(knotVectorXi(elemUpperBoundXi)-knotVectorXi(elemLowerBoundXi))/2, 0;
 0, (knotVectorEta(elemUpperBoundEta)-knotVectorEta(elemLowerBoundEta))/2];
det_J_xi_u = det(J_xi_u);


% Read material properties and thickness of the membrane.
E = materialProp.E;
nu = materialProp.nu;
t = materialProp.t;

% Create the material matrix. 
D = (E/(1-nu^2)).*[1,  nu,  0;
                   nu, 1,   0;
                   0,  0,   (1-nu)/2];


% Initialize the stiffness matrix, B-matrix, Derivative of B-matrix and 
% sensitivity cell matrix.
elStiffnessMtx = zeros(length(elDOFTable));
elStiffnessSensitivity = cell(length(controlPointsSurface)*2,1);
for ii=1:length(elStiffnessSensitivity)
    elStiffnessSensitivity{ii} = zeros(length(elDOFTable));
end

BMtx = zeros(3,2*height(elControlPoints));
BMtx_sens = zeros(3,2*height(elControlPoints));

% Initialize B-matrices for second order partial derivatives
BMtx_sensSi = zeros(3,2*height(elControlPoints));
BMtx_sensSj = zeros(3,2*height(elControlPoints));

elConsistentForceVec = zeros(length(xiActive)*length(etaActive), 1);
elForceSensitivity = cell(length(controlPointsSurface)*2,1);
for ii=1:length(elForceSensitivity)
    elForceSensitivity{ii} = zeros(length(elDOFTable)/2, 2);
end



% Initialize the first derivative matrix.
%elementBasisFunctionVector = zeros(length(xiActive)*length(etaActive), 1);
elementBasisFunctionFirstDerivs = zeros(2, length(xiActive)*...
                                                length(etaActive));
elementBasisFunctionVector = zeros(1, length(xiActive)*length(etaActive));

% Create gradient matrices for the sensitivity calculation.
dX_ds = cell(height(controlPointsSurface)*2,1);
for ii=1:height(controlPointsSurface)*2
    dX_ds{ii} = zeros(length(elCPTable),2);
end

% Fill in the gradient matrices for the active control points in the
% element.
for ii=1:height(elCPTable)
    dX_ds{2*elCPTable(ii)}(ii, 2) = 1;
    dX_ds{2*elCPTable(ii)-1}(ii,1) = 1;
end

% Get gauss points. The polynomial order values considered in the 
% gaussWeightsandPoints function are the integrand polynomial orders. 
[weightsU, quadraturePointsU] = ...
    gaussWeightsandPoints(2*(polOrderXi+polOrderEta)-2);
[weightsV, quadraturePointsV] = ...
    gaussWeightsandPoints(2*(polOrderXi+polOrderEta)-2);
% [weightsU, quadraturePointsU] = ...
%     gaussWeightsandPoints(2*polOrderXi);
% [weightsV, quadraturePointsV] = ...
%     gaussWeightsandPoints(2*polOrderEta);

% Sensitivity filtering vector initialization.
ax = zeros(width(xiActive),1);

% Perform gauss integration for element stiffness and sensitivity matrices.
% Loop through the gauss quadrature(s).
for ii=1:length(weightsU)
    for jj=1:length(weightsV)

        % As the quadrature points are in the Gaussian space, and the
        % basis functions are computed in the parametric space, the
        % quadrature points are mapped from the Gaussian space to the
        % parametric space.
        xi_num =  quadraturePointsU(ii)*(knotVectorXi(elemUpperBoundXi)...
                  - knotVectorXi(elemLowerBoundXi))/2 +...
                  (knotVectorXi(elemUpperBoundXi) + ...
                  knotVectorXi(elemLowerBoundXi))/2;

        eta_num =quadraturePointsV(jj)*(knotVectorEta(elemUpperBoundEta)...
                  - knotVectorEta(elemLowerBoundEta))/2 ...
                  + (knotVectorEta(elemUpperBoundEta) + ...
                  knotVectorEta(elemLowerBoundEta))/2;

        % These line of functions can be used in order to check if the
        % gauss points are on the boundaries of the elements! Since the
        % derivatives do not exist in these values, integration cannot be
        % exact if the error is raised!
        if xi_num == knotVectorXi(elemLowerBoundXi) || xi_num == knotVectorXi(elemUpperBoundXi) ||...
                eta_num == knotVectorEta(elemLowerBoundEta) || eta_num == knotVectorEta(elemUpperBoundEta)
            fprintf('Gauss point found on the derivative with no seat! Possible that matrices are wrong! \n')
        end

        % Fill in the shape function derivative matrix.
        counter = 1;
        for kk=1:length(xiActive)

            [xiBasisFunc, xiBasisFuncDerivative] = computeBSplineBasisFunctionAndDerivatives(xiActive(kk), polOrderXi, xi_num, knotVectorXi);

            for zz=1:length(etaActive)
    
                [etaBasisFunc, etaBasisFuncDerivative] = computeBSplineBasisFunctionAndDerivatives(etaActive(zz), polOrderEta, eta_num, knotVectorEta);
                % the numbering goes as R_i,j = [R_1,1 , R_1,2, R_1,3, R_2,1 .....
                % R_m,n]. So first the x-directional (or xi-) CPs are filled.
                elementBasisFunctionVector(counter) = xiBasisFunc*etaBasisFunc;
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

        J_inv = inv(J_x_xi);

        % Compute the derivatives with respect to the CP coordinates.
        dR_dx_Mtx = J_inv*elementBasisFunctionFirstDerivs;

        % Rearrange the derivatives in B-matrix.
        for ll=1:width(dR_dx_Mtx)

            BMtx(1,2*ll-1) = dR_dx_Mtx(1,ll);
            BMtx(1,2*ll) = 0;

            BMtx(2,2*ll-1) = 0;
            BMtx(2,2*ll) = dR_dx_Mtx(2,ll);

            BMtx(3,2*ll-1) = dR_dx_Mtx(2,ll);
            BMtx(3,2*ll) = dR_dx_Mtx(1,ll);

        end

        % Perform gauss integration for the stiffness matrix.
        elStiffnessMtx = elStiffnessMtx + transpose(BMtx)*D*BMtx*...
            det(J_x_xi).*weightsU(ii).*weightsV(jj).*det_J_xi_u.*t;
        

        % Perform gauss integration for the sensitivities.
        for ss=1:length(elStiffnessSensitivity)
            % Calculate the derivatives of the elements in B matrix.
            dR_dx_Mtx_sens = -J_inv*elementBasisFunctionFirstDerivs*...
                dX_ds{ss}*J_inv*elementBasisFunctionFirstDerivs;

            % Rearrange the derivatives for B Matrix.
            for ll=1:width(dR_dx_Mtx_sens)

                BMtx_sens(1,2*ll-1) = dR_dx_Mtx_sens(1,ll);
                BMtx_sens(1,2*ll) = 0;

                BMtx_sens(2,2*ll-1) = 0;
                BMtx_sens(2,2*ll) = dR_dx_Mtx_sens(2,ll);

                BMtx_sens(3,2*ll-1) = dR_dx_Mtx_sens(2,ll);
                BMtx_sens(3,2*ll) = dR_dx_Mtx_sens(1,ll);

            end

            % Perform gauss integration for the sensitivities.
            elStiffnessSensitivity{ss} = elStiffnessSensitivity{ss} + ...
                (transpose(BMtx_sens)*D*BMtx.*det(J_x_xi) + transpose(BMtx)...
                *D*BMtx_sens.*det(J_x_xi) + transpose(BMtx)*D*BMtx....
                *det(J_x_xi).*trace(J_inv*elementBasisFunctionFirstDerivs...
                *dX_ds{ss}))*weightsU(ii).*weightsV(jj).*det_J_xi_u.*t;

        end

        % Initialize the element stiffness hessians.
        elStiffnessHessian = cell(height(controlPointsSurface)*2);
        for hh=1:height(elStiffnessHessian)
            for hz=1:width(elStiffnessHessian)

                elStiffnessHessian{hh,hz} = zeros(2*height(elControlPoints));

            end
        end

        % If these hessians are wished to be computed.
        if quadraticProgramming

            % Loop through the sensitivities
            for si=1:length(elStiffnessSensitivity)
                for sj=1:length(elStiffnessSensitivity)

                    % Calculate the derivatives of the elements in B matrix.
                    dR_dx_Mtx_sensSi = -J_inv*elementBasisFunctionFirstDerivs*...
                        dX_ds{si}*J_inv*elementBasisFunctionFirstDerivs;

                    % Calculate the derivatives of the elements in B matrix.
                    dR_dx_Mtx_sensSj = -J_inv*elementBasisFunctionFirstDerivs*...
                        dX_ds{sj}*J_inv*elementBasisFunctionFirstDerivs;


                    % Rearrange the derivatives for B Matrix.
                    for ll=1:width(dR_dx_Mtx_sensSi)
                        % For i-th parameter
                        BMtx_sensSi(1,2*ll-1) = dR_dx_Mtx_sensSi(1,ll);
                        BMtx_sensSi(1,2*ll) = 0;

                        BMtx_sensSi(2,2*ll-1) = 0;
                        BMtx_sensSi(2,2*ll) = dR_dx_Mtx_sensSi(2,ll);

                        BMtx_sensSi(3,2*ll-1) = dR_dx_Mtx_sensSi(2,ll);
                        BMtx_sensSi(3,2*ll) = dR_dx_Mtx_sensSi(1,ll);

                        % For j-th parameter
                        BMtx_sensSj(1,2*ll-1) = dR_dx_Mtx_sensSj(1,ll);
                        BMtx_sensSj(1,2*ll) = 0;

                        BMtx_sensSj(2,2*ll-1) = 0;
                        BMtx_sensSj(2,2*ll) = dR_dx_Mtx_sensSj(2,ll);

                        BMtx_sensSj(3,2*ll-1) = dR_dx_Mtx_sensSj(2,ll);
                        BMtx_sensSj(3,2*ll) = dR_dx_Mtx_sensSj(1,ll);

                    end
                    
                    % Derivatives of the jacobian w.r.t. i-th and j-th
                    % design variables
                    dJ_x_xi_dsi = det(J_x_xi).*trace(J_inv*elementBasisFunctionFirstDerivs*dX_ds{si});
                    dJ_x_xi_dsj = det(J_x_xi).*trace(J_inv*elementBasisFunctionFirstDerivs*dX_ds{sj});

                    % Compute the hessian.
                    elStiffnessHessian{si, sj} = elStiffnessHessian{si, sj} + (det_J_xi_u*weightsU(ii)*weightsV(jj)*t)*...
                        (transpose(BMtx_sensSi)*D*BMtx_sensSj*det(J_x_xi) + transpose(BMtx_sensSi)*D*BMtx*dJ_x_xi_dsj + ...
                        transpose(BMtx_sensSj)*D*BMtx_sensSi*det(J_x_xi) + transpose(BMtx)*D*BMtx_sensSi*dJ_x_xi_dsj + ...
                        transpose(BMtx_sensSj)*D*BMtx*dJ_x_xi_dsi + transpose(BMtx)*D*BMtx_sensSj*dJ_x_xi_dsi);
                    



                end
            end



        end

    end

    % Sensitivity Filtering
    for aa=1:length(xiActive)
        for bb=1:length(xiActive)
            ax(aa) = ax(aa) + computeBSplineBasisFunctionAndDerivatives...
                (xiActive(bb),polOrderXi,xi_num,knotVectorXi)*J_x_xi(1,1)*J_xi_u(1,1);
        end
    end

    % Numerical integration for element force vector
    elConsistentForceVec = elConsistentForceVec + ...
        transpose(elementBasisFunctionVector)*distributedLoadVector...
        *weightsU(ii)*J_x_xi(1,1)*J_xi_u(1,1)*t;


    % Compute the element force vector derivatives w.r.t. the design
    % variables.
    for bb=1:length(elForceSensitivity)

        J_x_xi_sens = elementBasisFunctionFirstDerivs*dX_ds{bb};
        elForceSensitivity{bb} = elForceSensitivity{bb} + ...
            transpose(elementBasisFunctionVector)*distributedLoadVector...
            *weightsU(ii)*J_x_xi_sens(1,1)*J_xi_u(1,1)*t;

    end

end
n = length(xiActive);
ax = repmat(ax',n);
ax = ax(1,:);
a = diag(ax);

% Initialize the force sensitivities in matrix-vector multiplication
% format.
elForceSensitivityVoigt = cell(height(controlPointsSurface)*2,1);
for ii=1:length(elForceSensitivityVoigt)

    elForceSensitivityVoigt{ii} = zeros(length(elDOFTable),1);

end

% Assign the force vector values in matrix-vector multiplication format.
elConsistentForceVecVoigt = zeros(height(elConsistentForceVec)*2,1);
for ii=1:length(elConsistentForceVec)

    elConsistentForceVecVoigt(2*ii-1) = elConsistentForceVec(ii,1);
    elConsistentForceVecVoigt(2*ii) = elConsistentForceVec(ii,2);


end

% Assign the force sensitivity values in matrix-vector multiplication format.
for ii=1:length(elForceSensitivity)

    for jj=1:length(elForceSensitivity{ii})

        elForceSensitivityVoigt{ii}(2*jj-1) = elForceSensitivity{ii}(jj,1);
        elForceSensitivityVoigt{ii}(2*jj) = elForceSensitivity{ii}(jj,2);

    end

end

end
% End of code.