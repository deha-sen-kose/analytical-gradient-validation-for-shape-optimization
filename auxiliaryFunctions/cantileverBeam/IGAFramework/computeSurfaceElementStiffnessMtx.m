function [elStiffnessMtx, elStiffnessSensitivity, elDOFTable] = computeSurfaceElementStiffnessMtx(...
         controlPointsSurface, knotVectorXi, knotVectorEta, polOrderXi, ...
         polOrderEta, elementIndex, designVariables, designVariablesNum, materialProp, symbolicFlag)
%
%   Outdated SYMBOLIC function. It is not as efficient as the numerical
%   computations. Please use computeSurfaceElStffnssMtxFrcVctNum.m instead.
%

numCPY = length(knotVectorEta) - polOrderEta - 1;

[xiActive, etaActive] = findActiveCPs(knotVectorXi, knotVectorEta, polOrderXi, polOrderEta);

elemBoundXi = findNonZeroKnotSpans(knotVectorXi);
elemLowerBoundXi = elemBoundXi(elementIndex(1));
elemUpperBoundXi = elemBoundXi(elementIndex(1)) + 1;

elemBoundEta = findNonZeroKnotSpans(knotVectorEta);
elemLowerBoundEta = elemBoundEta(elementIndex(2));
elemUpperBoundEta = elemBoundEta(elementIndex(2)) + 1;

xiActive = xiActive(elementIndex(1),:);
etaActive = etaActive(elementIndex(2),:);


if symbolicFlag

    if any(isnan(symvar(designVariables)))

        error("The operations are considered symbolic, however the given desing variables are numeric!")

    end


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

    elDOFTable = zeros(length(elCPTable)*2,1);
    for ii=1:length(elCPTable)

        elDOFTable(2*ii-1) = 2*elCPTable(ii)-1;
        elDOFTable(2*ii) = 2*elCPTable(ii);

    end

    id = [];
    for ii=1:length(elCPTable)

        strx = append("x", string(elCPTable(ii)));
        stry = append("y", string(elCPTable(ii)));

        if any(ismember(designVariables, strx))

            id = [id, find(ismember(designVariables, strx))];

        end
        if any(ismember(designVariables, stry))

            id = [id, find(ismember(designVariables, stry))];

        end


    end

    elDesignVariables = designVariables(sort(id));

    elControlPoints = sym(zeros(length(xiActive)*length(etaActive),2));
    counterX = 1;
    counterY = 1;
    for ii=1:length(elCPTable)

        strx = append("x", string(elCPTable(ii)));
        stry = append("y", string(elCPTable(ii)));

        if ~any(ismember(elDesignVariables, strx))

            elControlPoints(counterX,1) = controlPointsSurface(elCPTable(ii),1);
            counterX = counterX + 1;
        else
            idx = ismember(elDesignVariables, strx);
            elControlPoints(counterX,1) = elDesignVariables(idx);
            counterX = counterX + 1;
        end

        if ~any(ismember(elDesignVariables, stry))

            elControlPoints(counterY,2) = controlPointsSurface(elCPTable(ii),2);
            counterY = counterY + 1;
        else
            idy = ismember(elDesignVariables, stry);
            elControlPoints(counterY,2) = elDesignVariables(idy);
            counterY = counterY + 1;
        end



    end


    designVariablesFlat = zeros(length(designVariablesNum)*2,1);%%%%%
    for ii=1:length(designVariablesNum)

        designVariablesFlat(2*ii-1) = designVariablesNum(ii,1);
        designVariablesFlat(2*ii)   = designVariablesNum(ii,2);

    end
    
    designVariablesFlat = transpose(designVariablesFlat);

    elDesignVariablesFlat = designVariablesFlat(id);


else


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

    elDOFTable = zeros(length(elCPTable)*2,1);
    for ii=1:length(elCPTable)

        elDOFTable(2*ii-1) = 2*elCPTable(ii)-1;
        elDOFTable(2*ii) = 2*elCPTable(ii);

    end

    elControlPoints = zeros(length(xiActive)*length(etaActive),2);

    for ii=1:length(elCPTable)

        elControlPoints(ii,:) = controlPointsSurface(elCPTable(ii),:);

    end


end



% compute the jacobian of transformation between real coordinates (x,y) and
% parametric space (xi,eta)
J_xi_u = [(knotVectorXi(elemUpperBoundXi)-knotVectorXi(elemLowerBoundXi))/2, 0;
          0, (knotVectorEta(elemUpperBoundEta)-knotVectorEta(elemLowerBoundEta))/2];
det_J_xi_u = det(J_xi_u);


% read material properties and thickness of the membrane
E = materialProp.E;
nu = materialProp.nu;
t = materialProp.t;

D = (E/(1-nu^2)).*[1, nu , 0;
                 nu, 1,  0;
                 0,  0, (1-nu)/2];


% Get gauss points. The polynomial order values are considered to be the integrand values! 
[weightsU, quadraturePointsU] = gaussWeightsandPoints(2*polOrderXi+2);
[weightsV, quadraturePointsV] = gaussWeightsandPoints(2*polOrderEta+2);


% Initialize the stiffness matrix, B matrix, N vector and matrix.
elStiffnessMtx = zeros(length(elDOFTable));
if ~symbolicFlag
    elStiffnessSensitivity = cell(length(controlPointsSurface)*2,1);
    for ii=1:length(elStiffnessSensitivity)

        elStiffnessSensitivity{ii} = zeros(length(elDOFTable));

    end
else
    elStiffnessSensitivity = cell(length(elDesignVariablesFlat),1);

    for ii=1:length(elStiffnessSensitivity)

        elStiffnessSensitivity{ii} = zeros(length(elDOFTable));

    end

end

if symbolicFlag
    BMtx = sym(zeros(3,2*height(elControlPoints)));
else
    BMtx = zeros(3,2*height(elControlPoints));
    BMtx_sens = zeros(3,2*height(elControlPoints));
end

%elementBasisFunctionVector = zeros(length(xiActive)*length(etaActive), 1);
elementBasisFunctionFirstDerivs = zeros(2, length(xiActive)*length(etaActive));


% create gradient matrices for the sensitivity calculation
if ~symbolicFlag
    dX_ds = cell(height(controlPointsSurface)*2,1);
    for ii=1:height(controlPointsSurface)*2
        dX_ds{ii} = zeros(length(elCPTable),2);
    end

    for ii=1:height(elCPTable)

        dX_ds{2*elCPTable(ii)}(ii, 2) = 1;
        dX_ds{2*elCPTable(ii)-1}(ii,1) = 1;
    end
end


for ii=1:length(weightsU)
    for jj=1:length(weightsV)

        xi_num =  quadraturePointsU(ii)*(knotVectorXi(elemUpperBoundXi) - knotVectorXi(elemLowerBoundXi))/2 ...
                  + (knotVectorXi(elemUpperBoundXi) + knotVectorXi(elemLowerBoundXi))/2;

        eta_num = quadraturePointsV(jj)*(knotVectorEta(elemUpperBoundEta) - knotVectorEta(elemLowerBoundEta))/2 ...
                  + (knotVectorEta(elemUpperBoundEta) + knotVectorEta(elemLowerBoundEta))/2;

        % These line of functions can be used in order to check if the
        % gauss points are on the boundaries of the elements!
        if xi_num == knotVectorXi(elemLowerBoundXi) || xi_num == knotVectorXi(elemUpperBoundXi) ||...
                eta_num == knotVectorEta(elemLowerBoundEta) || eta_num == knotVectorEta(elemUpperBoundEta)
            fprintf('Gauss point found on the derivative with no seat! Possible that matrices are wrong! \n')
        end

        counter = 1;
        for kk=1:length(xiActive)
            for zz=1:length(etaActive)

                [xiBasisFunc, xiBasisFuncDerivative] = computeBSplineBasisFunctionAndDerivatives(xiActive(kk), polOrderXi, xi_num, knotVectorXi);
                [etaBasisFunc, etaBasisFuncDerivative] = computeBSplineBasisFunctionAndDerivatives(etaActive(zz), polOrderEta, eta_num, knotVectorEta);
                % the numbering goes as R_i,j = [R_1,1 , R_1,2, R_1,3, R_2,1 .....
                % R_m,n]. So first the x-directional (or xi-) CPs are filled.
                elementBasisFunctionVector(counter) = xiBasisFunc*etaBasisFunc;
                elementBasisFunctionFirstDerivs(1, counter) = xiBasisFuncDerivative*etaBasisFunc;
                elementBasisFunctionFirstDerivs(2, counter) = etaBasisFuncDerivative*xiBasisFunc;
                counter = counter +1;

            end
        end


        J_x_xi = elementBasisFunctionFirstDerivs*elControlPoints;


        J_inv = inv(J_x_xi);

        dR_dx_Mtx = J_inv*elementBasisFunctionFirstDerivs;



        % sensitivities
        % for ss=1:height(controlPointsSurface)*2
        % 
        %     dX_ds = zeros(2, length(elCPTable));
        %     dX_ds(ss) = 
        % 
        %     dR_dx_Mtx_sens = -J_inv*elementBasisFunctionFirstDerivs*[]
        % end

        for ll=1:width(dR_dx_Mtx)

            BMtx(1,2*ll-1) = dR_dx_Mtx(1,ll);
            BMtx(1,2*ll) = 0;

            BMtx(2,2*ll-1) = 0;
            BMtx(2,2*ll) = dR_dx_Mtx(2,ll);

            BMtx(3,2*ll-1) = dR_dx_Mtx(2,ll);
            BMtx(3,2*ll) = dR_dx_Mtx(1,ll);

        end

        % if symbolicFlag
        %     BMtx_num = double(subs(BMtx, [xi, eta, reshape(activeControlPointsSurfaceSym,[1, height(activeControlPointsSurfaceSym)*2])],...
        %         [xi_num, eta_num, reshape(activeControlPointsSurface, [1, height(activeControlPointsSurface)*2])]));
        %     nana = isnan(BMtx_num);
        % 
        %     if any(nana==1)
        %         error("Some of the B-matrix values are NaN! Terminating!")
        %     end
        % else
        % 
        % 
        %     nana = isnan(BMtx);
        %     if any(nana==1)
        %         error("Some of the B-matrix values are NaN! Terminating!")
        %     end
        % 
        % 
        % end



        %numberOfZeros = nnz(~double(subs(BMtx,[xi, eta], [xi_num, eta_num])))
        %numberOfZeros = nnz(~BMtx)

        elStiffnessMtx = elStiffnessMtx + transpose(BMtx)*D*BMtx*det(J_x_xi).*weightsU(ii).*weightsV(jj).*det_J_xi_u.*t;

        if ~symbolicFlag

            elementBasisFunctionMtx = zeros(2, 2*length(xiActive)*length(etaActive));

            for zz=1:width(elementBasisFunctionVector)

                elementBasisFunctionMtx(1,2*zz-1) = elementBasisFunctionVector(zz);
                elementBasisFunctionMtx(1,2*zz) = 0;

                elementBasisFunctionMtx(2,2*zz-1) = 0;
                elementBasisFunctionMtx(2,2*zz) = elementBasisFunctionVector(zz);

            end



            for ss=1:length(elStiffnessSensitivity)
                dR_dx_Mtx_sens = -J_inv*elementBasisFunctionFirstDerivs*dX_ds{ss}*J_inv*elementBasisFunctionFirstDerivs;


                for ll=1:width(dR_dx_Mtx_sens)

                    BMtx_sens(1,2*ll-1) = dR_dx_Mtx_sens(1,ll);
                    BMtx_sens(1,2*ll) = 0;

                    BMtx_sens(2,2*ll-1) = 0;
                    BMtx_sens(2,2*ll) = dR_dx_Mtx_sens(2,ll);

                    BMtx_sens(3,2*ll-1) = dR_dx_Mtx_sens(2,ll);
                    BMtx_sens(3,2*ll) = dR_dx_Mtx_sens(1,ll);

                end

                elStiffnessSensitivity{ss} = elStiffnessSensitivity{ss} + ...
                    (transpose(BMtx_sens)*D*BMtx.*det(J_x_xi) + transpose(BMtx)*D*BMtx_sens.*det(J_x_xi) + transpose(BMtx)*D*BMtx.*det(J_x_xi).*trace(J_inv*elementBasisFunctionFirstDerivs*dX_ds{ss}))...
                    *weightsU(ii).*weightsV(jj).*det_J_xi_u.*t;

            end
        end


    end
end

if symbolicFlag

    for ss=1:length(elDesignVariables)


        snew = double(subs(diff(elStiffnessMtx, elDesignVariables(ss)), elDesignVariables, elDesignVariablesFlat));
        elStiffnessSensitivity{ss} = elStiffnessSensitivity{ss} + snew;


    end

end




if ~symbolicFlag
    nonsym = [];
    for ii=1:height(elStiffnessMtx)
        for jj=1:width(elStiffnessMtx)

            if abs(elStiffnessMtx(ii,jj) - elStiffnessMtx(jj,ii)) > 1e-6
                nonsym = [nonsym; ii, jj];

            end
        end
    end


    try chol(elStiffnessMtx)
        disp('Matrix is symmetric positive definite.')
    catch ME
        disp('Matrix is not symmetric positive definite')
    end
end



end

