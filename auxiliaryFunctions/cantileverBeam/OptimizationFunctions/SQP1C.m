function opt = SQP1C(objFunctionHandle, conFuncHandle, designVariables, lagrangeMultp, hessianCalc)
%SQP1C Summary of this function goes here
%   Detailed explanation goes here

%% Think about the lagrange multipliers!!!!!! lagrangeMultp

% Read the number of design variables
numDesignVariables = length(designVariables);

% check for hessian input, if not specified initialize it with identity
% matrix
if hessianCalc == "BFGS"

    lagrangeHessian = eye(numDesignVariables);

elseif hessianCalc == "Analytical"

    [~ , ~, objHessian] =  objFunctionHandle(designVariables);
    [~, ~, conHessian] = conFuncHandle(designVariables);
    lagrangeHessian = objHessian + lagrangeMultp*conHessian;

else
    error("Given method of hessian calculation as a string does not " + ...
        "correspond to any impelented method! " + ...
        "Available methods are: BFGS and Analytical")
end




% Initialize the incumbent step vector
step = ones(numDesignVariables,1);
iteration = 1;
while norm(step) >= 1e-6 && iteration < 100

    if hessianCalc == "Analytical"

        % Calculate the function and constraint values and gradients
        [objValue, objFuncGradient, objHessian] =  objFunctionHandle(designVariables);
        [constraintValues, conFuncGradient, conHessian] = conFuncHandle(designVariables);
        lagrangeHessian = objHessian + lagrangeMultp*conHessian;
        

    else 

        % Calculate the function and constraint values and gradients
        [objValue, objFuncGradient, ~] =  objFunctionHandle(designVariables);
        [constraintValues, conFuncGradient, ~] = conFuncHandle(designVariables);
    end



    violatedRealCon = find(constraintValues < 0);
    
    if isempty(violatedRealCon)

        penaltyFunc = objValue;

    else
        
        penaltyFunc = objValue + transpose(lagrangeMultp)*abs(constraintValues(violatedRealCon));

        
    end

    % Read the number of constraints
    numConstraints = length(constraintValues);

    % Initialize the LHS-matrix for the linear system of equations
    lhsEqMtx = zeros(numDesignVariables + numConstraints);

    %% DONT FORGET TO CHANGE THE HESSIAN!!!

    % Fill in the data of the LHS-matrix
    lhsEqMtx(1:numDesignVariables, 1:numDesignVariables) = lagrangeHessian;
    lhsEqMtx(1:numDesignVariables, numDesignVariables+1:numConstraints+numDesignVariables) = - conFuncGradient;
    lhsEqMtx(numDesignVariables+1:numConstraints+numDesignVariables, 1:numDesignVariables) =  conFuncGradient';

    % Create RHS-vector of the linear system of equations
    rhsEqVec = -[objFuncGradient; constraintValues];
    % systemRank = rank(lhsEqMtx,1e-6);
    
    % if systemRank < numDesignVariables+numConstraints
    % 
    %     lhsEqMtx =  lhsEqMtx + eye(size(lhsEqMtx)).*1e-5;
    % 
    % 
    % end

    % Solve the system for the step size(s) and lagrange multiplier(s).
    %infIDs = find(abs(sol(numDesignVariables+1:numDesignVariables+numConstraints))== Inf);

    % if systemRank < numDesignVariables+numConstraints
    % 
    %     lhsEqMtx(1:numDesignVariables, 1:numDesignVariables) = lagrangeHessian;
    %     lhsEqMtx(1:numDesignVariables, numDesignVariables+1:numConstraints+numDesignVariables) = - (conFuncGradient + 1e-2);
    %     lhsEqMtx(numDesignVariables+1:numConstraints+numDesignVariables, 1:numDesignVariables) =  (conFuncGradient' + 1e-2);
    % 
    %     % Create RHS-vector of the linear system of equations
    %     rhsEqVec = -[objFuncGradient; constraintValues];
    % 
    %     % if numRows <= numConstraints
    %     %     rowIDs = numDesignVariables+numConstraints-numRows+1:numDesignVariables+numConstraints;
    %     %     % Fill in the data of the LHS-matrix
    %     %     lhsEqMtx(rowIDs,:) = [];
    %     %     lhsEqMtx(:,rowIDs) = [];
    %     %     % Create RHS-vector of the linear system of equations
    %     %     rhsEqVec(rowIDs) = [];
    %     % 
    %     %     sol = lhsEqMtx \ rhsEqVec;
    %     % 
    %     %     solTemp = zeros(length(sol)+length(rowIDs),1);
    %     % 
    %     %     solTemp(1:rowIDs(1)-1) = sol(1:rowIDs(1)-1);
    %     % 
    %     % 
    %     %     sol = solTemp;
    %     % end
    % 
    % end

    sol = lhsEqMtx\rhsEqVec;

    % Check if the constraint in the next step can be neglected.
    lagrangeMultp = sol(numDesignVariables+1 : numDesignVariables+numConstraints);
    negativeLagrangeID = find(lagrangeMultp < 0);
    posLagrangeID = find(lagrangeMultp > 0);

    if ~isempty(negativeLagrangeID)

        lhsEqMtx = zeros(numDesignVariables + numConstraints-length(negativeLagrangeID));

        % Fill in the data of the LHS-matrix
        lhsEqMtx(1:numDesignVariables, 1:numDesignVariables) = lagrangeHessian;
        lhsEqMtx(1:numDesignVariables, numDesignVariables+1:numConstraints+numDesignVariables-length(negativeLagrangeID)) = - conFuncGradient(posLagrangeID);
        lhsEqMtx(numDesignVariables+1:numConstraints+numDesignVariables-length(negativeLagrangeID), 1:numDesignVariables) = conFuncGradient(posLagrangeID)';

        % Create RHS-vector of the linear system of equations
        rhsEqVec = -[objFuncGradient; constraintValues(posLagrangeID)];

        sol = lhsEqMtx \ rhsEqVec;
        
        solTemp = zeros(numDesignVariables+numConstraints,1);

        for ii=1:numDesignVariables

            solTemp(ii) = sol(ii);

        end
        for ii=numDesignVariables+1:numDesignVariables+numConstraints

            if any(negativeLagrangeID == (ii - (numDesignVariables)))
                
                solTemp(ii) = 0;

            else

                solTemp(ii) = sol(ii);

            end

        end
        sol = solTemp;
        
    end


    % Check if the propesed step lowers the penalty function.
    designVarTemp = designVariables + sol(1:numDesignVariables);
    [objValueTemp, objGradientTemp, ~] =  objFunctionHandle(designVarTemp);
    [conValuesTemp, conGradientTemp, ~] = conFuncHandle(designVarTemp);


    % Check if any of the constraints are violated
    violatedConID = find(conValuesTemp < 0);


    if  ~isempty(violatedConID)

        lagrangeMultpTemp = sol(numDesignVariables+1:numDesignVariables+numConstraints);
        penaltyFuncTemp = objValueTemp + transpose(lagrangeMultpTemp(violatedConID))*abs(conValuesTemp(violatedConID));

        if penaltyFuncTemp < penaltyFunc

            step = sol(1:numDesignVariables);
            lagrangeMultp = sol(numDesignVariables+1 : numDesignVariables+numConstraints);

        else

            step = sol(1:numDesignVariables);
            while penaltyFuncTemp >= penaltyFunc

                step = .5.*step;
                designVarTemp = designVariables + step;
                [objValueTemp, objGradientTemp, ~] =  objFunctionHandle(designVarTemp);
                [conValuesTemp, conGradientTemp, ~] = conFuncHandle(designVarTemp);
                penaltyFuncTemp = objValueTemp + transpose(lagrangeMultpTemp(violatedConID))*abs(conValuesTemp(violatedConID));
                fprintf(penaltyFuncTemp + "\n")

            end
        end

    else

        lagrangeMultpTemp = sol(numDesignVariables+1:numDesignVariables+numConstraints);
        penaltyFuncTemp = objValueTemp;

        if penaltyFuncTemp < penaltyFunc

            step = sol(1:numDesignVariables);
            lagrangeMultp = sol(numDesignVariables+1 : numDesignVariables+numConstraints);

        else

            while penaltyFuncTemp >= penaltyFunc

                step = .5.*step;
                designVarTemp = designVariables + step;
                [objValueTemp, objGradientTemp, ~] =  objFunctionHandle(designVarTemp);
                [~, conGradientTemp, ~] = conFuncHandle(designVarTemp);
                penaltyFuncTemp = objValueTemp;
                fprintf("3")

            end
        end

    end


    designVariables = designVariables + step;
    lagrangeMultp = lagrangeMultpTemp;


    if hessianCalc == "BFGS"


        a = objFuncGradient - lagrangeMultp*conFuncGradient;
        b = objGradientTemp - lagrangeMultp*conGradientTemp;

        gamma = b - a;
        delta = step;

        lagrangeHessian = lagrangeHessian + (gamma*gamma')/(gamma'*delta)...
            - (lagrangeHessian*(delta*delta')*lagrangeHessian)/(delta'*...
            lagrangeHessian*delta);
    end

iteration = iteration +1;

end

opt = designVariables;


end

