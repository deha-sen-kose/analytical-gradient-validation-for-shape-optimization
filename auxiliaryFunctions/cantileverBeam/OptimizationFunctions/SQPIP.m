function [opt, step, obj, con, objGrad] = SQPIP(objFunHandle, conFunHandle, designVariables, hessianCalcFlag)
%% FUNCTION SQPIP 
%   This function is an internal point constrained optimizer for non-linear
%   objective function and non-linear constraints.

% con are column vector
% obj are column vector
% conGrad are column vector of each row being a gradient_1 _2 _3

numDV = length(designVariables);

if hessianCalcFlag == "BFGS"


    lagHessian = eye(numDV);



    % Compute the objective and constraint values and gradients
    [~, objGrad, ~] = objFunHandle(designVariables);
    [con, ~, ~] = conFunHandle(designVariables);

    numCon = height(con);
    e = ones(numCon,1);


    % Initialize the slack as the constraint value at the current point
    slack = diag(con);
    slackVec = con;

    % Guess the system's size
    mux = norm(objGrad, Inf);
    nu = mux;
    %mux = 2;
    step = 1;
    muRed = mux;

    % Solve for lambda values such that the S*Lambda*e - mu*e = 0.
    % Approximate values here also work.
    lambda = transpose(transpose(slack) \ (0.01 + mux));
    lambdaVec = zeros(height(lambda),1);
    for ii=1:height(lambda)
        lambdaVec(ii) = lambda(ii,ii);
    end



    while norm(meritFun) > 1e-6

        lambda = diag(lambdaVec);
        slack = diag(slackVec);

        % Compute the objective and constraint values and gradients
        [obj, objGrad, ~] = objFunHandle(designVariables);
        [con, conGrad, ~] = conFunHandle(designVariables);

        % Get neccessary data to create left hand matrix
        j = conGrad';

        % Have an initial guess of mu
        mu = mux.*ones(numCon,1);


        %%%%%%%
        % lambda = 2;
        % lambdaVec = lambda;
        %%%%%%%%%
        % Put in the values for the left hand side matrix.
        lshMtxR1 = [lagHessian, zeros(numDV,numCon), transpose(-j)];

        lshMtxR2 = [zeros(numDV,numCon)', lambda, slack];

        lshMtxR3 = [j, -eye(numCon,numCon), zeros(numCon,numCon)];

        lshMtx = cat(1,cat(1,lshMtxR1,lshMtxR2),lshMtxR3);


        % Create the right hand side vector.
        rhsVec1 = objGrad - lambda*conGrad;
        rhsVec2 = slack*lambda*e - mu*e;

        rhsVec3 = con - slackVec;

        rhsVec = -cat(1,cat(1,rhsVec1,rhsVec2),rhsVec3);

        sol = lshMtx \ rhsVec;

        step = sol(1:numDV);
        slackStep = sol(numDV+1:numDV+numCon);
        lambdaStep = sol(numDV+numCon+1:numDV+2*numCon);

        % Check if the slack is positive after the step and plan the step size
        % accordingly.
        slackVecTemp = slackVec + slackStep;
        lambdaVecTemp = lambdaVec + lambdaStep;
        if any(slackVecTemp < 0)

            [minVal, id] = min(slackVec);
            alpha = (0.001 - minVal)/slackStep(id);

            slackStep = alpha*slackStep;
            step = alpha*step;

        end
        if any(lambdaVecTemp < 0)

            [minVal, id] = min(lambdaVec);
            alpha = (0.001 - minVal)/lambdaStep(id);

            lambdaStep = alpha*lambdaStep;

        end
        %stepTemp = step;
        %slackVecTemp = slackVec + slackStep;
        dvTemp = designVariables + step;
        lambdaVecTemp = lambdaVec + lambdaStep;

        % Compute the merit function
        [objTemp, objGradTemp, ~] = objFunHandle(dvTemp);
        [conTemp, conGradTemp, ~] = conFunHandle(dvTemp);


        [r, ~] = find(con < 0);
        if ~isempty(r)
            meritFun = obj + nu*abs(con);
        else
            meritFun = obj;
        end

        [r, ~] = find(conTemp < 0);
        if ~isempty(r)
            meritFunTemp = objTemp + nu*abs(conTemp(r));
        else
            meritFunTemp = objTemp;
        end

        
        while meritFun < meritFunTemp



            slackStep = .5*slackStep;
            step = .5*step;
            lambdaStep = .5.*lambdaStep;
            lambdaVecTemp = lambdaVec + lambdaStep;
            dvTemp = designVariables + step;

            [objTemp, objGradTemp, ~] = objFunHandle(dvTemp);
            [conTemp, conGradTemp, ~] = conFunHandle(dvTemp);

            [r, ~] = find(conTemp < 0);
            if ~isempty(r)
                meritFunTemp = objTemp + nu*abs(conTemp(r));
            else
                meritFunTemp = objTemp;
            end
            %fprintf(meritFunTemp + "\n")
        end


        designVariables = designVariables + step;
        lambdaVec = lambdaVec + lambdaStep;
        slackVec = slackVec + slackStep;
        %mux = mux*log(norm(slackVec, Inf));



        if hessianCalcFlag == "BFGS"


            a = objGrad - lambdaVec*conGrad;
            b = objGradTemp - lambdaVec*conGradTemp;

            gamma = b - a;
            delta = step;
            if abs(norm(gamma, Inf)) > 1e-6 && abs(norm(delta)) > 1e-6
            

            % if abs(norm(gamma, Inf)) < mux*1e-6
            %     fprintf("Function ended in optimality measure!")
            %     break
            % end

            lagHessian = lagHessian + (gamma*gamma')/(gamma'*delta)...
                - (lagHessian*(delta*delta')*lagHessian)/(delta'*...
                lagHessian*delta);
            end

        end

        mux = mux/muRed;


    end

    opt = designVariables;

end

[obj, objGrad, ~] = objFunHandle(opt);
[con, conGrad, ~] = conFunHandle(opt);


end

