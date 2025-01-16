function [displacementField] = solveIGASystem(masterStiffnessMtx, forceVector)
%%  Function solveIGASystem
%
%   This function solves the Ku=F FEA equation and returns the displacement
%   field.
%
%   Author(s)               : Deha Şen Köse, deha.koese@tum.de
%
%%  Input(s):
%
%   masterStiffnessMtx      : The reduces master stiffness matrix.
%
%   forceVector             : The reduces master force vector.
%
%%  Output(s):
%
%   displacementField       : The displacement vector.
%
%% End of function definition - Code.


% Check for input
if height(masterStiffnessMtx) ~= width(masterStiffnessMtx)

    error("Given stiffness matrix is not square!")

elseif length(forceVector) ~= height(masterStiffnessMtx)

    error("Linear system of equations are inconsistent by size!")

end

% Check for symmetry of the stiffness matrix
% If not symmetric in Matlab's precision, then
if ~issymmetric(masterStiffnessMtx)
    
    % Initialize a counter for the non-symmetrical values.
    nonSymmetricCounter = 0;

    % Loop through the columns and rows of the stiffness matrix.
    for ii=width(masterStiffnessMtx)
        for jj=1:height(masterStiffnessMtx)
            
            % If the difference of the symmetric values are greater than
            % zero
            if abs(masterStiffnessMtx(ii,jj) - masterStiffnessMtx(jj,ii)) > 1e-6
               
                % Increase the counter.
                nonSymmetricCounter = nonSymmetricCounter +1;

            end

        end
    end

    % If the counter is bigger than zero, than something has gone wrong.
    if nonSymmetricCounter > 0

        error("Provided Stiffness Matrix is not symmetric in structural ..." + ...
            "engineering precision!! Please check your calculations!")

    else

        % To speed up the solver.
        masterStiffnessMtx = (masterStiffnessMtx + transpose(masterStiffnessMtx))/2;

    end

end


% Solve the system
displacementField = masterStiffnessMtx\forceVector;

end
% End of code.