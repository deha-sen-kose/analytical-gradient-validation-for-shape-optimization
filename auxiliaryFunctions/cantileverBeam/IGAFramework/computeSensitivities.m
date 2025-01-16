function [desingVariableSensitivities, designHessian] = ...
         computeSensitivities...
         (displacementField, designVariablesDOFs, stiffnessSensitivities, ...
         forceSensitivities, stiffnessHessian)
%% FUNCTION computeSensitivities
%  
%   Computes ovearall sensitivities of the strain energy function using the
%   derivatives of the stiffness matrix and force vector.
%
%   Author(s)               : Deha Şen Köse, deha.koese@tum.de
%
%% Input(s):
%
%   displacementField       : Solution field of the IGA linear system of 
%                             equations. Namely the displacement vector.
%
%   designVariablesDOFs     : DOF numbering of the design variables.
%
%   stiffnessSensitivities  : Sensitivities of the stiffness matrix with
%                             respect to the ALL control point coordinates,
%                             and not just the design variables. In the
%                             future this implementation can be changed in
%                             order to reduce computational time and memory
%                             efficiency.
%
%   forceSensitivities      : Sensitivities of the force vector with
%                             respect to the ALL control point coordinates,
%                             and not just the design variables. In the
%                             future this implementation can be changed in
%                             order to reduce computational time and memory
%                             efficiency.
%
%   stiffnessHessian        : Hessian of the master stiffness matrix with
%                             respect to the design variables. It is a 
%                             numeric array  
%
%% Output(s):
%
%desingVariableSensitivities: Gradients of the strain energy function with
%                             respect to the design variables as stored in 
%                             increasing DOF order.
%
%
%   designHessian           : Hessian of the objective function with
%                             respect to the design variables. In this 
%                             thesis it is not validated and used.
%
%% End of function definition - Code


% Memory allocation for the gradient vector.
desingVariableSensitivities = zeros(length(designVariablesDOFs), 1);

% Looping through the desing variables and computing overall sensitivities
% of the strain energy function with respect to nodal coordinates.
for ii=1:length(desingVariableSensitivities)
    desingVariableSensitivities(ii) = 0.5*transpose(displacementField)*forceSensitivities{designVariablesDOFs(ii),1}...
        - 0.5*transpose(displacementField)*(stiffnessSensitivities{designVariablesDOFs(ii),1}*displacementField...
        - forceSensitivities{designVariablesDOFs(ii),1});

end

% Return the hessian if necessary. If the quadraticProgramming in 
% assembleMasterSystemAndSensitivitiesNum.m was set to false, the code in
% below will return zero matrices.
if nargin == 5

    designHessian = zeros(length(designVariablesDOFs));

    for ii=1:height(designHessian)
        for jj=1:width(designHessian)
            
            designHessian(ii,jj) =  designHessian(ii,jj) + transpose(displacementField)*stiffnessHessian{designVariablesDOFs(ii),designVariablesDOFs(jj)}...
                *displacementField;


        end
    end



end


end
% End of code