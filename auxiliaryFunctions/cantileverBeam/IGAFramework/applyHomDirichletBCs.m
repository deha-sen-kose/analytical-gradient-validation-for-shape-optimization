function [reducedStiffnessMtx, reducedForceVector] = applyHomDirichletBCs...
    (homDOFs, masterStiffnessMtx, forceVector)
%% FUNCTION applyHomDirichletBCs
%   
%   This function reduces the finite element system using the homogeneous
%   dirichlet boundary conditions.
%
%   Author(s)       : Deha Şen Köse, dehasenkose@gmail.com
%   
%% Input(s):
%
%   homDOFs             : Numbering of the homogeneous degrees of freedoom. In 
%                       this framework, the degrees of freedoom are in x- and
%                       y-displacement pairs. Also, the global numbering of
%                       dofs or rather control points are first incremented
%                       in y- and then x-direction. e.g.,
%                     
%                             *5 ---------*10 --------- *15 --------- *20
%                             *4 --------- *9 --------- *14 --------- *19
%                             *3 --------- *8 --------- *13 --------- *18
%                             *2 --------- *7 --------- *12 --------- *17
%                             *1 --------- *6 --------- *11 --------- *16
%
%   masterStiffnessMtx  : The master stiffness matrix that is wished to be
%                         reduced.
%
%   forceVector         : The force vector that is wished to be reduced.
%
%% Output(s):
%
%  reducedStiffnessMtx  : The reduced stiffness matrix. Size;
%                         ((N-homDOF)x(N-homDOF)) 
%
%  forceVector          : The reduced force vector. Size: ((N-homDOF)x1)   
%
%% End of function definition - Code

% Check for input 
if length(homDOFs) >= height(masterStiffnessMtx) || length(homDOFs) >= length(forceVector)

    error("Given number of homogeneous degrees of freedoom (%d) is larger than the system!", length(homDOFs))

elseif height(masterStiffnessMtx) ~= width(masterStiffnessMtx)

    error("Provided stffness matrix is not square!")
    
end

% Delete the rows and columns of master stiffness matrix that are related
% to homogeneous degrees of freedoom
masterStiffnessMtx(homDOFs, :) = [];
masterStiffnessMtx(:, homDOFs) = [];
reducedStiffnessMtx = masterStiffnessMtx;


% Delete the rows of force vector that are related to homogeneous degrees 
% of freedoom
forceVector(homDOFs) = [];
reducedForceVector = forceVector;


end
%% End of code