function cpmtxnew = symmetryFun(cpmtx, cpmtxk, numCPy, symtrynum)
%%  Function symmetryFun
%
%   This function symmetrically moves below the mid x-line control points
%   together with the above mid x-line control points (parameter linking).
%   Done for regularizing the shape and obtaining symmetrical shapes after
%   optimization.
%
%   Author(s)               : Deha Şen Köse, dehasenkose@gmail.com
%
%%  Input(s):
%
%   cpmtx                   : (nx2) array of the coordinates of the control
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
%
%   cpmtxk                  : The control point matrix in k-th iteration 
%                             step. Or can be used to mae symmetric shapes 
%                             outside of the optimization loop.    
%
%   numCPy                  : The number of control points in y-direction.
%
%   symtrynum               : The vector containing the numbering of
%                             control points that are wished to be
%                             symmetrically moved with the corresponding
%                             upper ones.
%
%%  Output(s):
%
%   cpmtxnew                : The new control point matrix after symmetry
%                             is applied. The indexing is kept the same as
%                             explained above.
%
%% End of function definition - Code.

% Check if there are even number of control points in y-direction.
% Depending on this, create a pattern vector.
if mod(numCPy,2) == 0
    ptrn = [numCPy-1,numCPy-2-1];
else
    ptrn = numCPy-1;
end

% Get the difference between the original control points and their current
% step.
delta = gsubtract(cpmtxk,cpmtx);

% Initialize the new cp matrix using the original one.
cpmtxnew = cpmtxk;

% Loop through the control point indices that are wihed to be made
% symmetric.
for ii=1:length(symtrynum)

    % Switch for the pattern.
    if length(ptrn) > 1
        % Apply the symmetries.
        if mod(symtrynum(ii),2) == 0
            cpmtxnew(symtrynum(ii),1) = cpmtxnew(symtrynum(ii),1) + delta(symtrynum(ii)+ptrn(2),1);
            cpmtxnew(symtrynum(ii),2) = cpmtxnew(symtrynum(ii),2) - delta(symtrynum(ii)+ptrn(2),2);
        else
            cpmtxnew(symtrynum(ii),1) = cpmtxnew(symtrynum(ii),1) + delta(symtrynum(ii)+ptrn(1),1);
            cpmtxnew(symtrynum(ii),2) = cpmtxnew(symtrynum(ii),2) - delta(symtrynum(ii)+ptrn(1),2);

        end
    else
            cpmtxnew(symtrynum(ii),1) = cpmtxnew(symtrynum(ii),1) + delta(symtrynum(ii)+ptrn(1),1);
            cpmtxnew(symtrynum(ii),2) = cpmtxnew(symtrynum(ii),2) - delta(symtrynum(ii)+ptrn(1),2);
    end

end

end
% End of code.