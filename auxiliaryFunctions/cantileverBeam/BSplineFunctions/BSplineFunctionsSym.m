function symBSplineBasisFunc = BSplineFunctionsSym(i, p, knotVector, direction)
%% FUNCTION BSplineFunctionsSym
%   This function computes the B-Spline basis functions. All values are 
%   symbolically computed. This piece of code is not used in generation of
%   the results in the thesis due to efficiency reasons.
%
%   Author(s): Deha Şen Köse, deha.koese@tum.de
%
%% Reference(s):
%
%   Prochazkova J. (n.d.). Derivative of B-Spline Function.
%   25. KONFERENCE O GEOMETRII A POCITACOVEGRAFICE
%   https://mat.fsv.cvut.cz/gcg/sbornik/prochazkova.pdf
%
%% Input(s):
%
%   i               : The index of the B-Spline function. In this framework, i
%                     starts from 1. i has to be integer.
%
%   p               : The polynomial order of the B-Spline basis function.
%
%   direction       : A string for the knot vector variable, "xi" or "eta".
%
%   knotVector      : Knot vector as a non-descending array.
%
%% Output(s):
%
%   symBSplineBasisFunc  : Symbolic evaluatin of the B-Spline
%                          basis function.
%
%% End of function definition - Code

% Check for data validity
if length(knotVector) < i + p + 1

    error("Too many shape functions (%d) for the given configuration!", i)

elseif i <= 0

    error("Shape function numbering (%d) is wrong. It must start from 1!", i)

elseif p < 0

    error("Polynomial degree (%d) cannot be lower than zero!", p)

end

% Create xi and eta symbolic parameters.
if direction == "xi"

    syms xi
    dir = xi;

elseif direction == "eta"

    syms eta
    dir = eta;
    
end

% Define a tolerance for zero denominator.
tolerance = 1e-5;
%tolerance2 = 0.01;

% Initialize outputs.
symBSplineBasisFunc = 0.0;

% Recursive formulae for both function evaluation and derivatives.
if p == 0

    symBSplineBasisFunc = piecewise((knotVector(i) <= dir) & (dir < knotVector(i+1)), 1.0, 0);
    
else

    if abs(knotVector(i+p)-knotVector(i)) > tolerance

        symBSplineBasisFunc = symBSplineBasisFunc +...
            (dir - knotVector(i))/(knotVector(i+p)-knotVector(i))*BSplineFunctionsSym(i, p-1, knotVector, direction);

    end

    if abs(knotVector(i+p+1)-knotVector(i+1)) > tolerance

        symBSplineBasisFunc = symBSplineBasisFunc +...
            (knotVector(i+p+1)- dir)/(knotVector(i+p+1)-knotVector(i+1))*BSplineFunctionsSym(i+1, p-1, knotVector, direction);

    end

end


end
%% End of code.

