function [BSplineBasisFunctionEvaluation, BSplineDerivativeEvaluation] = ...
    computeBSplineBasisFunctionAndDerivatives(i, p, curveParameter, knotVector)
%% FUNCTION computeBSplineBasisFunctionAndDerivatives
%   This function computes the B-Spline basis functions and their first 
%   derivatives. All values are numerically computed.      
%
%   Author(s)       : Deha Şen Köse, dehasenkose@gmail.com
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
%   curveParameter  : The curve parameter where the B-Spline function is
%                     wished to be evaluated or differentiated.
%
%   knotVector      : Knot vector as a non-descending array.
%
%% Output(s):
%
%   BSplineBasisFunctionEvaluation  : Numerical evaluatin of the B-Spline
%                                     basis function.
%
%   BSplineDerivativeEvaluation     : Numerical evaluation of the B-Spline
%                                     basis function derivative.
%
%% End of function definition - Code

% Check for data validity
if length(knotVector) < i + p + 1

    error("Shape function index (%d) is out of range for the given configuration!", i)

elseif i <= 0

    error("Shape function numbering (%d) is wrong. It must start from 1!", i)

elseif p < 0

    error("Polynomial degree (%d) cannot be lower than zero!", p)

elseif curveParameter < knotVector(1) || curveParameter > knotVector(end)

    error("Given curve parameter (%d) is outside of the knot vector!", curveParameter)

elseif ~issorted(knotVector)

    error("Given knot vector does not satisfy the non-descending order character!")

end


% Define a tolerance for zero denominator.
tolerance = 1e-6;

% Initialize outputs.
BSplineBasisFunctionEvaluation = 0.0;

BSplineDerivativeEvaluation = 0.0;

% Recursive formulae for both function evaluation and derivatives.
if p == 0
    
    if (curveParameter >= knotVector(i) && curveParameter < knotVector(i+1))...
            || (abs(curveParameter - knotVector(i+1)) < tolerance && abs(curveParameter - knotVector(end)) < tolerance)
    
        BSplineBasisFunctionEvaluation = 1.0;

        BSplineDerivativeEvaluation = 0.0;

    end

else

    if abs(knotVector(i+p)-knotVector(i)) > tolerance
        
        BSplineBasisFunctionEvaluation = BSplineBasisFunctionEvaluation +...
        (curveParameter - knotVector(i))/(knotVector(i+p)-knotVector(i))*computeBSplineBasisFunctionAndDerivatives(i, p-1, curveParameter, knotVector);

        BSplineDerivativeEvaluation = BSplineDerivativeEvaluation +...
        p/(knotVector(i+p)-knotVector(i))*computeBSplineBasisFunctionAndDerivatives(i, p-1, curveParameter, knotVector);
           
    end

    if abs(knotVector(i+p+1)-knotVector(i+1)) > tolerance
        
        BSplineBasisFunctionEvaluation = BSplineBasisFunctionEvaluation +...
        (knotVector(i+p+1)-curveParameter)/(knotVector(i+p+1)-knotVector(i+1))*computeBSplineBasisFunctionAndDerivatives(i+1, p-1, curveParameter, knotVector);

        BSplineDerivativeEvaluation = BSplineDerivativeEvaluation -...
        p/(knotVector(i+p+1)-knotVector(i+1))*computeBSplineBasisFunctionAndDerivatives(i+1, p-1,curveParameter, knotVector);

    end

end


end
%% End of code.
