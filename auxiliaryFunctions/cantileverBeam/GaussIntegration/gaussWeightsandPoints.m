function [weights, quadraturePoints] = gaussWeightsandPoints(integrandOrder)
%% FUNCTION gaussWeightsandPoints
%   This function returns the Gauss quadrature points in 1D and their 
%   respective weights for a given integrand polynomial order.
%
%   Author(s): Deha Şen Köse, deha.koese@tum.de
%
%% Reference(s):
%
%   Carrera E., Cinefra M., Zappino E., Petrolo M. (2014). 
%   Finite Element Analysis of Structures Through Unified Formulation, 
%   Appendix A: Numerical Integration. John Wiley & Sons Ltd.
%   DOI:10.1002/9781118536643
%   
%
%% Input(s):
%
%   integrandOrder      : Polynomial order of the complete integrand.
%
%% Output(s):
%
%   weights             : The vector that contains the weights to the
%                         corresponding quadrature points.
%
%   quadraturePoints    : Locations of the quadrature points in 1D.
%
%% End of function definition - Code

% Find the correct number of required Gauss points to perform exact
% numerical integration.
nGP = ceil((integrandOrder+1)/2);

% Select the correct values from the dictionary below and return them.

if nGP == 1

    weights = 2;
    quadraturePoints = 0;

elseif nGP == 2

    weights = [1, 1];
    quadraturePoints = [-0.5773502691896257, 0.5773502691896257];

elseif nGP ==3

    weights = [0.8888888888888888, 0.5555555555555556, 0.5555555555555556];
    quadraturePoints = [0, -0.7745966692414834, 0.7745966692414834];

elseif nGP == 4

    weights = [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
    quadraturePoints = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];

elseif nGP == 5

    weights = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
    quadraturePoints = [0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];

elseif nGP == 6

    weights = [0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704];
    quadraturePoints = [0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521];

elseif nGP == 7
    
    weights = [0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697];
    quadraturePoints = [0, 0.4058451513773972, -0.4058451513773972, -0.7415311855993945, 0.7415311855993945, -0.9491079123427585, 0.9491079123427585];

else 

    error("The dictionary does not have quadrature points above 7! Please update it!")

end

end
%% End of code.