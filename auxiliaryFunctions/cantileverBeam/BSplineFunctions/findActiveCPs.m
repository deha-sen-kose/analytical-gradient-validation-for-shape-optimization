function [xiActive,etaActive] = findActiveCPs(knotVectorXi, knotVectorEta,...
    polOrderXi, polOrderEta)
%%  FUNCTION findActiveCPs
%
%   Function that finds the active basis functions (control points)
%   in both xi- and eta-directions, for a given configuration of a surface.
%
%   Author(s): Deha Şen Köse, deha.koese@tum.de
%
%% Input(s):
%
%     knotVectorXi: (m_xi,1) array of non-descending knots. 
%                   (m_xi = n_xi + p_xi + 1)
%
%    knotVectorEta: (m_eta,1) array of non-descending knots. 
%                   (m_eta = n_eta + p_eta + 1)
%
%       polOrderXi: polynomial order in xi-direction. must be >= 0.
%
%      polOrderEta: polynomial order in eta-direction. must be >=0.
%
%% Output(s):
%
%         xiActive: (Num_e_xi, active_id) array of active CP id(s) as 
%                   columns.
%
%        etaActive: (Num_e_eta, active_id) array of active CP id(s) as 
%                   columns.
%
%% End of function definition - Code

% Find the index of the elements in xi- and eta-directions.
elXi = findNonZeroKnotSpans(knotVectorXi);
elEta = findNonZeroKnotSpans(knotVectorEta);

% Get the number of elements in xi- and eta-directions.
numElXi = length(elXi);
numElEta = length(elEta);
%numElem = numElXi*numElEta;

% Compute the number of active contol points (B-Spline basis functions) in
% one element.
numActiveControlPointsXi = polOrderXi+1;
numActiveControlPointsEta = polOrderEta+1;
%numActiveControlPointsEl = numActiveControlPointsXi*numActiveControlPointsEta;

% Get the number of control points.
numControlPointsXi = length(knotVectorXi) - polOrderXi - 1;
numControlPointsEta = length(knotVectorEta) - polOrderEta - 1;

% Initialize output.
xiActive = zeros(numElXi,polOrderXi+1);

% Loop over the elements of xi.
for jj=1:numElXi

    % Create a counter for array indexing.
    xiActiveCounter = 1;

    % Loop over the control points. 
    for ii=1:numControlPointsXi
        
        % Check if the control point is active.
        if  ii <= elXi(jj) && ii+polOrderXi+1 >= elXi(jj)+1
            
            % Add the control point to the active set.
            xiActive(jj, xiActiveCounter) = ii;

            xiActiveCounter = xiActiveCounter + 1;

        end
    end
end

% Check for output validity.
if xiActiveCounter-1 ~= numActiveControlPointsXi

    error("Determining the active control points in xi-direction  went wrong!")
    
end


% Initialize output.
etaActive = zeros(numElEta,polOrderEta+1);

% Loop over the elements of eta.
for jj=1:numElEta

    % Create a counter for array indexing.
    EtaActiveCounter = 1;
    
    % Loop over the control points.
    for ii=1:numControlPointsEta
        
        % Check if the control point is active.
        if  ii <= elEta(jj) && ii+polOrderEta+1 >= elEta(jj)+1

            % Add the control point to the active set.
            etaActive(jj, EtaActiveCounter) = ii;

            EtaActiveCounter = EtaActiveCounter + 1;

        end
    end
end

% Check for output validity.
if EtaActiveCounter-1 ~= numActiveControlPointsEta

    error("Determining the active control points in eta-direction went wrong!")
    
end
   
end
%% End of code
