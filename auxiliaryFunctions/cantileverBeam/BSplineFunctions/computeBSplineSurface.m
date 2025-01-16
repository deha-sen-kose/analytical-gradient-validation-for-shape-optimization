function surface = computeBSplineSurface(knotVectorXi, polOrderXi,...
    knotVectorEta, polOrderEta, controlPointsSurface, samplePoints,...
    plotFlag, dvDOF)
%% FUNCTION computeBSplineSurface
%
%   Computes B-Spline surface by the tensor product of two 2D curves. Can
%   also be used to plot the computed 2D surface.
%
%   Author(s)               : Deha Şen Köse, deha.koese@tum.de
%
%% Reference(s):
%
%   Kollmannsberger S. (2023). Computation in Engineering II Lecture Notes.
%
%% Input(s):
%
%   knotVectorXi            : Knot vector in Xi direction. It has to be in 
%                             the intervals of [0, 1].
%
%   knotVectorEta           : Knot vector in Eta direction. It has to be in
%                             the intervals of [0, 1].
%
%   polOrderXi              : Polynomial order in xi direction.
%
%   polOrderEta             : Polynomial order in eta direction.
%
%   controlPointsSurface    : (nx2) array of the coordinates of the control
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
%   samplePoints            : (1x2) array containing the number of sampling
%                             points in x- and y-direction.
%
%   plotFlag                : A boolean for the plotting of the B-spline
%                             surface. If true, the surface is plotted, if
%                             false only the surface points are computed
%                             and returned.
%
%   dvDOF                   : Degrees of freedom of the design variables as
%                             a vector. This way the design variable
%                             control points are plotted in green and the
%                             ones that are not desing variables are
%                             plotted red as also shown in the thesis. If
%                             there is no design variable, this can be
%                             given as NaN value, so that all of the
%                             control points are drawn in red.
%
%% Output(s):
%
%   surface                 : Surface points as a 2D cell array where first
%                             matrix contains the x- and second matrix 
%                             contains y-coordinates of the surface. Since 
%                             coordinate points are in the same plane in 
%                             z-axis, the surface is assumed to be 2D.
%
%% End of function definition - Code


% Check for the input knot vectors
if knotVectorXi(1) ~=0 || knotVectorXi(end) ~=1 || knotVectorEta(1) ~=0 || knotVectorEta(end) ~=1

    error("Provided knot vectors are not between 0 and 1. Please update the vector(s) or this code!!!")

end

% Get the number of control points
numControlPointsX = length(knotVectorXi) - polOrderXi - 1;
numControlPointsY = length(knotVectorEta) - polOrderEta - 1;

% Allocate memory for surface points
surface = cell(2,1);

% Loop through the dimensionality. In membrane action the problem is
% assumed to be 2D.
for iField=1:2

    surface{iField} = zeros(samplePoints(1),samplePoints(2));

    % Loop through the realization points.
    for iSample=1:samplePoints(1)

        for jSample=1:samplePoints(2)

            % Get the realization points in xi and eta.
            s = (iSample -1) / (samplePoints(1) - 1);
            r = (jSample -1) / (samplePoints(2) - 1);


            % Loop for the double Summation
            for ii=1:numControlPointsX
                for jj=1:numControlPointsY

                    % Calculate the basis functions
                    [Ns, ~] = computeBSplineBasisFunctionAndDerivatives(ii,...
                        polOrderXi, s, knotVectorXi);
                    [Nr, ~] = computeBSplineBasisFunctionAndDerivatives(jj,...
                        polOrderEta, r, knotVectorEta);

                    % Compute the surface points
                    surface{iField}(iSample, jSample) = surface{iField}...
                        (iSample, jSample) + Ns*Nr*controlPointsSurface...
                        ((ii-1)*numControlPointsY+jj,iField);

                end
            end


        end
    end

end

% Plot the surface if wanted
if plotFlag
    if isnan(dvDOF)
        figure
        surf(surface{1},surface{2}, zeros(size(surface{1})),FaceColor=[0,0.396,0.741])
        hold on
        scatter3(controlPointsSurface(:,1),controlPointsSurface(:,2),zeros(size(controlPointsSurface(:,1))),"filled",MarkerFaceColor="r")
        text(controlPointsSurface(:,1)+0.25,controlPointsSurface(:,2)+0.25,zeros(size(controlPointsSurface(:,1)))+0.25,...
            cellstr(num2str(transpose(1:height(controlPointsSurface)))));
    else
        numDV = length(dvDOF);
        r = zeros(numDV,1);
        % Loop through the design variables to update the surface control points
        % matrix respecting the control polygons.
        for ii=1:numDV

            % Get the DOF of the current desing variable.
            dof = dvDOF(ii);
            % Find the row of that design variable in the surface control point
            % matrix.
            r(ii) = ceil(dof/2);
        end
        figure
        surf(surface{1},surface{2}, zeros(size(surface{1})),FaceColor=[0,0.396,0.741])
        hold on
        scatter3(controlPointsSurface(r,1),controlPointsSurface(r,2),zeros(size(controlPointsSurface(r,1))),"filled",MarkerFaceColor="g")
        text(controlPointsSurface(:,1)+0.25,controlPointsSurface(:,2)+0.25,zeros(size(controlPointsSurface(:,1)))+0.25,...
        cellstr(num2str(transpose(1:height(controlPointsSurface)))));
        hold on
        controlPointsSurfaceCopy = controlPointsSurface;
        controlPointsSurfaceCopy(r,:) = [];
        scatter3(controlPointsSurfaceCopy(:,1),controlPointsSurfaceCopy(:,2),zeros(size(controlPointsSurfaceCopy(:,1))),"filled",MarkerFaceColor="r")

    end
end

end
%% End of code
