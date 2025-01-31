function controlPointsSurface = generateSurfaceCPs(controlPointsX, ...
    controlPointsY, reshapeFlag)
%% FUNCTION generateSurfaceCPs
%   Generates symmetrical surface control points from control points of
%   given curves.
%   
%   Author(s)       : Deha Şen Köse, dehasenkose@gmail.com
%
%% Input(s):
%
%   controlPointsX      : Control points in x-direction. (nx2) array where
%                         first column is x- and second column is 
%                         y-coordinate of the curve. For a rectangular 
%                         surface, second column must be all zeros.
%
%   controlPointsY      : Control points in y-direction. (nx2) array where
%                         first column is x- and second column is 
%                         y-coordinate of the curve. For a rectangular 
%                         surface, first column must be all zeros.
%
%   reshapeFlag         : A boolean for the type of the output. If true,
%                         output is a (n_cpx2) array and if false a cell
%                         array as the shape of the control points on the
%                         surface.
%
%% Output(s):       
%
% controlPointsSurface  : Either cell (as the shape of the control points 
%                         on the surface) or a (n_CPx2) array of the 
%                         control point coordinates where indexing first 
%                         increases in y- or eta-direction. If reshapeFlag
%                         is false then it is a cell array, if true then it
%                         is a (n_CPx2) array of the control point 
%                         coordinates where indexing first increases in y- 
%                         or eta-direction. 
%
%% End of function definition - Code

% Initialize output.
controlPointsSurface = cell(length(controlPointsX),length(controlPointsY));

% Loop over the input and add them.
for ii=1:height(controlPointsX)
    for jj=1:height(controlPointsY)
     
            controlPointsSurface{ii,jj} = [controlPointsX(ii,1), controlPointsY(jj,2)];
    end
end


% Check for the flag.
if reshapeFlag

    % Initialize output.
    controlPointsSurfaceReshape = zeros(width(controlPointsSurface)*height(controlPointsSurface),2);
    counter = 1;

    % Loop over the surface control points.
    for ii=1:height(controlPointsSurface)
        for jj=1:width(controlPointsSurface)
            
            % Add the point coordinates to the rows of a matrix.
            controlPointsSurfaceReshape(counter,:) = controlPointsSurface{ii,jj};
            counter = counter + 1;

        end
    end
    
    controlPointsSurface = controlPointsSurfaceReshape;

end

end
%% End of code.