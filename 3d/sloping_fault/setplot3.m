%  setplot3.m
%  called in plotclaw3 before plotting to set various parameters

OutputDir = '_output';
setviews  % set viewpoints so that view(xSlice), for example, can be used.
%setopengl;

PlotType = 1;                % type of plot to produce:
			     % 1 = pcolor on slices (with optional contours)
			     % 2 = contour lines in 3d on transparent slices
			     % 3 = Schlieren plot on slices
			     % 4 = scatter plot of q vs. r
			     % 5 = isosurface plot (at given levels)

mq = 9;                      % which component of q to plot
UserVariable = 0;            % set to 1 to specify a user-defined variable
UserVariableFile = '';      % name of m-file mapping data to q
MappedGrid = 1;              % set to 1 if mapc2p.m exists for nonuniform grid
MaxFrames = 1000;            % max number of frames to loop over
MaxLevels = 6;               % max number of AMR levels

PlotData =  [1 1 1 1 1 1];       % Data on refinement level k is plotted only
			         % if k'th component is nonzero
PlotGrid =  [1 1 1 1 1 1];       % Plot grid lines on each level?
PlotGridEdges =  [0 0 0 0 0 0];  % Plot edges of patches of each grid at
                                 % this level on slices?
PlotCubeEdges = [0 0 0 0 0 0];   % Plot edges of cube of refinement patch at
                                 % this level?


% ContourValues is a vector of contour lines that can be used with
% PlotType = 1,2,3.   Empty ==> no contour lines drawn

ContourValues = [];   % draw contour lines from 'afterframe.m'

% The next three parameters are vectors of x,y,z coordinates of 2d slices
% to be displayed for PlotType = 1,2,3.

  xSliceCoords = [];
  ySliceCoords = [0.0];
  zSliceCoords = [0.0];

  IsosurfValues    =  [];     % Plot surfaces at q = surfValue(i).

  IsosurfAlphas    =  [];     % Transparency of each surface
                                          % (0=clear; 1=opaque)
                                          % NOTE: Your system must be able to
                                          % use the OpenGL Renderer.

  IsosurfColors    = 'q';      % Colors for each surface.
