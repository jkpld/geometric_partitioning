function options = setDeclumpOptions(varargin)
% SETDECLUMPOPTIONS  Set options needed for declumping.
%
% Input can be structure array or parameter value pairs. Options not set
% will be given default values. To see the default values, look at the
% output with no inputs:
%   defaultValues = setDeclumpOptions()
%
% Available options:
%
% Max_Radius: (0, Inf) , [pixels]
%   The maximum distance a boundry vertex can be from a center to be
%   assigned to it.
%
% Min_Angle:  [-1, 1] , [unitless]
%   The minimum *dot product* between the line formed between a boundry
%   vertex and a center and the unit normal (pointing inward) at the
%   boundary vertex. -- The angle could be found by acos(Min_Angle).
%
% Wigner_Seitz_Radius: (0, Inf) , [pixels]
%   The effective size of each particle in the simulation. This sets the
%   density of the particles used. The approximate number of particles used
%   in a particular object will be area(object)/(pi* r_s^2), where r_s is
%   the Wigner_Seitz_Radius.
%
% Initial_Speed: (-Inf, Inf) , [pixels/time]
%   The initial speed of the particles in the simulation.
%
% Point_Selection_Method: {'Random','UniformRandom','curvatureRandom', 
%                          'curvatureUniformRandom'}
%   The method used to initialize the particle locations.
%
% Potential_Depth: (-Inf,0] , [arb. units]
%   The depth of the potential.
%
% Potential_Minimum_Location: (0, Inf) , [pixels]
%   The location of the potential minimum.
%
% Potential_Extent: (0, Inf) & > Potential_Minimum_Location , [pixels]
%   The radius at which the potential goes from attractive to repulsive.
%
% Min_Interior_Angle: [0, 180] , [degrees]
%   The minimum interior angle for the triangulation of the centers.
%
% Max_Interior_Angle: [0, 180] & > Min_Interior_Angle , [degrees]
%   The maximum interior angle for the triangulation of the centers.
%
% Search_Radius: [0, Inf) & integer , [pixels]
%   The radius of the region over which each cut vertex is optimized over.
%
% Minimum_Hole_Size: [0, Inf) , [pixels^2]
%   The minimum hole size allowed in the mask.
%
% Curvature_Smoothing_Size: [0, Inf) & integer , [pixels]
%   The standard deviation of the gaussians used for smoothing the boundary
%   and calculating the curvature.
%
% Use_GPU: boolean
%   Determines if a GPU will be used to speed up calculation.
%
% Use_Parallel: boolean
%   Determines if multiple CPUs will be sued to speed up calculation.
%
% Debug: boolean
%   Determines if additional information will be returned from each
%   function.
%
% See also DECLUMPNUCLEI COMPUTEOBJECTCENTERS

% James Kapaldo
% 2016-10-28

ip = inputParser;
ip.FunctionName = 'set_declumpOptions';

% Boundary Assignment
addParameter(ip,'Max_Radius', 35, @(x)validateattributes(x,{'double'},{'scalar','positive','real','finite'}));
addParameter(ip,'Min_Angle', 0.5, @(x)validateattributes(x,{'double'},{'scalar','>=',-1,'<=',1,'real'}));

% Object center calculation (particle simulation)
addParameter(ip,'Wigner_Seitz_Radius', 10, @(x)validateattributes(x,{'double'},{'positive','scalar','real','finite'}));
addParameter(ip,'Initial_Speed', 0.01, @(x)validateattributes(x,{'double'},{'scalar','real','finite'}));
addParameter(ip,'Point_Selection_Method', 'curvatureUniformRandom', @(x)validatestring(x,{'Random','UniformRandom','curvatureRandom','curvatureUniformRandom'}));
addParameter(ip,'Potential_Depth', -1, @(x)validateattributes(-x,{'double'},{'scalar','nonnegative','real','finite'}));
addParameter(ip,'Potential_Minimum_Location', 2, @(x)validateattributes(x,{'double'},{'scalar','positive','real','finite'}));
addParameter(ip,'Potential_Extent', 15, @(x)validateattributes(x,{'double'},{'scalar','positive','real','finite'}));

% Triangulation
addParameter(ip,'Min_Interior_Angle', 20, @(x)validateattributes(x,{'double'},{'scalar','>=',0,'<=',180}));
addParameter(ip,'Max_Interior_Angle', 110, @(x)validateattributes(x,{'double'},{'scalar','>=',0,'<=',180}));

% Optimization
addParameter(ip,'Search_Radius', 7, @(x)validateattributes(x,{'double'},{'scalar','integer','nonnegative','real','finite'}));

% Mask
addParameter(ip,'Minimum_Hole_Size', 50, @(x)validateattributes(x,{'double'},{'scalar','nonnegative','real','finite'}));

% Curvature
addParameter(ip,'Curvature_Smoothing_Size', 2, @(x)validateattributes(x,{'double'},{'scalar','nonnegative','integer','real','finite'}));

% Computation
addParameter(ip,'Use_GPU', false, @(x) (x==1) || (x==0));
addParameter(ip,'Use_Parallel', false, @(x) (x==1) || (x==0));

% Debug
addParameter(ip,'Debug', false, @(x) (x==1) || (x==0));

parse(ip,varargin{:})


options.Max_Radius                  = ip.Results.Max_Radius;                 % [pixels]
options.Min_Angle                   = ip.Results.Min_Angle;                  % [dot product] : angle would be acos(Min_Angle)
options.Wigner_Seitz_Radius         = ip.Results.Wigner_Seitz_Radius;        % [pixels]
options.Initial_Speed               = ip.Results.Initial_Speed;              % [pixels/time]
options.Point_Selection_Method      = ip.Results.Point_Selection_Method;     % {'Random','UniformRandom','curvatureRandom','curvatureUniformRandom'}
options.Potential_Depth             = ip.Results.Potential_Depth;            % [arb. units.]
options.Potential_Minimum_Location  = ip.Results.Potential_Minimum_Location; % [pixels]
options.Potential_Extent            = ip.Results.Potential_Extent;           % [pixels]
options.Min_Interior_Angle          = ip.Results.Min_Interior_Angle;         % [degrees]
options.Max_Interior_Angle          = ip.Results.Max_Interior_Angle;         % [degrees]
options.Search_Radius               = ip.Results.Search_Radius;              % [pixels]
options.Minimum_Hole_Size           = ip.Results.Minimum_Hole_Size;          % [pixels^2]
options.Curvature_Smoothing_Size    = ip.Results.Curvature_Smoothing_Size;   % [pixels]
options.Use_GPU                     = ip.Results.Use_GPU;                    % [boolean]
options.Use_Parallel                = ip.Results.Use_Parallel;               % [boolean]
options.Debug                       = ip.Results.Debug;                      % [boolean]

if options.Potential_Extent <= options.Potential_Minimum_Location
    error('setDeclumpOptions:badParameterValues', 'Potential_Extent must be larger than Potential_Minimum_Location.')
end

if options.Max_Interior_Angle <= options.Min_Interior_Angle
    error('setDeclumpOptions:badParameterValues', 'Max_Interior_Angle must be larger than Min_Interior_Angle.')
end


% % Create default options structure
% % Parameters -------------------------------------------------------------
%                   
% % Boundary Assignment
% options.Max_Radius = 35;                                    % [pixels]
% options.Min_Angle = 0.5;                                    % [dot product] : angle would be acos(Min_Angle)
% 
% % Object center calculation (particle simulation)
% options.Wigner_Seitz_Radius = 10;                           % [pixels]
% options.Initial_Speed = 0.01;                               % [pixels/time]
% options.Point_Selection_Method = 'curvatureUniformRandom';  % {'Random','UniformRandom','curvatureRandom','curvatureUniformRandom'}
% options.Potential_Depth = -1;                               % [arb. units.]
% options.Potential_Minimum_Location = 2;                     % [pixels]
% options.Potential_Extent = 15;                              % [pixels]
% 
% % Triangulation
% options.Min_Interior_Angle = 20;                            % [degrees]
% options.Max_Interior_Angle = 110;                           % [degrees]
% 
% % Optimization
% options.Search_Radius = 7;                                  % [pixels]
% 
% % Mask
% options.Minimum_Hole_Size = 50;                             % [pixels^2]
% 
% % Curvature
% options.Curvature_Smoothing_Size = 2;                       % [pixels]
% 
% % Computation
% options.Use_GPU = true;                                     % [boolean]
% options.Use_Parallel = false;                               % [boolean]
% 
% % Debug : debuging outputs extra potentially useful information, this takes
% % extra memory
% options.Debug = true;                                       % [boolean]


% validFieldNames = {'Max_Radius';
%                    'Min_Angle';
%                    'Wigner_Seitz_Radius';
%                    'Initial_Speed';
%                    'Point_Selection_Method';
%                    'Potential_Depth';
%                    'Potential_Minimum_Location';
%                    'Potential_Extent';
%                    'Min_Interior_Angle';
%                    'Max_Interior_Angle';
%                    'Search_Radius';
%                    'Minimum_Hole_Size';
%                    'Curvature_Smoothing_Size';
%                    'Use_GPU';
%                    'Use_Parallel';
%                    'Debug'};

end