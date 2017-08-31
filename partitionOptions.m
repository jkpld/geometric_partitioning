classdef partitionOptions
% PARTITIONOPTIONS  Set options needed for partitioning.
%
% Input can be structure array or parameter value pairs. Options not set
% will be given default values. To see the default values, look at the
% output with no inputs:
%   defaultValues = partitionOptions()
%
% paritionOptions Properties:
%
% Max_Radius - The maximum distance a boundry vertex can be from a center
%   to be assigned to it.
%   (0, Inf) , [pixels]
%
% Min_Angle - The minimum *dot product* between the line formed between a
%   boundry vertex and a center and the unit normal (pointing inward) at 
%   the boundary vertex. -- The angle could be found by acos(Min_Angle).
%   [-1, 1] , [unitless]
%
% Min_Interior_Angle - The minimum interior angle for the triangulation of
%   the centers.
%   [0, 180] , [degrees]
%
% Max_Interior_Angle - The maximum interior angle for the triangulation of
%   the centers.
%   [0, 180] & > Min_Interior_Angle , [degrees]
%
% Search_Radius - The radius of the region over which each cut vertex is
%   optimized over.
%   [0, Inf) & integer , [pixels]
%
% Minimum_Hole_Size - The minimum hole size allowed in the mask.
%   [0, Inf) , [pixels^2]
%
% Use_GPU - Determines if a GPU will be used to speed up calculation.
%   logical
%
% Use_Parallel - Determines if multiple CPUs will be used to speed up
%   calculation.
%   logical
%
% Debug -
%   Determines if additional information will be returned from each
%   function.
%   logical
%
% Object_Of_Interest - The index of an object of interest to declump. If an
%   object is giving an error. Then set this property to the object index
%   and set Debug to true. Plots showing intermediate steps of the
%   calculation will be shown that should help you debug the problem. Leave
%   this property empty for normal use.
    
    properties
        Max_Radius                   = 35;
        Min_Angle                    = 0.5;
                
        Min_Interior_Angle           = 20;
        Max_Interior_Angle           = 110;
        
        Search_Radius                = 7;
        
        % 2D mask clean up
        Minimum_Hole_Size            = 50;
        
%         Curvature_Smoothing_Size     = 2;
        
        Use_GPU                      = false;
        Use_Parallel                 = false;
        
        Debug                        = false;
        
        Object_Of_Interest           = [];
    end
    
    properties (Hidden)

        Use_ConvexHull = true;
        
        % These next two parameters would make apearence in functions
        % relating to computing boundary curvature.
        Curvature_Smoothing_Size     = 2; % Not used
        Curvature_Max_Radius         = 35; % Not used 
        
        Debug_plot = false;
    end
    
    methods
        function options = partitionOptions(varargin)
            
            % Now assign any properties given.
            if nargin > 0
                if numel(varargin) == 1 && isstruct(varargin{1})
                    op = varargin{1};
                    fd = fieldnames(op)';
                    for fieldname = fd
                        options.(char(fieldname)) = op.(char(fieldname));
                    end
                else
                    if ~mod(nargin+1,2)
                        error('paritionOptions:badInput','Input must be structure with properties as fiels, or property/value list with even number of elements.')
                    end
                    for i = 1:2:numel(varargin)
                        options.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end
        
        function obj = set.Max_Radius(obj,value)
            validateattributes(value,{'double'},{'scalar','positive','real','finite'})
            obj.Max_Radius = value;
        end
        
        function obj = set.Min_Angle(obj,value)
            validateattributes(value,{'double'},{'scalar','>=',-1,'<=',1,'real'})
            obj.Min_Angle = value;
        end
        
        function obj = set.Min_Interior_Angle(obj,value)
            validateattributes(value,{'double'},{'scalar','>=',0,'<=',180})
            if value > obj.Max_Interior_Angle %#ok<MCSUP>
                error('paritionOptions:badSet','Min_Interior_Angle must be smaller than Max_Interior_Angle.')
            end
            obj.Min_Interior_Angle = value;
        end
        
        function obj = set.Max_Interior_Angle(obj,value)
            validateattributes(value,{'double'},{'scalar','>=',0,'<=',180})
            if value < obj.Min_Interior_Angle %#ok<MCSUP>
                error('paritionOptions:badSet','Max_Interior_Angle must be larger than Min_Interior_Angle.')
            end
            obj.Max_Interior_Angle = value;
        end
        
        function obj = set.Search_Radius(obj,value)
            validateattributes(value,{'double'},{'scalar','integer','nonnegative','real','finite'})
            obj.Search_Radius = value;
        end
        
        function obj = set.Curvature_Smoothing_Size(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','integer','real','finite'})
            obj.Curvature_Smoothing_Size = value;
        end
        
        function obj = set.Use_GPU(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('paritionOptions:badInput','Expected input to be logical.')
            end
            obj.Use_GPU = value;
        end
        
        function obj = set.Use_Parallel(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('paritionOptions:badInput','Expected input to be logical.')
            end
            obj.Use_Parallel = value;
        end
        
        function obj = set.Debug(obj,value)
            if (value ~= 0) && (value ~= 1)
                error('paritionOptions:badInput','Expected input to be logical.')
            end
            obj.Debug = value;
        end
        
        function obj = set.Object_Of_Interest(obj,value)
            if isempty(value)
                obj.Object_Of_Interest = value;
            else
                validateattributes(value,{'double'},{'scalar','integer','nonnegative','real','finite'})
                obj.Object_Of_Interest = value;
            end
        end
        
        function obj = set.Minimum_Hole_Size(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative','real','finite'})
            obj.Minimum_Hole_Size = value;
        end
    end
end





































