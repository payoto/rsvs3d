classdef RSVS3D_log < RSVS3D_DerivativeStudy
    
  
    
    methods(Access=protected)
        
        function BuildStructure(obj,varargin)
            % Function which builds the data structure at the core of the
            % object.
            % This function should be overriden in subclasses to specialise
            % how the structure is built
            [obj.resout] = ParseConvergenceRSVS3D(obj.derivDirectory, ...
                obj.derivFileName);
        end
    end
    methods
        function obj = RSVS3D_log(derivDirectory)
            obj.derivDirectory = derivDirectory;
            obj.derivFileName = 'convergence_';
        end
    end
end