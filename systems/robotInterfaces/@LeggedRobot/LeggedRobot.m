classdef LeggedRobot
  % Interface class for legged robots. Currently, this doesn't do much more than
  % occupy the top of the heirarchy containing Biped, Quadruped, etc., but we'd
  % eventually like to move general-purpose multi-leg walking planning into
  % this class.

  properties(Abstract,SetAccess = protected, GetAccess = public)
    foot_names
  end
  properties (Abstract)
    fixed_point_file
  end

  methods
    function obj = LeggedRobot(obj)
      typecheck(obj,{'RigidBodyManipulator','TimeSteppingRigidBodyManipulator'});
    end
  end

  methods(Abstract)
    planFootsteps(obj, goal, params);
  end

  methods
    function xstar = loadFixedPoint(obj)
      load(obj.fixed_point_file, 'xstar');
    end
    function foot_idx = getFootIdx(obj)
      % Get the link index for each foot of the legged robot
      foot_idx = zeros(1,length(obj.foot_names));
      for i = 1:length(obj.foot_names)
          foot_idx(i) = findLinkId(obj,obj.foot_names{i});
      end
    end 
  end
end

