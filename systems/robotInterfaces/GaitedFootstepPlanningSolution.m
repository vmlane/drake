classdef GaitedFootstepPlanningSolution
  properties
    robot; 
    t;
    pose;
    full_gait;
    safe_regions;
    region_assignments;
  end

  methods
    function [xtraj, v] = getSimpleGaitTrajectory(obj)
      x = [obj.pose.body];
      feet = fieldnames(obj.full_gait)';
      for f = feet
        foot = f{1};
        x = [x; obj.pose.(foot)];
      end
      xtraj = PPTrajectory(foh(obj.t, x));
      if nargout > 1
        v = SimpleGaitVisualizer(feet);
        xtraj = setOutputFrame(xtraj. v.inputFrame);
      end
    end

    function plan = getBipedFootstepPlan(obj, seed_plan)
      plan = seed_plan;
      nsteps = length(plan.footsteps);
      feet = fieldnames(obj.full_gait)';
      for j = 1:nsteps
        for f = feet
          foot = f{1};
          if plan.footsteps(j).frame_id == obj.robot.foot_frame_id.(foot)
            plan.footsteps(j).pos([1,2,3,6]) = obj.pose.(foot)(:,j);
            plan.footsteps(j).pos([4,5]) = 0;
            plan.region_order(j) = obj.region_assignments(j).(foot);
          end
        end
      end
    end
  end
end
