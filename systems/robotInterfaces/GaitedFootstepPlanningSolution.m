classdef GaitedFootstepPlanningSolution
% This class is a helper object to hold the solution to a GaitedFootstepPlanningProblem. 
% It allows the user to retrieve a FootstepPlan object compatible with the DRC stack or
% to visualize the plan using a simple general-purpose Matlab visualizer. 
  properties
    t;
    pose;
    full_gait;
    safe_regions;
    region_assignments;
  end

  methods
    function [xtraj, v] = getSimpleGaitTrajectory(obj)
      % Return a trajectory of [body_pose; foot1_pose; foot2_pose; etc...] and,
      % optionally, a visualizer to display it with.
      x = [obj.pose.body];
      feet = fieldnames(obj.full_gait)';
      for f = feet
        foot = f{1};
        x = [x; obj.pose.(foot)];
      end
      xtraj = PPTrajectory(foh(obj.t, x));
      if nargout > 1
        v = SimpleGaitVisualizer(feet);
        for j = 1:length(obj.safe_regions)
          normal = obj.safe_regions(j).eq.A(1,1:3)';
          point = repmat(obj.safe_regions(j).eq.b / sum(normal), 3, 1);
          assert(normal' * point == obj.safe_regions(j).eq.b);
          v.safeRegions(j) = struct('A', obj.safe_regions(j).ineq.A, 'b', obj.safe_regions(j).ineq.b, 'point', point, 'normal', normal);
        end
        xtraj = setOutputFrame(xtraj, v.inputFrame);
      end
    end

    function plan = getBipedFootstepPlan(obj, biped, seed_plan)
      plan = seed_plan;
      nsteps = length(plan.footsteps);
      feet = fieldnames(obj.full_gait)';
      for j = 1:nsteps
        for f = feet
          foot = f{1};
          if plan.footsteps(j).frame_id == biped.foot_frame_id.(foot)
            plan.footsteps(j).pos([1,2,3,6]) = obj.pose.(foot)(:,j);
            plan.footsteps(j).pos([4,5]) = 0;
            plan.region_order(j) = obj.region_assignments(j).(foot);
          end
        end
      end
    end
  end
end
