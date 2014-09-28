classdef GaitedFootstepPlanningProblem < FootstepPlanningProblem
  % This class represents a general category of fixed-gait footstep planning problems. A
  % fixed-gait planning problem involves optimizing the poses of a number of "frames", each
  % of which consists of a pose for the robot's body and all of its feet. The gait is a binary
  % structure, which indicates whether a given foot is in contact at a given frame. Feet that
  % are in contact may not move between frames. Additionally, we use the same convex safe
  % region technique used in, e.g. footstepMIQP to assign foot poses to convex regions.
  % 
  % This class is demonstrated in bipedGaitedMISOCP.m and testQuadrupedPlanner.m
  properties
    gait = struct('right', {0, 1}, 'left', {1, 0});
  end

  methods
    function sol = solveYalmip(obj)
      if isempty(obj.start_pose) || isempty(obj.goal_pose)
        error('start_pose and goal_pose should be set before solving');
      end
      checkDependency('yalmip');
      checkDependency('gurobi');

      % A general upper limit on the distance covered in a single plan (in meters);
      MAX_DISTANCE = 30;

      % Which indices of [x, y, z, roll, pitch, yaw] do we actually use
      POSE_INDICES = [1,2,3,6];

      nregions = length(obj.safe_regions);
      pose = struct();
      pose.body = sdpvar(4, obj.nframes, 'full');
      region = struct();
      for f = obj.feet
        foot = f{1};
        pose.(foot) = sdpvar(4, obj.nframes, 'full');
        if nregions > 0
          region.(foot) = binvar(nregions, obj.nframes, 'full');
        end
      end

      body_yaw = pose.body(4,:);
      cos_yaw = sdpvar(1, obj.nframes, 'full');
      sin_yaw = sdpvar(1, obj.nframes, 'full');

      dt = sdpvar(1, obj.nframes, 'full');

      angle_boundaries = (-pi-pi/8):(pi/4):(pi-pi/8);
      yaw_sector = binvar(length(angle_boundaries) - 1, obj.nframes, 'full');

      min_yaw = obj.start_pose.body(6) - pi;
      max_yaw = obj.start_pose.body(6) + pi;

      % Replicate the gait (so we can hand it back as an output)
      full_gait = repmat(obj.gait, 1, ceil(obj.nframes / length(obj.gait)));
      full_gait = full_gait(1:obj.nframes);


      % Set up the basic constraints, including some general ranges on body pose (to make the mixed-integer formulation easier for yalmip)
      constraints = [...
        pose.body(:,1) == obj.start_pose.body(POSE_INDICES),...
        min_yaw <= body_yaw <= max_yaw,...
        obj.start_pose.body(1) - MAX_DISTANCE <= pose.body(1,:) <= obj.start_pose.body(1) + MAX_DISTANCE,...
        obj.start_pose.body(2) - MAX_DISTANCE <= pose.body(2,:) <= obj.start_pose.body(2) + MAX_DISTANCE,...
        obj.start_pose.body(3) - MAX_DISTANCE <= pose.body(3,:) <= obj.start_pose.body(3) + MAX_DISTANCE,...
        sum(yaw_sector, 1) == 1,...
        -1 <= sin_yaw <= 1,...
        -1 <= cos_yaw <= 1,...
        0 <= dt <= MAX_DISTANCE / obj.swing_speed,...
        ];

      % Constrain the initial conditions
      start_fields = fieldnames(obj.start_pose)';
      for f = start_fields
        field = f{1};
        constraints = [constraints, pose.(field)(1:2, 1) == obj.start_pose.(field)(1:2)];
      end

      % Set up general bounds on foot poses and foot region assignments
      for f = obj.feet
        foot = f{1};
        constraints = [constraints,...
          obj.start_pose.body(1) - MAX_DISTANCE <= pose.(foot)(1,:) <= obj.start_pose.body(1) + MAX_DISTANCE,...
          obj.start_pose.body(2) - MAX_DISTANCE <= pose.(foot)(2,:) <= obj.start_pose.body(2) + MAX_DISTANCE,...
          obj.start_pose.body(3) - MAX_DISTANCE <= pose.(foot)(3,:) <= obj.start_pose.body(3) + MAX_DISTANCE,...
          min_yaw <= pose.(foot)(4,:) <= max_yaw,...
          ];
        if nregions > 0
          constraints = [constraints, sum(region.(foot), 1) == 1];
        end
      end

      for j = 1:obj.nframes
        for f = obj.feet
          foot = f{1};
          
          % Enforce reachability
          constraints = [constraints,...
            obj.body_to_feet_constraints.(foot).A(:,POSE_INDICES) * (pose.(foot)(:,j) - pose.body(:,j)) <= obj.body_to_feet_constraints.(foot).b];

          for k = 1:length(obj.foci.(foot))
            c = obj.foci.(foot)(k);
            constraints = [constraints, cone(pose.(foot)(1:2,j) - (pose.body(1:2,j) + [cos_yaw(j), -sin_yaw(j); sin_yaw(j), cos_yaw(j)] * c.v), c.r)];
          end
          
          if nregions > 0
            if j > length(obj.gait) || any(~[full_gait(1:j).(foot)]) % If the foot has been allowed to move yet
              for r = 1:nregions
                Ar_ineq = obj.safe_regions(r).ineq.A(:,POSE_INDICES);
                br_ineq = obj.safe_regions(r).ineq.b;
                Ar_eq = obj.safe_regions(r).eq.A(:,POSE_INDICES);
                br_eq = obj.safe_regions(r).eq.b;
                constraints = [constraints,...
                  implies(region.(foot)(r,j), Ar_ineq * pose.(foot)(:,j) <= br_ineq),...
                  implies(region.(foot)(r,j), Ar_eq * pose.(foot)(:,j) == br_eq),...
                  ];
              end
            end
          end
          
          if j < obj.nframes
            % enforce gait and swing speed
            if full_gait(j).(foot)
              constraints = [constraints, pose.(foot)(:,j) == pose.(foot)(:,j+1)];
              if nregions > 0
                constraints = [constraints, region.(foot)(:,j) == region.(foot)(:,j+1)];
              end
            else
              constraints = [constraints, ...
                  abs(pose.(foot)(1:2,j+1) - pose.(foot)(1:2,j)) <= dt(j) * obj.swing_speed,...
      %           cone(pos.(foot)(1:2,j+1) - pos.(foot)(1:2,j), (dt(j)) * SWING_SPEED),...
                  ];
            end
          end
        end
        if j < obj.nframes
          constraints = [constraints,...
                         abs(pose.body(4,j+1) - pose.body(4,j)) <= dt(j) * obj.yaw_speed,...
                        cone(pose.body(1:2,j+1) - pose.body(1:2,j), dt(j) * obj.body_speed),...
                         ];
        end
        
        for s = 1:length(angle_boundaries) - 1
          th0 = angle_boundaries(s);
          th1 = angle_boundaries(s+1);

          th = (th0 + th1)/2;
          ct = cos(th);
          st = sin(th);
          k = tan((th1 - th0)/2) / ((th1 - th0) / 2);
          constraints = [constraints,...
                         implies(yaw_sector(s,j), th0 <= body_yaw(j) <= th1)...
                         % implies(yaw_sector(s,j), [ct; st]' * [cos_yaw(j); sin_yaw(j)] == 1),...
                         implies(yaw_sector(s,j), [cos_yaw(j); sin_yaw(j)] == [ct; st] + (body_yaw(j) - th) * k * [-st; ct])];
          if j < obj.nframes
            constraints = [constraints,...
                           sum(yaw_sector(max(s-1, 1):min(s+1, size(yaw_sector,1)), j+1)) >= yaw_sector(s, j)];
          end
          if j < obj.nframes - 1
            constraints = [constraints,...
                           sum(yaw_sector(max(s-1, 1):min(s+1, size(yaw_sector,1)), j+2)) >= yaw_sector(s, j)];
          end
        end
      end

      objective = sum(dt);
      w_goal = diag([100,100,100,100]);

      fnames = fieldnames(obj.goal_pose)';
      for f = fnames
        field = f{1};
        objective = objective + (pose.(field)(:,end) - obj.goal_pose.(field)(POSE_INDICES))' * w_goal * (pose.(field)(:,end) - obj.goal_pose.(field)(POSE_INDICES));
      end
      for j = 1:obj.nframes
        for f = obj.feet
          foot = f{1};
          objective = objective + (pose.(foot)(:,j) - pose.body(:,j))' * (pose.(foot)(:,j) - pose.body(:,j));
          if j < obj.nframes
            objective = objective + (pose.(foot)(:,j+1) - pose.(foot)(:,j))' * (pose.(foot)(:,j+1) - pose.(foot)(:,j));
          end
        end
      end

      optimize(constraints, objective, sdpsettings('solver', 'gurobi'));

      % Extract the result
      fnames = fieldnames(pose)';
      for f = fnames
        field = f{1};
        pose.(field) = double(pose.(field));
      end
      region_assignments = struct();

      for f = obj.feet
        foot = f{1};
        region.(foot) = logical(round(double(region.(foot))));
        for j = 1:obj.nframes
          r = find(region.(foot)(:,j));
          assert(length(r) == 1);
          region_assignments(j).(foot) = r;
        end
      end

      dt = double(dt);
      t = [0, cumsum(dt(1:end-1))];

      sol = GaitedFootstepPlanningSolution();
      sol.t = t;
      sol.pose = pose;
      sol.full_gait = full_gait;
      sol.safe_regions = obj.safe_regions;
      sol.region_assignments = region_assignments;
    end
  end
end


