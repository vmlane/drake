classdef NonGaitedFootstepPlanningProblem < FootstepPlanningProblem
  properties
    dt = 0.25;
    g = 9.81; % accel due to gravity (m/s^2)
    foot_force = 60; % N = kg m/s^2
    body_mass = 5; % kg
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
      velocity.body = sdpvar(3, obj.nframes, 'full');
      acceleration.body = sdpvar(3, obj.nframes, 'full');

      region = struct();
      gait = struct();
      gait_sum = zeros(1, obj.nframes);
      for f = obj.feet
        foot = f{1};
        pose.(foot) = sdpvar(4, obj.nframes, 'full');
        gait.(foot) = binvar(1, obj.nframes, 'full');
        gait_sum = gait_sum + gait.(foot);
        if nregions > 0
          region.(foot) = binvar(nregions, obj.nframes, 'full');
        end
      end

      body_yaw = pose.body(4,:);
      cos_yaw = sdpvar(1, obj.nframes, 'full');
      sin_yaw = sdpvar(1, obj.nframes, 'full');

      angle_boundaries = (-pi-pi/8):(pi/4):(pi-pi/8);
      yaw_sector = binvar(length(angle_boundaries) - 1, obj.nframes, 'full');

      min_yaw = obj.start_pose.body(6) - pi;
      max_yaw = obj.start_pose.body(6) + pi;

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
        ];

      % Require the final velocity to be zero in the z direction (to avoid the solution where
      % we just plummet through the ground forever). 
      constraints = [constraints, velocity.body(3,end) == 0];

      % Constrain the initial conditions
      start_fields = fieldnames(obj.start_pose)';
      for f = start_fields
        field = f{1};
        constraints = [constraints, pose.(field)(1:2, 1) == obj.start_pose.(field)(1:2),...
          gait_sum(1) == length(obj.feet),...
          ];
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
          constraints = [constraints, sum(region.(foot), 1) == gait.(foot)];
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

          if j < obj.nframes
            % Enforce ballistic dynamics
            constraints = [constraints, pose.body(1:3,j+1) - pose.body(1:3,j) == velocity.body(1:3,j),...
                           velocity.body(1:3,j+1) == velocity.body(1:3,j) + acceleration.body(1:3,j) * obj.dt + [0; 0; -obj.g * obj.dt],...
                           cone(velocity.body(1:3,j), obj.body_speed),...
                           ];
            % Enforce bounded impulse from each foot
            constraints = [constraints, cone(acceleration.body(1:3,j), gait_sum(j) * obj.foot_force / obj.body_mass)];
          end 
          
          if nregions > 0
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
          
          if j < obj.nframes
            % enforce gait and swing speed
            constraints = [constraints,...
              implies(gait.(foot)(j), pose.(foot)(:,j) == pose.(foot)(:,j+1)),...
              cone(velocity.body(1:3,j) - (pose.(foot)(1:3,j+1) - pose.(foot)(1:3,j)), obj.dt * obj.swing_speed)];

            % if nregions > 0
            %   constraints = [constraints, implies(gait.(foot)(j), region.(foot)(:,j) == region.(foot)(:,j+1))];
            % end
          end
        end
        if j < obj.nframes
          constraints = [constraints,...
                         abs(pose.body(4,j+1) - pose.body(4,j)) <= obj.dt * obj.yaw_speed,...
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

      objective = 0;
      w_goal = diag([100,100,100,100]);

      fnames = fieldnames(obj.goal_pose)';
      for f = fnames
        field = f{1};
        objective = objective + (pose.(field)(:,end) - obj.goal_pose.(field)(POSE_INDICES))' * w_goal * (pose.(field)(:,end) - obj.goal_pose.(field)(POSE_INDICES));
      end
      for j = 1:obj.nframes
        for f = obj.feet
          foot = f{1};
          % objective = objective + (pose.(foot)(1:2,j) - pose.body(1:2,j))' * (pose.(foot)(1:2,j) - pose.body(1:2,j));
          if j < obj.nframes
            objective = objective + norm((pose.(foot)(1:3,j+1) - pose.(foot)(1:3,j)) - velocity.body(1:3,j));
            objective = objective + norm(acceleration.body(1:3,j) - [0;0;9.8]);
            % objective = objective + (pose.(foot)(:,j+1) - pose.(foot)(:,j))' * (pose.(foot)(:,j+1) - pose.(foot)(:,j));
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
        gait.(foot) = logical(round(double(gait.(foot))));
        for j = 1:obj.nframes
          r = find(region.(foot)(:,j));
          assert(length(r) <= 1);
          region_assignments(j).(foot) = r;
        end
      end

      t = 0:obj.dt:((obj.nframes - 1) * obj.dt);

      sol = GaitedFootstepPlanningSolution();
      sol.t = t;
      sol.pose = pose;
      sol.full_gait = gait;
      sol.safe_regions = obj.safe_regions;
      sol.region_assignments = region_assignments;

      figure(5)
      clf
      hold on
      subplot(311)
      plot(t, double(acceleration.body(3,:)));
      subplot(312)
      plot(t, double(velocity.body(3,:)));
      subplot(313)
      plot(t, double(pose.body(3,:)));
    end
  end
end
