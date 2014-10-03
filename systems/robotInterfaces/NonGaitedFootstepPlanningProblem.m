classdef NonGaitedFootstepPlanningProblem < FootstepPlanningProblem
  properties
    dt = 0.25;
    g = 9.81; % accel due to gravity (m/s^2)
    % foot_force = [60, 60, 90, 90]; % N = kg m/s^2
    foot_force = 30;
    body_mass = 5; % kg
    nominal_com_height = 0.2; % m
    max_angular_momentum = 3; 
  end

  methods
    function sol = solveYalmip(obj)
      if isempty(obj.start_pose) || isempty(obj.goal_pose)
        error('start_pose and goal_pose should be set before solving');
      end
      checkDependency('yalmip');
      % checkDependency('gurobi');

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
      angular_momentum = struct();
      angular_momentum.body = sdpvar(2, obj.nframes, 'full');
      contact_force = struct('total', 0);
      gait_sum = zeros(1, obj.nframes);
      for j = 1:length(obj.feet)
        foot = obj.feet{j};
        pose.(foot) = sdpvar(4, obj.nframes, 'full');
        gait.(foot) = binvar(1, obj.nframes, 'full');
        gait_sum = gait_sum + gait.(foot);
        contact_force.(foot) = sdpvar(3, obj.nframes, 'full');
        contact_force.total = contact_force.total + contact_force.(foot);
        foot_nominal_pose = [mean([obj.foci.(foot).v], 2); 0];
        angular_momentum.(foot) = cross(foot_nominal_pose, [0;0;1]);
        angular_momentum.(foot) = angular_momentum.(foot)(1:2);
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
        yaw_sector(5,:) == 1,... % Disable yaw
        ];

      % Require the final velocity to be zero (to avoid the solution where
      % we just plummet through the ground forever). 
      % constraints = [constraints, velocity.body(3,end) == 0];

      % Enforce continuity
      constraints = [constraints, velocity.body(:,end) == velocity.body(:,1)];
      constraints = [constraints, acceleration.body(:,end) == acceleration.body(:,1)];
      constraints = [constraints, velocity.body(1,:) >= 1.105];
      constraints = [constraints, pose.body(1:2,1) == 0];
      constraints = [constraints, pose.body(2:3,1) == pose.body(2:3,end)];
      for f = obj.feet
        foot = f{1};
        constraints = [constraints, gait.(foot)(1) == gait.(foot)(end)];
      end

      % Constrain the initial conditions
      % start_fields = fieldnames(obj.start_pose)';
      % for f = start_fields
      %   field = f{1};
      %   constraints = [constraints, pose.(field)(1:2, 1) == obj.start_pose.(field)(1:2),...
      %     velocity.body(:,1) == 0,...
      %   ];
      % end

      for f = obj.feet
        foot = f{1};
      % Feet on the ground at the start and end
        % constraints = [constraints, gait.(foot)([1, end]) == 1];
        % complementarity conditions
        constraints = [constraints,...
          -contact_force.(foot)(3,:) <= contact_force.(foot)(1,:) <= contact_force.(foot)(3,:),...
          -contact_force.(foot)(3,:) <= contact_force.(foot)(2,:) <= contact_force.(foot)(3,:),...
          0 <= contact_force.(foot)(3,:) <= obj.foot_force * gait.(foot),...
          % gait.(foot)(1:2:end) == gait.(foot)(2:2:end),...
          ];
      end

      constraints = [constraints,...
          % angular_momentum.body(:,1) == 0,...
          % angular_momentum.body(:,end) == 0,...
          angular_momentum.body(:,1) == angular_momentum.body(:,end),...
          abs(angular_momentum.body) <= obj.max_angular_momentum,...
          ];

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
          constraints = [constraints, sum(region.(foot), 1) >= gait.(foot)];
        end
      end

      for j = 1:obj.nframes
        if j < obj.nframes
          % Enforce ballistic dynamics
          constraints = [constraints, ...
                         pose.body(1:3,j+1) == pose.body(1:3,j) + velocity.body(1:3,j) * obj.dt,...
                         velocity.body(1:3,j+1) == velocity.body(1:3,j) + acceleration.body(1:3,j) * obj.dt,...
                         acceleration.body(1:3,j) == contact_force.total(1:3,j)/obj.body_mass + [0;0;-obj.g],...
                         ];

          new_momentum = angular_momentum.body(:,j);

          for f = obj.feet
            foot = f{1};
            new_momentum = new_momentum + angular_momentum.(foot) * contact_force.(foot)(3,j);
          end
          constraints = [constraints, ...
            angular_momentum.body(:,j+1) == new_momentum,...
            ];

          % Enforce yaw rate
          constraints = [constraints,...
                         abs(pose.body(4,j+1) - pose.body(4,j)) <= obj.dt * obj.yaw_speed,...
                         ];
        end 

        % Enforce angle approximations
        for s = 1:length(angle_boundaries) - 1
          th0 = angle_boundaries(s);
          th1 = angle_boundaries(s+1);

          th = (th0 + th1)/2;
          ct = cos(th);
          st = sin(th);
          k = tan((th1 - th0)/2) / ((th1 - th0) / 2);
          constraints = [constraints,...
                         implies(yaw_sector(s,j), th0 <= body_yaw(j) <= th1)...
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

        % Foot-specific constraints
        for f = obj.feet
          foot = f{1};

          % Enforce reachability
          constraints = [constraints,...
            obj.body_to_feet_constraints.(foot).A(:,POSE_INDICES) * (pose.(foot)(:,j) - pose.body(:,j)) <= obj.body_to_feet_constraints.(foot).b];
          for k = 1:length(obj.foci.(foot))
            c = obj.foci.(foot)(k);
            constraints = [constraints, ...
              obj.pcone(pose.(foot)(1:2,j) - (pose.body(1:2,j) + [cos_yaw(j), -sin_yaw(j); sin_yaw(j), cos_yaw(j)] * c.v), c.r, 8)];
          end

          % Enforce region membership
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
              implies(gait.(foot)(j), region.(foot)(:,j) == region.(foot)(:,j+1)),...
              obj.pcone((pose.(foot)(1:2,j+1) - pose.(foot)(1:2,j))/obj.dt - velocity.body(1:2,j), obj.swing_speed + MAX_DISTANCE / obj.dt * gait.(foot)(j), 4),...
              ];
          end
        end
      end

      objective = 0;

      % fnames = fieldnames(obj.goal_pose)';
      % for f = fnames
      %   field = f{1};
      %   for k = 1:length(POSE_INDICES)
      %     if ~isnan(obj.goal_pose.(field)(POSE_INDICES(k)))
      %       constraints = [constraints, pose.(field)(k,end) == obj.goal_pose.(field)(POSE_INDICES(k))];
      %     end
      %   end
      % end

      % Penalize contact force
      objective = objective + 0.05 * (sum(sum(contact_force.total.^2,1)) - (obj.g * obj.body_mass)^2 * (obj.nframes-1));

      % Penalize angular momentum
      objective = objective + 0.01 * sum(sum(angular_momentum.body.^2,1));

      % Keep the legs near the body if possible
      for f = obj.feet
        foot = f{1};
        objective = objective + 0.01 * sum(sum(abs(pose.(foot)-pose.body), 1));
        % objective = objective + 0.01 * sum(sum((pose.(foot)-pose.body).^2,1));
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

      figure(6);
      clf
      hold on
      subplot(221)
      plot(t, double(sqrt(sum(contact_force.lf.^2 ,1))));
      title('lf')
      subplot(222)
      plot(t, double(sqrt(sum(contact_force.rf.^2 ,1))));
      title('rf')
      subplot(223)
      plot(t, double(sqrt(sum(contact_force.lh.^2 ,1))));
      title('lh')
      subplot(224)
      plot(t, double(sqrt(sum(contact_force.rh.^2 ,1))));
      title('rh')

      figure(7);
      clf
      hold on
      subplot(221)
      plot(t, double(contact_force.lf(3,:)));
      title('lf')
      subplot(222)
      plot(t, double(contact_force.rf(3,:)));
      title('rf')
      subplot(223)
      plot(t, double(contact_force.lh(3,:)));
      title('lh')
      subplot(224)
      plot(t, double(contact_force.rh(3,:)));
      title('rh')

      figure(8);
      clf
      hold on
      subplot(221)
      plot(t(1:end-1), sqrt(sum(double((pose.lf(1:2,2:end) - pose.lf(1:2,1:end-1))/obj.dt - velocity.body(1:2,1:end-1)).^2,1)));
      title('lf')
      subplot(222)
      plot(t(1:end-1), sqrt(sum(double((pose.rf(1:2,2:end) - pose.rf(1:2,1:end-1))/obj.dt - velocity.body(1:2,1:end-1)).^2,1)));
      title('rf')
      subplot(223)
      plot(t(1:end-1), sqrt(sum(double((pose.lh(1:2,2:end) - pose.lh(1:2,1:end-1))/obj.dt - velocity.body(1:2,1:end-1)).^2,1)));
      title('lh')
      subplot(224)
      plot(t(1:end-1), sqrt(sum(double((pose.rh(1:2,2:end) - pose.rh(1:2,1:end-1))/obj.dt - velocity.body(1:2,1:end-1)).^2,1)));
      title('rh')

    end
  end
  methods(Static)
    function constraint = pcone(v, r, num_pieces)
      % Polynomial approximation of a conic constraint norm(v) <= r
      A = zeros(num_pieces, 2);
      b = repmat(r, num_pieces, 1);
      ths = linspace(0, 2*pi, num_pieces);
      for j = 1:num_pieces
        th = ths(j);
        c = cos(th);
        s = sin(th);
        R = [c, -s; s, c];
        A(j,:) = (R * [1;0])';
      end
      constraint = A * v <= b;
    end
        

  end
end
