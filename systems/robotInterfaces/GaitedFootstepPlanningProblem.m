classdef GaitedFootstepPlanningProblem
  properties
    robot;
    feet = {'left', 'right'};
    foci = struct('right', struct('v', {[0; 0], [0; -0.25]},...
                                             'r', {0.19, 0.16}),...
                                 'left', struct('v', {[0; 0], [0; 0.25]},...
                                             'r', {0.19, 0.16}));
    body_to_feet_constraints; % currently must be linear
    gait = struct('right', {0, 1}, 'left', {1, 0});
    nframes = 10; % number of frames of the gait to plan. nframes/length(gait) gives the number of complete gait cycles
    swing_speed = 1; % m/s
    body_speed = 0.25; % m/s
    yaw_speed = 1; % rad/s
    start_pose = [];
    goal_pose = [];
    safe_regions = struct('ineq', {}, 'eq', {});
  end

  methods
    function obj = GaitedFootstepPlanningProblem(robot)
      obj.robot = robot;
    end

    function obj = addIRISRegions(obj, iris_regions)
      for j = 1:length(iris_regions)
        A = iris_regions(j).A;
        A = [A(:,1:2), zeros(size(A, 1), 3), A(:,3)];
        b = iris_regions(j).b;
        point = iris_regions(j).point;
        normal = iris_regions(j).normal;
        n = size(A, 1);
        m = size(A, 2);
        inframe = drakeFunction.frames.realCoordinateSpace(m);
        outframe = drakeFunction.frames.realCoordinateSpace(n);
        obj.safe_regions(end+1).ineq = drakeFunction.Affine(inframe, outframe, A, b);

        A = [point', zeros(1, 3)];
        b = point' * normal;
        outframe = drakeFunction.frames.realCoordinateSpace(1);
        obj.safe_regions(end).eq = drakeFunction.Affine(inframe, outframe, A, b);
      end
    end


    function sol = solveYalmip(obj)
      if isempty(start_pose) || isempty(goal_pose)
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
      for f = obj.feet
        foot = f{1};
        pos.(foot) = sdpvar(4, obj.nframes, 'full');
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

      min_yaw = obj.start_pose.(body)(6) - pi;
      max_yaw = obj.start_pose.(body)(6) + pi;

      % Replicate the gait (so we can hand it back as an output)
      full_gait = repmat(obj.gait, 1, ceil(length(obj.gait) / obj.nframes));
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
        constraints = [constraints, pose.(field)(:, 1) == obj.start_pose.(field)(POSE_INDICES)];
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
            if j > length(obj.gait) || any(~full_gait.(foot)(1:j)) % If the foot has been allowed to move yet
              for r = 1:nregions
                Ar_ineq = obj.safe_regions(r).ineq.A(:,POSE_INDICES);
                br_ineq = obj.safe_regions(r).ineq.b;
                Ar_eq = obj.safe_regions(r).eq.A(:,POSE_INDICES);
                br_eq = obj.safe_regions(r).eq.b;
                constraints = [constraints,...
                  implies(region.(foot)(r,j), Ar_ineq * pose.(foot)(:,j) <= br_ineq),...
                  implies(region.(foot)(r,j), Ar_eq * pose.(foot) == br_eq),...
                  ];
              end
            end
          end
          
          if j < obj.nframes
            % enforce gait and swing speed
            if full_gait(j).(foot)
              constraints = [constraints, pose.(foot)(:,j) == pose.(foot)(:,j+1),...
                             region.(foot)(:,j) == region.(foot)(:,j+1)];
            else
              constraints = [constraints, ...
                  abs(pose.(foot)(1:2,j+1) - pose.(foot)(1:2,j)) <= dt(j) * SWING_SPEED,...
      %           cone(pos.(foot)(1:2,j+1) - pos.(foot)(1:2,j), (dt(j)) * SWING_SPEED),...
                  ];
            end
          end
        end
        if j < obj.nframes
          constraints = [constraints,...
                         abs(pose.body(4,j+1) - pose.body(4,j)) <= dt(j) * ROTATION_RATE,...
                        cone(pose.body(1:2,j+1) - pose.body(1:2,j), dt(j) * BODY_SPEED),...
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
          region_assignments.(foot)(j) = r;
        end
      end

      dt = double(dt);
      t = [0, cumsum(dt(1:end-1))];

      sol = GaitedFootstepPlanningSolution();
      sol.robot = robot;
      sol.t = t;
      sol.pose = pose;
      sol.full_gait = full_gait;
      sol.safe_regions = obj.safe_regions;
      sol.region_assignments = region_assignments;
    end
  end
end



