function [poses] = footstepGaitedMISOCP(robot, feet_description, start_pose, goal_pose, safe_regions, nframes)

checkDependency('yalmip');
checkDependency('gurobi');

MAX_DISTANCE = 30;

nregions = length(seed_plan.safe_regions);

pos = struct();
pos.body = sdpvar(4, nframes, 'full');
cos_yaw = sdpvar(1, nframes, 'full');
sin_yaw = sdpvar(1, nframes, 'full');
body_yaw = pos.body(4,:);
dt = sdpvar(1, nframes, 'full');
region = struct();
for f = walk_info.feet
  foot = f{1};
  pos.(foot) = sdpvar(4, nframes, 'full');
  region.(foot) = binvar(nregions, nframes, 'full');
end
angle_boundaries = (-pi-pi/8):(pi/4):(pi-pi/8);
yaw_sector = binvar(length(angle_boundaries) - 1, nframes, 'full');

seed_steps = [seed_plan.footsteps.pos];

min_yaw = seed_steps(6,1) - pi;
max_yaw = seed_steps(6,1) + pi;

constraints = [...
               % TODO: generalize this
               pos.right(:,1) == seed_steps([1,2,3,6], 1),...
               pos.left(:,1) == seed_steps([1,2,3,6], 2),...
               body_yaw(1) == seed_steps(6,1),...
               min_yaw <= body_yaw <= max_yaw,...
               seed_steps(1,1) - MAX_DISTANCE <= pos.body(1,:) <= seed_steps(1,1) + MAX_DISTANCE,...
               seed_steps(2,1) - MAX_DISTANCE <= pos.body(2,:) <= seed_steps(2,1) + MAX_DISTANCE,...
               seed_steps(3,1) - MAX_DISTANCE <= pos.body(3,:) <= seed_steps(3,1) + MAX_DISTANCE,...
               sum(yaw_sector, 1) == 1,...
               -1 <= sin_yaw <= 1,...
               -1 <= cos_yaw <= 1,...
               0 <= dt <= MAX_DISTANCE / SWING_SPEED,...
%                yaw_sector(5,:) == 1,...
               ];
             
for f = walk_info.feet
  foot = f{1};
  constraints = [constraints,...
                 sum(region.(foot), 1) == 1,...
                 seed_steps(1,1) - MAX_DISTANCE <= pos.(foot)(1,:) <= seed_steps(1,1) + MAX_DISTANCE,...
                 seed_steps(2,1) - MAX_DISTANCE <= pos.(foot)(2,:) <= seed_steps(2,1) + MAX_DISTANCE,...
                 seed_steps(3,1) - MAX_DISTANCE <= pos.(foot)(3,:) <= seed_steps(3,1) + MAX_DISTANCE,...
                 min_yaw <= pos.(foot)(4,:) <= max_yaw,...
                 ];
end

for j = 1:nframes
  for f = walk_info.feet
    foot = f{1};
    
    % Enforce reachability
    constraints = [constraints, walk_info.lin_con.(foot).A(:,[1,2,3,6]) * (pos.(foot)(:,j) - pos.body(:,j)) <= walk_info.lin_con.(foot).b];
    for k = 1:length(walk_info.foci.(foot))
      c = walk_info.foci.(foot)(k);
      constraints = [constraints, cone(pos.(foot)(1:2,j) - (pos.body(1:2,j) + [cos_yaw(j), -sin_yaw(j); sin_yaw(j), cos_yaw(j)] * c.v), c.r)];
    end
    
    % TODO: generalize number of fixed steps
    if j > 2
      % Enforce region membership
      for r = 1:nregions
        Ar = [seed_plan.safe_regions(r).A(:,1:2), zeros(size(seed_plan.safe_regions(r).A, 1), 1), seed_plan.safe_regions(r).A(:,3)];
        constraints = [constraints, implies(region.(foot)(r,j), Ar * pos.(foot)(:,j) <= seed_plan.safe_regions(r).b),...
                       implies(region.(foot)(r,j), seed_plan.safe_regions(r).normal' * pos.(foot)(1:3,j) == seed_plan.safe_regions(r).normal' * seed_plan.safe_regions(r).point),...
                       ];
      end                   
    end
    
    if j < nframes
      % enforce gait and swing speed
      if walk_info.gait(mod(j-1, length(walk_info.gait)) + 1).(foot)
        constraints = [constraints, pos.(foot)(:,j) == pos.(foot)(:,j+1),...
                       region.(foot)(:,j) == region.(foot)(:,j+1)];
      else
        constraints = [constraints, ...
            abs(pos.(foot)(1:2,j+1) - pos.(foot)(1:2,j)) <= dt(j) * SWING_SPEED,...
%           cone(pos.(foot)(1:2,j+1) - pos.(foot)(1:2,j), (dt(j)) * SWING_SPEED),...
            ];
      end
    end
  end
  if j < nframes
    constraints = [constraints,...
                   % abs(pos.body(1:2,j+1) - pos.body(1:2,j)) <= dt(j) * BODY_SPEED,...
                   abs(pos.body(4,j+1) - pos.body(4,j)) <= dt(j) * ROTATION_RATE,...
                  cone(pos.body(1:2,j+1) - pos.body(1:2,j), dt(j) * BODY_SPEED),...
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
    if j < nframes
      constraints = [constraints,...
                     sum(yaw_sector(max(s-1, 1):min(s+1, size(yaw_sector,1)), j+1)) >= yaw_sector(s, j)];
    end
    if j < nframes - 1
      constraints = [constraints,...
                     sum(yaw_sector(max(s-1, 1):min(s+1, size(yaw_sector,1)), j+2)) >= yaw_sector(s, j)];
    end
  end
end

objective = sum(dt);
w_goal = diag([100,100,100,100]);

objective = objective + (pos.right(:,end) - goal_pos.right([1,2,3,6]))' * w_goal * (pos.right(:,end) - goal_pos.right([1,2,3,6]));
objective = objective + (pos.left(:,end) - goal_pos.left([1,2,3,6]))' * w_goal * (pos.left(:,end) - goal_pos.left([1,2,3,6]));
for j = 1:nframes
  for f = walk_info.feet
    foot = f{1};
    objective = objective + (pos.(foot)(:,j) - pos.body(:,j))' * (pos.(foot)(:,j) - pos.body(:,j));
    if j < nframes
      objective = objective + (pos.(foot)(:,j+1) - pos.(foot)(:,j))' * (pos.(foot)(:,j+1) - pos.(foot)(:,j));
    end
  end
end

solvesdp(constraints, objective, sdpsettings('solver', 'gurobi'));

pos.body = double(pos.body);
for f = walk_info.feet
  foot = f{1};
  pos.(foot) = double(pos.(foot));
  region.(foot) = logical(round(double(region.(foot))));
end
dt = double(dt);
t = [0, cumsum(dt(1:end-1))];

x = [pos.body];
for f = walk_info.feet
  foot = f{1};
  x = [x; pos.(foot)];
end
xtraj = PPTrajectory(foh(t, x));
v = SimpleGaitVisualizer(walk_info.feet);
% v.safeRegions = seed_plan.safe_regions;
xtraj = setOutputFrame(xtraj, v.inputFrame);

v.playback(xtraj, struct('slider', true));

plan = seed_plan;

for j = 1:nsteps
  for f = walk_info.feet
    foot = f{1};
    if plan.footsteps(j).frame_id == robot.foot_frame_id.(foot)
      plan.footsteps(j).pos([1,2,3,6]) = pos.(foot)(:,j);
      plan.region_order(j) = find(region.(foot)(:,j));
    end
  end
  plan.footsteps(j).pos([4,5]) = 0;
end


end