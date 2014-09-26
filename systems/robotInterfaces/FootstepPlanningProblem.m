classdef FootstepPlanningProblem
  properties
    feet = {'left', 'right'};
    foci = struct('right', struct('v', {[0; 0], [0; -0.25]},...
                                             'r', {0.19, 0.16}),...
                                 'left', struct('v', {[0; 0], [0; 0.25]},...
                                             'r', {0.19, 0.16}));
    body_to_feet_constraints = struct('right', struct('A', [0,0,1,0,0,0;
                                                      0,0,-1,0,0,0;
                                                      0,0,0,0,0,1;
                                                      0,0,0,0,0,-1],...
                                                'b', [0.1;
                                                      0.1;
                                                      0.01;
                                                      pi/8]),...
                                     'left', struct('A', [0,0,1,0,0,0;
                                                      0,0,-1,0,0,0;
                                                      0,0,0,0,0,1;
                                                      0,0,0,0,0,-1],...
                                                'b', [0.1;
                                                      0.1;
                                                      pi/8;
                                                      0.01])); % currently must be linear
    nframes = 10; % number of frames of the gait to plan. nframes/length(gait) gives the number of complete gait cycles
    swing_speed = 1; % m/s
    body_speed = 0.25; % m/s
    yaw_speed = 1; % rad/s
    start_pose = [];
    goal_pose = [];
    safe_regions = struct('ineq', {}, 'eq', {});
  end

  methods
    function obj = addIRISRegions(obj, iris_regions)
      for j = 1:length(iris_regions)
        A = iris_regions(j).A;
        A = [A(:,1:2), zeros(size(A, 1), 3), A(:,3)];
        b = iris_regions(j).b;
        point = iris_regions(j).point;
        normal = iris_regions(j).normal;
        obj.safe_regions(end+1).ineq = struct('A', A, 'b', b);

        A = [normal', zeros(1, 3)];
        b = normal' * point;
        obj.safe_regions(end).eq = struct('A', A, 'b', b);
      end
    end
  end
end


