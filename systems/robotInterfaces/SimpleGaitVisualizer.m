classdef SimpleGaitVisualizer < Visualizer
% This is a *very* basic visualizer for legged robots. The purpose of this class
% is to make it easy to display footstep plans for robots with arbitrary numbers
% of legs using the same visualizer. No information other than the pose of the 
% robot's body and feet is used to display, so every robot will just look like
% a collection of black sticks. 
  properties
    feet = {};
    safeRegions = struct('A', {}, 'b', {}, 'point', {}, 'normal', {});
    inputFrame;
  end
  
  methods
    function obj = SimpleGaitVisualizer(feet)
      coords = {'body_x', 'body_y', 'body_z', 'body_yaw'};
      for f = feet
        foot = f{1};
        coords = [coords, [foot, '_in_contact'], [foot, '_x'], [foot, '_y'], [foot, '_z'], [foot, '_yaw']];
      end
      inputFrame = CoordinateFrame('SimpleGaitVisualizerInput', 4 + 5 * length(feet), 'x', coords);
      obj = obj@Visualizer(inputFrame);
      obj.inputFrame = inputFrame;
      obj.feet = feet;
    end
    
    function draw(obj, t, x)
      persistent hFig;
      
      if isempty(hFig)
        hFig = figure(25);
        set(hFig,'DoubleBuffer', 'on');
      end
      
      figure(25); 
      cla; hold on;
      
      for j = 1:length(obj.safeRegions)
        reg = obj.safeRegions(j);
        V = iris.thirdParty.polytopes.lcon2vert(reg.A(:,1:2), reg.b);
        V = V';
        V = V(1:2, convhull(V(1,:), V(2,:)));
        
        V(3,:) = 1/(reg.normal(3)) * (reg.normal' * reg.point - (reg.normal(1:2)' * V(1:2,:)));
        
        patch(V(1,:), V(2,:), V(3,:), 'k', 'FaceColor', [0.8,0.8,0.8]);
      end
      
      scale = 0.1;
      p = Point(obj.inputFrame, x);
      plot3(p.body_x, p.body_y, p.body_z, 'ko');
      quiver3(p.body_x, p.body_y, p.body_z, scale * cos(p.body_yaw), scale * sin(p.body_yaw), 0, 'k', 'AutoScale', 'off');
      
      for f = obj.feet
        foot = f{1};
        if p.([foot, '_in_contact'])
          color = 'r';
        else
          color = 'k';
        end
        plot3([p.body_x, p.([foot, '_x'])],...
              [p.body_y, p.([foot, '_y'])],...
              [p.body_z, p.([foot, '_z'])], 'k-', 'LineWidth', 3, 'Color', color);
        quiver3(p.([foot, '_x']), p.([foot, '_y']), p.([foot, '_z']),...
                scale * cos(p.([foot, '_yaw'])), scale * sin(p.([foot, '_yaw'])), 0,...
                'r', 'AutoScale', 'off')
      end
      axis equal
      view(-140,44);
      drawnow();
    end
  end
end
      