classdef Trajectory
    properties
        name
        k_trajectory
        imageSize 
        FoV
        TR
        pixel_width
        spacial_coordinates
        scanTime
        dwellTime
        t
        sw
        spatialPoints
        spectralLength
    end
    methods
        function obj = Trajectory(name, k_trajectory, imageSize, spectralLength, FoV, dwellTime, ...
            sw, t, spatialPoints, TR)
            arguments
                name (1,:) char
                k_trajectory (:, :) double
                imageSize (1,2) double
                spectralLength (1, 1) double
                FoV (1,2) double
                dwellTime (1,1) double
                sw (1,1) double
                t (1,:) double
                spatialPoints (1,1) double
                TR (1,1) double = 1
            end
            obj.name = name;
            obj.k_trajectory  = k_trajectory;
            obj.pixel_width.y = FoV(1)/imageSize(1);
            obj.pixel_width.x = FoV(2)/imageSize(2);
            obj.imageSize = imageSize;
            obj.spectralLength = spectralLength;
            obj.FoV.y = FoV(2);
            obj.FoV.x = FoV(1);
            obj.t = t;
            obj.dwellTime = dwellTime;
            obj.spacial_coordinates.x = -obj.FoV.x/2+obj.pixel_width.x/2:obj.pixel_width.x:obj.FoV.x/2-obj.pixel_width.x/2;
            obj.spacial_coordinates.y = -obj.FoV.y/2+obj.pixel_width.y/2:obj.pixel_width.y:obj.FoV.y/2-obj.pixel_width.y/2;
            obj.scanTime = size(k_trajectory,1)*(TR + size(k_trajectory,2)*dwellTime);
            obj.TR = TR;
            obj.sw = sw;
            obj.spatialPoints = spatialPoints;
        end
        
        function writeTrajectoryToFile(obj, fileName)
            arguments
                obj (1,1) Trajectory
                fileName (1,:) char 
            end
            
            writematrix(["TR", "K_x", "K_y", "time"], fileName, 'WriteMode', 'overwrite');
            spatial_slice = obj.k_trajectory(:, 1:obj.spatialPoints);
            linear_traj = reshape(spatial_slice.', [], 1);
            linear_traj = [real(linear_traj) imag(linear_traj)];
            tr = 1:size(obj.k_trajectory,1);
            tr = repmat(tr, [obj.spatialPoints, 1]);
            tr = tr(:);
            time = obj.t(1:obj.spatialPoints);
            time = repmat(time, [1, size(obj.k_trajectory,1)]);
            time = time(:);
            svg = horzcat(tr, linear_traj, time);
            
            writematrix(svg, fileName, 'WriteMode', 'append');
            
        end
    end
end