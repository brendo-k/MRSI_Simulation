classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture({'tools', '../MRSI_Trajectory_Simulation/'})})...
        testMRSILoadTrajectory < matlab.unittest.TestCase
    properties
        gMax = 30;
        gamma = getGamma("overTwoPi", true);
    end

    methods(Test)
        function testNoGradient(testCase)
            trajectory = phaseEncoded('spectralWidth', 2000, 'Fov' ,[200, 200], ...
                'imageSize', [1, 1], 'spectralSize', 512);
            [gradient, time] = MRSI_load_ktrajectory(trajectory, testCase.gMax);
            testCase.verifyEqual(gradient, zeros(size(gradient)));
            testCase.verifyEqual(time, [0 repmat(trajectory.dwellTime, 1, length(time) - 1)], ...
                'relTol', 1e-10);
        end

        function testPhaseEncoded(testCase)
            %ASSIGN
            trajectory = phaseEncoded('spectralWidth', 2000, 'Fov' ,[200, 200], ...
                'imageSize', [3, 3], 'spectralSize', 512);
            trajectorySize = size(trajectory.k_trajectory);
            timeToFirstPoint = abs(trajectory.k_trajectory(:, 1))/(testCase.gamma*testCase.gMax/10^6);
            maxTime = max(timeToFirstPoint);
            expectedGradient = trajectory.k_trajectory(:, 1)/(testCase.gamma*maxTime);
            expectedGradient = [expectedGradient, zeros(trajectorySize(1), trajectorySize(2) - 1)];
            expectedTime = [repmat(maxTime, trajectorySize(1), 1), repmat(trajectory.dwellTime, trajectorySize(1), trajectorySize(2) - 1)];
            
            
            %APPLY
            [gradient, time] = MRSI_load_ktrajectory(trajectory, testCase.gMax);
            
            %ASSERT
            testCase.verifyEqual(gradient, expectedGradient, 'AbsTol', 1e-8);
            testCase.verifyEqual(time, expectedTime,  'AbsTol', 1e-8);
        end

    end
end
