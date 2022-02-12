classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture({'tools', '../MRSI_Trajectory_Simulation/'})})...
        testMRSILoadTrajectory < matlab.unittest.TestCase
    properties
        gMax = 30;
        gamma = 42577000;
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
            trajectory = phaseEncoded('spectralWidth', 2000, 'Fov' ,[200, 200], ...
                'imageSize', [3, 3], 'spectralSize', 512);
            [gradient, time] = MRSI_load_ktrajectory(trajectory, testCase.gMax);
            realGradient = trajectory.k_trajectory(:, 1)/(testCase.gamma*trajectory.dwellTime);
            
            testCase.verifyEqual(gradient(:, 1), realGradient);
            %testing if points after first is zero.
            testCase.verifyEqual(gradient(:, 2:end), zeros(size(gradient(:, 2:end))));
            testCase.verifyEqual(time, repmat(trajectory.dwellTime, size(time)));
        end

    end
end
