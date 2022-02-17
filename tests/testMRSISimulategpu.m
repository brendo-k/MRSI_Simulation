classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture({'./tools', '../MRSI_Trajectory_Simulation/'})})...
        testMRSISimulategpu < matlab.unittest.TestCase
    properties
        metabolite
        phantom
        phantomSize
        resolution
        b0
        relativeTolerance
        gamma
        scalingFactor
        T2
    end

    methods(TestMethodSetup)
        function createPhantom(testCase)
            h2o = load('H2O.mat');
            testCase.metabolite = h2o.sysH2O;
            testCase.b0 = 3;
            metabolites = cell(3, 3);
            metabolites(:, :) = {testCase.metabolite};
            testCase.phantomSize = [200, 200];
            testCase.phantom = MRSI_build_phantom(testCase.phantomSize, metabolites, testCase.b0);
            testCase.resolution = size(metabolites);
            testCase.gamma = -getGamma("overTwoPi", true);
            testCase.scalingFactor = 2*prod(testCase.resolution);
            testCase.T2 = testCase.phantom.T2(1);
        end
    end
    methods(Test)
        function testNoGradient(testCase)
            %ARRANGE
            trajectory = phaseEncoded('spectralWidth',2000, 'Fov', [200, 200], ...
                'imageSize', [1, 1], 'spectralSize', 512);
            frequency = testMRSISimulategpu.getFrequency(testCase);
            time = trajectory.t;
            points = testMRSISimulategpu.readout(frequency, time, testCase.T2, testCase.scalingFactor);

            %ACT 
            output = MRSI_simulate_gpu(trajectory, testCase.phantom, 'spinEcho', true);

            %ASSERT
            testCase.verifyEqual(output.data, points', 'AbsTol', 1e-3);
        end

        function testYGradient(testCase)
            trajectory = phaseEncoded('spectralWidth', 2000, 'Fov', [201, 201], ...
                'imageSize', [1, 3], 'spectralSize', 512);

            phaseMap = testMRSISimulategpu.getPhaseMapFromKSpace(testCase, trajectory.k_trajectory(:, 1));
            frequency = testMRSISimulategpu.getFrequency(testCase);
            time = trajectory.t;
            points = testMRSISimulategpu.readout(frequency, time, testCase.T2, 2);
            points = repmat(points', [1, prod(testCase.resolution), size(phaseMap, 2)]);
            points = permute(points, [2,3,1]);
            points = points .* phaseMap;
            points = permute(points, [3, 1, 2]);
            points = squeeze(sum(points, 2));
            points = reshape(points, [size(points, 1), 1, size(points, 2)]);

            output = MRSI_simulate_gpu(trajectory, testCase.phantom, 'spinEcho', true);

            testCase.verifyEqual(output.data, points, 'RelTol', 1e-4);
        end

        function testXGradient(testCase)
            %ARRANGE
            trajectory = phaseEncoded("spectralWidth", 2000, 'Fov',[201, 201], ...
                'imageSize', [3, 1], 'spectralSize', 512);

            
            phaseMap = testMRSISimulategpu.getPhaseMapFromKSpace(testCase, trajectory.k_trajectory(:, 1));
            frequency = testMRSISimulategpu.getFrequency(testCase);
            time = trajectory.t;
            points = testMRSISimulategpu.readout(frequency, time, testCase.T2, 2);
            points = repmat(points', [1, size(phaseMap)]);
            points = permute(points, [2, 3, 1]);
            points = points .* phaseMap;
            points = permute(points, [3, 1, 2]);
            points = squeeze(sum(points, 2));
            points = reshape(points, [size(points, 1), size(points, 2)]);

            output = MRSI_simulate_gpu(trajectory, testCase.phantom, 'spinEcho', true);
            testCase.verifyEqual(output.data, points, 'RelTol', 1e-4);
        end
    end
    methods (Static)
        %returns the frequency (in Hertz) of the rotation of the metabolite with
        %no gradients.
        function frequency = getFrequency(testCase)
            frequency = (testCase.metabolite.shifts(1) - 4.65)*testCase.gamma*testCase.b0/1e6;
        end

        %This funciton is used to get phase maps from phase encoding.
        function phaseMap = getPhaseMapFromKSpace(testCase, kTrajectory)
            %load first k-space points
            kStart = kTrajectory;
            %copy the kTrajectory so it is the same resolution as phantom
            kStart = repmat(kStart, [1, testCase.resolution]);
            %permute so y and x dimensions are first
            kStart = permute(kStart, [2, 3, 1]);
            %get kx and ky
            kX = real(kStart);
            kY = imag(kStart);
            [xMesh, yMesh] = meshgrid(testCase.phantom.x, testCase.phantom.y);

            phaseMap = exp(1i*(xMesh.*kX + yMesh.*kY)*2*pi);
            %returning phaseMap of dimensions y*x, kSpacePoints
            phaseMap = reshape(phaseMap, prod(testCase.resolution), length(kTrajectory));
        end

        %compute signal with T2 decay and scaling
        function signal = readout(frequency, time, T2, scaling)
            %create signal based on frequency and time lenght
            signal = exp(frequency*2*pi*time*1i);
            %scale signal based on number of points
            signal = signal .* scaling;
            %apply T2 decay 
            signal = signal .* exp(-time/T2);
        end

    end
end
