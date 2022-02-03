classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture({'./tools', '../MRSI_Trajectory_Simulation/'})})...
        testMRSISimulate < matlab.unittest.TestCase
    properties
        metabolite
        phantom
        phantomSize
        resolution
        b0
        relativeTolerance
        gamma
        scalingFactor
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
            testCase.relativeTolerance = single(0.0001);
            testCase.gamma = -42577000;
            testCase.scalingFactor = 2*prod(testCase.resolution);
        end
    end
    methods(Test)
        function testNoGradient(testCase)
            trajectory = phaseEncoded('spectralWidth',2000, 'Fov', [200, 200], ...
                'imageSize', [1, 1], 'spectralSize', 512);
            output = MRSI_simulate(trajectory, testCase.phantom, 30, 3, 'spinEcho', true);
            frequency = (testCase.metabolite.shifts(1) - 4.65)*testCase.gamma*testCase.b0/1e6;
            time = 0:trajectory.dwellTime:trajectory.dwellTime*(trajectory.imageSize(3) - 1);

            points = exp(frequency*2*pi*time*1i)*testCase.scalingFactor .* exp(-time/testCase.phantom.T2);

            testCase.verifyEqual(output.data, points', 'RelTol', 0.0001);
        end

        function testYGradient(testCase)
            trajectory = phaseEncoded('spectralWidth', 2000, 'Fov', [201, 201], ...
                'imageSize', [3, 1], 'spectralSize', 512);
            output = MRSI_simulate(trajectory, testCase.phantom, 30, 3, 'spinEcho', true);

            kStart = trajectory.k_trajectory(:, 1);
            kStart = repmat(kStart, [1, testCase.resolution]);
            kStart = permute(kStart, [2, 3, 1]);
            kX = real(kStart);
            kY = imag(kStart);
            [xMesh, yMesh] = meshgrid(testCase.phantom.x, testCase.phantom.y);

            phaseMap = exp(1i*(xMesh.*kX + yMesh.*kY)*2*pi);
            phaseMap = reshape(phaseMap, prod(testCase.resolution), []);
            
            
            frequency = (testCase.metabolite.shifts(1) - 4.65)*testCase.gamma*testCase.b0/1e6;
            time = 0:trajectory.dwellTime:trajectory.dwellTime*(trajectory.imageSize(3) - 1);
            points = exp(frequency*2*pi*time*1i)*2 .* exp(-time/testCase.phantom.T2);
            points = repmat(points', [1, prod(testCase.resolution), size(phaseMap, 2)]);
            points = permute(points, [2,3,1]);
            points = points .* phaseMap;
            points = permute(points, [3, 1, 2]);
            points = squeeze(sum(points, 2));
            points = reshape(points, [size(points, 1), 1, size(points, 2)]);
            testCase.verifyEqual(output.data, points, 'RelTol', 0.0001);
        end

        function testXGradient(testCase)
            trajectory = phaseEncoded("spectralWidth", 2000, 'Fov',[201, 201], ...
                'imageSize', [1, 3], 'spectralSize');
            output = MRSI_simulate(trajectory, testCase.phantom, 30, 3, 'spinEcho', true);

            kStart = trajectory.k_trajectory(:, 1);
            kStart = repmat(kStart, [1, testCase.resolution]);
            kStart = permute(kStart, [2, 3, 1]);
            kX = real(kStart);
            kY = imag(kStart);
            [xMesh, yMesh] = meshgrid(testCase.phantom.x, testCase.phantom.y);

            phaseMap = exp(1i*(xMesh.*kX + yMesh.*kY)*2*pi);
            phaseMap = reshape(phaseMap, prod(testCase.resolution), []);
            
            
            frequency = (testCase.metabolite.shifts(1) - 4.65)*testCase.gamma*testCase.b0/1e6;
            time = 0:trajectory.dwellTime:trajectory.dwellTime*(trajectory.imageSize(3) - 1);
            points = exp(frequency*2*pi*time*1i)*2 .* exp(-time/testCase.phantom.T2);
            points = repmat(points', [1, prod(testCase.resolution), size(phaseMap, 2)]);
            points = permute(points, [2,3,1]);
            points = points .* phaseMap;
            points = permute(points, [3, 1, 2]);
            points = squeeze(sum(points, 2));
            points = reshape(points, [size(points, 1), size(points, 2)]);
            testCase.verifyEqual(output.data, points, 'RelTol', 0.0001);
        end
    end
end
