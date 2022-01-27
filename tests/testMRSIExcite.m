classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture({'tools', '../MRSI_Trajectory_Simulation/'})})...
        testMRSIExcite < matlab.unittest.TestCase
    properties
        metabolite
        phantom
        phantomSize
        resolution
        b0
        relativeTolerance
    end
    methods(TestMethodSetup)
        function createPhantom(testCase)
            h2o = load('H2O.mat');
            testCase.metabolite = h2o.sysH2O;
            testCase.b0 = 3;
            metabolites = {testCase.metabolite};
            testCase.phantomSize = [200, 200];
            testCase.phantom = MRSI_build_phantom(testCase.phantomSize, metabolites, testCase.b0);
            testCase.resolution = size(metabolites);
            
            testCase.relativeTolerance = single(0.01);
        end
    end
    methods(Test)
        function test_90_flip_x(testCase)
            evolved_phan = MRSI_excite(testCase.phantom, 90, 'x');
            sig = calculateSignal(testCase.phantom.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(0 - 1i), 'RelTol', testCase.relativeTolerance);
        end

        function test_no_flip(testCase)
            evolved_phan_x = MRSI_excite(testCase.phantom, 0, 'x');
            evolved_phan_y = MRSI_excite(testCase.phantom, 0, 'y');
            sig_x = calculateSignal(testCase.phantom.met(1), evolved_phan_x.spins{1});
            sig_y = calculateSignal(testCase.phantom.met(1), evolved_phan_y.spins{1});
            testCase.verifyEqual(sig_x, single(0));
            testCase.verifyEqual(sig_y, single(0));
        end

        function test_90_flip_y(testCase)
            evolved_phan = MRSI_excite(testCase.phantom, 90, 'y');
            sig = calculateSignal(testCase.phantom.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(1 + 1i*0), 'RelTol', testCase.relativeTolerance);
        end

        function test_180(testCase)
            excite_phan_x = MRSI_excite(testCase.phantom, 180, 'x');
            excite_phan_y = MRSI_excite(testCase.phantom, 180, 'y');
            sig_x = calculateSignal(testCase.phantom.met(1), excite_phan_x.spins{1});
            sig_y = calculateSignal(testCase.phantom.met(1), excite_phan_y.spins{1});
            testCase.verifyEqual(sig_x, single(0), 'AbsTol', testCase.relativeTolerance);
            testCase.verifyEqual(sig_y, single(0), 'AbsTol', testCase.relativeTolerance);
        end
    end
end

