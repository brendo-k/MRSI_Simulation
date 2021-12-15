%% Main function to generate tests
classdef testMRSIBuildPhantom < matlab.unittest.TestCase
    properties
        phantom
        resolution
        phantomSize
        metaboliteRegions
        metabolite 
    end

    methods(TestMethodSetup)
        function createPhantom(testCase)
            metabolites = cell(128,64);
            h2o = load('H2O.mat');
            testCase.metabolite = h2o.sysH2O;
            
            metabolites(1:32, 1:32) = {h2o.sysH2O};
            testCase.phantomSize = [400, 200];
            testCase.phantom = MRSI_build_phantom(testCase.phantomSize, metabolites);
            testCase.resolution = size(metabolites);
        end
    end
    methods(Test)

        function test_points(testCase)
            testCase.verifyEqual(size(testCase.phantom.spins{1}, [3,4]), [128,64])
        end

        function test_sizeX(testCase)
            sizeX = testCase.phantom.x(end) - testCase.phantom.x(1) + (testCase.phantom.x(2) - testCase.phantom.x(1));
            testCase.verifyEqual(sizeX, testCase.phantomSize(2))
        end

        function test_sizeY(testCase)

            sizeY = testCase.phantom.y(end) - testCase.phantom.y(1) + (testCase.phantom.y(2) - testCase.phantom.y(1));
            testCase.verifyEqual(sizeY, testCase.phantomSize(1))
        end

        function test_deltaX(testCase)
            delta_x = testCase.phantom.x(2) - testCase.phantom.x(1);
            testCase.verifyEqual(delta_x, testCase.phantomSize(2)/testCase.resolution(2))
        end

        function test_deltaY(testCase)
            delta_y = testCase.phantom.y(2) - testCase.phantom.y(1);
            testCase.verifyEqual(delta_y, testCase.phantomSize(1)/testCase.resolution(1))
        end

        %make sure density matrix and Hamiltonians are properly copied over
        function test_water_sig(testCase)
            met = testCase.phantom.met(1);
            den = testCase.phantom.spins{1};

            h2o = testCase.metabolite;
            h2o.shifts = h2o.shifts - 4.65;
            [H, d] = sim_Hamiltonian(h2o, 3);
            d = repmat(single(d{1}), [1, 1, 32, 32]);
            testCase.verifyEqual(H, met);
            testCase.verifyEqual(den(:, :, 1:32, 1:32), d);
            
        end
    end
end
