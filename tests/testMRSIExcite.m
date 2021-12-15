classdef testMRSIExcite < matlab.unittest.TestCase
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
            testCase.phantom = MRSI_excite(testCase.phantom, 90, 'y');  
            testCase.resolution = size(metabolites);
            
            testCase.relativeTolerance = single(0.01);
        end
    end
    methods(Test)
        function test_90_flip_x(testCase)
            met = {sysH2O};
            B0 = 3;
            phantom = MRSI_build_phantom([0.2, 0.2], met, B0);

            Fx = phantom(1,1).met(1).Fx;
            Fy = phantom(1,1).met(1).Fy;

            evolved_phan = MRSI_excite(phantom, 90, 'x');
            sig = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
            testCase.verifyEqual(sig, 0 - 1i, 'RelTol', testCase.relativeTolerance);
        end

        function test_no_flip
            try
                load H2O.mat
            catch
                moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
            end
            met = {sysH2O};
            B0 = 3;
            phantom = MRSI_build_phantom([0.2, 0.2], met, B0);

            Fx = phantom(1,1).met(1).Fx;
            Fy = phantom(1,1).met(1).Fy;

            evolved_phan_x = MRSI_excite(phantom, 0, 'x');
            evolved_phan_y = MRSI_excite(phantom, 0, 'y');
            sig_x = trace((Fx + 1i*Fy)*evolved_phan_x(1,1).d{1});
            sig_y = trace((Fx + 1i*Fy)*evolved_phan_y(1,1).d{1});
            assertEqual(sig_x, 0);
            assertEqual(sig_y, 0);
        end

        function test_90_flip_y
            try
                load H2O.mat
            catch
                moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
            end
            met = {sysH2O};
            B0 = 3;
            phantom = MRSI_build_phantom([0.2, 0.2], met, B0);

            Fx = phantom(1,1).met(1).Fx;
            Fy = phantom(1,1).met(1).Fy;

            evolved_phan = MRSI_excite(phantom, 90, 'y');
            sig = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
            assertElementsAlmostEqual(sig, 1 + 1i*0);
        end

        function test_180
            try
                load H2O.mat
            catch
                moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
            end
            met = {sysH2O};
            B0 = 3;
            phantom = MRSI_build_phantom([0.2, 0.2], met, B0);

            Fx = phantom(1,1).met(1).Fx;
            Fy = phantom(1,1).met(1).Fy;

            excite_phan_x = MRSI_excite(phantom, 180, 'x');
            excite_phan_y = MRSI_excite(phantom, 180, 'y');
            sig_x = trace((Fx + 1i*Fy)*squeeze(excite_phan_x.spins{1}(1,1,:,:)));
            sig_y = trace((Fx + 1i*Fy)*squeeze(excite_phan_y.spins{1}(1,1,:,:)));
            assertElementsAlmostEqual(sig_x, 0);
            assertElementsAlmostEqual(sig_y, 0);
        end
    end
end

