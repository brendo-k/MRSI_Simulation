function test_suite = MRSI_evolve_test()
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

%function setupOnce(testCase)
%import path fixtures class
%import matlab.unittest.fixtures.PathFixture
%apply path fixtures to temp add the simulation package
%testCase.applyFixture(PathFixture(['../']));
%end
%
%
%function test_no_spin(testCase)
%phantom = MRSI_build_phantom([0.2,0.2], 
%evolved_phan = MRSI_evolve(testCase.TestData.phantom, 0, 3, testCase.TestData.Ix);
%verifyEqual(testCase, evolved_phan, testCase.TestData.phantom);
%end
%
%function setup()

%end
