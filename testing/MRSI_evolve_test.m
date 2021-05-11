function tests = MRSI_evolve_test
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
%import path fixtures class
import matlab.unittest.fixtures.PathFixture
%apply path fixtures to temp add the simulation package
testCase.applyFixture(PathFixture(['../']));
end

function setup(testCase)
testCase.TestData.phantom = MRSI_build_phantom([32,32], [0.2,0.2], [10,10], [20,20]);
testCase.TestData.Ix = complex(0.5*[0 1;1 0]);

end

function test_no_spin(testCase)
evolved_phan = MRSI_evolve(testCase.TestData.phantom, 0, 3, testCase.TestData.Ix);
verifyEqual(testCase, evolved_phan, testCase.TestData.phantom);
end
