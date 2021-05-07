function tests = MRSI_build_phantom_test
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
%import path fixtures class
import matlab.unittest.fixtures.PathFixture
%apply path fixtures to temp add the simulation package
testCase.applyFixture(PathFixture(['../']));
end


function test_points(testCase)
phantom = MRSI_build_phantom([64,64], [0.2,0.2], [20,20], [45,45]);
verifyEqual(testCase, size(phantom), [64,64])
end

function test_sizeX(testCase)
phantom = MRSI_build_phantom([64,64], [0.2,0.4], [20,20], [45,45]);
sizeX = phantom(64).x - phantom(1).x + (phantom(2).x - phantom(1).x);
verifyEqual(testCase, sizeX, 0.2)
end

function test_sizeY(testCase)
phantom = MRSI_build_phantom([64,64], [0.2,0.4], [20,20], [45,45]);
sizeX = phantom(1,64).y - phantom(1,1).y + (phantom(1,2).y - phantom(1,1).y);
verifyEqual(testCase, sizeX, 0.4)
end

function test_sig(testCase)
phantom = MRSI_build_phantom([64,64], [0.2,0.4], [20,20], [45,45]);
equal = [phantom(20:45, 20:45).dI] == 1;
verifyEqual(testCase, sum(equal), 676)
end