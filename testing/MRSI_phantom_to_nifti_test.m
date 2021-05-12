function tests = MRSI_phantom_to_nifti_test
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
import matlab.unittest.fixtures.PathFixture
%apply path fixtures to temp add the simulation package
testCase.applyFixture(PathFixture(['../']));

testCase.TestData.phantom = MRSI_build_phantom([32,32], [0.2,0.2], [10,10], [20,20]);
end

function setup(testCase)
testCase.TestData.tempFile = 'test.nii';
MRSI_phantom_to_nifti(testCase.TestData.phantom, testCase.TestData.tempFile);
testCase.TestData.file = fopen('test.nii', 'r');
testCase.TestData.bytes = fread(testCase.TestData.file);
end

function teardown(testCase)
fclose(testCase.TestData.file);
delete test.nii
end

function test_header_size(testCase)
bytes = testCase.TestData.bytes;
header = typecast(int8(bytes(1:4))', 'int32');
verifyEqual(testCase, header, int32(348));
end

function test_dims(testCase)
bytes = testCase.TestData.bytes;
header = typecast(int8(bytes(41:56))', 'int16');
verifyEqual(testCase, header, int16([2,32,32,0,0,0,0,0]));
end

function test_data_type(testCase)
bytes = testCase.TestData.bytes;
header = typecast(int8(bytes(71:72))', 'int16');
verifyEqual(testCase, header, int16(512));
end

function test_bitpix(testCase)
bytes = testCase.TestData.bytes;
header = typecast(int8(bytes(73:74))', 'int16');
verifyEqual(testCase, header, int16(16));
end

function test_pixdim(testCase)
phantom = testCase.TestData.phantom;
deltaX = phantom(2,1).x - phantom(1,1).x;
deltaY = phantom(1,2).y - phantom(1,1).y;


bytes = testCase.TestData.bytes;
header = typecast(uint8(bytes(77:108))', 'single');
verifyEqual(testCase, header, single([-1, deltaX*1000, deltaY*1000, 0, 0, 0, 0, 0]));
end

function test_vox_offset(testCase)
bytes = testCase.TestData.bytes;
header = typecast(uint8(bytes(109:112))', 'single');
verifyEqual(testCase, header, single(352));
end

function test_xyzt_units(testCase)
bytes = testCase.TestData.bytes;
header = uint8(bytes(124));
verifyEqual(testCase, header, uint8(2));
end

function test_description(testCase)
bytes = testCase.TestData.bytes;
header = nonzeros(uint8(bytes(149:228)));
verifyEqual(testCase, char(header'), 'Simulated phantom')
end

function test_qform_code(testCase)
bytes = testCase.TestData.bytes;
header = typecast(uint8(bytes(253:254))', 'uint16');
verifyEqual(testCase, header, uint16(0))
end

function test_sform_code(testCase)
bytes = testCase.TestData.bytes;
header = typecast(uint8(bytes(255:256))', 'uint16');
verifyEqual(testCase, header, uint16(1))
end

function test_srow(testCase)
phantom = testCase.TestData.phantom;
bytes = testCase.TestData.bytes;
header = typecast(uint8(bytes(281:328))', 'single');
deltaX = phantom(2,1).x - phantom(1,1).x;
deltaY = phantom(1,2).y - phantom(1,1).y;

row_x = header(1:4);
row_y = header(5:8);
row_z = header(9:12);
verifyEqual(testCase, row_x, single([deltaX*1000, 0, 0, phantom(1,1).x*1000]))
verifyEqual(testCase, row_y, single([0, deltaY*1000, 0, phantom(1,1).y*1000]))
verifyEqual(testCase, row_z, single([0, 0, 1, 0]))
end

function test_magic_str(testCase)
bytes = testCase.TestData.bytes;
header = char(nonzeros(uint8(bytes(345:348)))');
verifyEqual(testCase, header, 'n+1')
end