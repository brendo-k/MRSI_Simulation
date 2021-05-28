function test_suite = MRSI_build_phantom_test
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_points()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([64,64], metabolites);
assertEqual(size(phantom), [64,64])
end

function test_sizeX()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
sizeX = phantom(64).x - phantom(1).x + (phantom(2).x - phantom(1).x);
assertEqual(sizeX, 0.2)
end

function test_sizeY()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
sizeY = phantom(1,64).y - phantom(1,1).y + (phantom(1,2).y - phantom(1,1).y);
assertEqual(sizeY, 0.4)
end

function test_deltaX()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
delta_x = phantom(2, 1).x - phantom(1,1).x;
assertElementsAlmostEqual(delta_x, 0.2/64)
end

function test_deltaY()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
delta_y = phantom(1,2).y - phantom(1,1).y;
assertElementsAlmostEqual(delta_y, 0.4/64)
end

function test_no_metabolites()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([64,64], metabolites, 3);
results = [phantom.met];
assertEqual(results, [])
end

function test_water_sig()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('Could not load H2O.mat, please add FID-A and subdirectories to path')
end
    metabolites = cell(64,64);
    metabolites(30:35, 30:35) = {sysH2O}; 
    phantom = MRSI_build_phantom([0.2,0.2], metabolites, 3);
    met = reshape([phantom.met], [6,6]);
    empty = phantom(setdiff(1:64, 30:35), setdiff(1:64, 30:35));
    empty_met = [empty.met];

    [H, d] = sim_Hamiltonian(sysH2O, 3);
    H.d = d;
    expected(1:6,1:6) = H;

    assertEqual(expected, met)

    %array of empty arrays becomes an empty array
    assertEqual(empty_met, [])
end
