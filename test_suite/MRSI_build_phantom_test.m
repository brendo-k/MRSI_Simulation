function test_suite = MRSI_build_phantom_test
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_points()
metabolites = cell(128,64);
phantom = MRSI_build_phantom([64,64], metabolites);
assertEqual(size(phantom.spins{1}, [1,2]), [128,64])
end

function test_sizeX()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([200,400], metabolites, 3);
sizeX = phantom.x(64) - phantom.x(1) + (phantom.x(2) - phantom.x(1));
assertEqual(sizeX, 400)
end

function test_sizeY()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
sizeY = phantom.y(64) - phantom.y(1) + (phantom.y(2) - phantom.y(1));
assertEqual(sizeY, 0.2)
end

function test_deltaX()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
delta_x = phantom.x(2) - phantom.x(1);
assertElementsAlmostEqual(delta_x, 0.4/64)
end

function test_deltaY()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([0.2,0.4], metabolites, 3);
delta_y = phantom.y(2) - phantom.y(1);
assertElementsAlmostEqual(delta_y, 0.2/64)
end

function test_no_metabolites()
metabolites = cell(64,64);
phantom = MRSI_build_phantom([64,64], metabolites, 3);
results = [phantom.met];
assertEqual(results, [])
end

%make sure density matrix and Hamiltonians are properly copied over
function test_water_sig()
try
    load H2O.mat sysH2O
catch
    moxunit_throw_test_skipped_exception('Could not load H2O.mat, please add FID-A and subdirectories to path')
end
    metabolites = cell(64,64);
    metabolites(30:35, 30:35) = {sysH2O}; 
    phantom = MRSI_build_phantom([0.2,0.2], metabolites, 3);
    met = phantom.met(1);
    den = phantom.spins{1};
 

    sysH2O.shifts = sysH2O.shifts - 4.65;
    [H, d] = sim_Hamiltonian(sysH2O, 3);
   
    assertEqual(H, met)
    expected2 = den(30:35,30:35,:,:);
    d = repmat(d{1}, [1,1,6,6]);
    d = permute(d, [3,4,1,2]);
    d = single(d);
    assertEqual(d, expected2);

end
