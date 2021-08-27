function test_suite = MRSI_excite_test
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_90_flip_x
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

evolved_phan = MRSI_excite(phantom, 90, 'x');
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).d{1});
assertElementsAlmostEqual(sig, 0 - 1i);
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
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).d{1});
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
sig_x = trace((Fx + 1i*Fy)*excite_phan_x(1,1).d{1});
sig_y = trace((Fx + 1i*Fy)*excite_phan_y(1,1).d{1});
assertElementsAlmostEqual(sig_x, 0);
assertElementsAlmostEqual(sig_y, 0);
end

