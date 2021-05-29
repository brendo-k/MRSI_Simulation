function test_suite = MRSI_evolve_test()
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_no_spin()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
phantom = MRSI_build_phantom([0.2,0.2], met, 3);
evolved_phan = MRSI_evolve(phantom, 0, 3);
assertEqual(evolved_phan, phantom);
end

function test_rotation_no_excite()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
v_diff = sysH2O.shifts(1)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period, B0);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).met(1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).met(1).d{1});
assertElementsAlmostEqual(sig, before_sig);
end

function test_full_rotation()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
phantom = MRSI_excite(phantom, 90, 'y');
v_diff = sysH2O.shifts(1)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period, B0);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).met(1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).met(1).d{1});
assertElementsAlmostEqual(sig, before_sig, 'absolute', 0.001);
end

function test_half_rotation()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
phantom = MRSI_excite(phantom, 90, 'y');
v_diff = sysH2O.shifts(1)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period/2, B0);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).met(1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).met(1).d{1});
assertElementsAlmostEqual(sig, -before_sig, 'absolute', 0.001);
end

function test_quater_rotation()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
phantom = MRSI_excite(phantom, 90, 'y');
v_diff = sysH2O.shifts(1)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period/4, B0);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).met(1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).met(1).d{1});
assertElementsAlmostEqual(sig, -1i*before_sig, 'absolute', 0.001);
end

function test_three_quater_rotation()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
phantom = MRSI_excite(phantom, 90, 'y');
v_diff = sysH2O.shifts(1)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, 3*period/4, B0);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).met(1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).met(1).d{1});
assertElementsAlmostEqual(sig, 1i*before_sig, 'absolute', 0.001);
end

function test_one_eight_rotation()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
phantom = MRSI_excite(phantom, 90, 'y');
v_diff = sysH2O.shifts(1)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period/8, B0);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).met(1).d{1});
assertElementsAlmostEqual(sig, sqrt(1/2)-1i*sqrt(1/2), 'absolute', 0.001);
end