function test_suite = MRSI_evolve_test()
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

%TESTING IF THERE IS NO ROTATION. EXPECTED TO BE NO SIGNAL
function test_no_spin()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O};
phantom = MRSI_build_phantom([0.2,0.2], met, 3);
evolved_phan = MRSI_evolve(phantom, 0);
assertEqual(evolved_phan, phantom);
end

%TESTING FOR ROTATION BUT NO EXCITATION. EXPECTED TO BE NO SIGNAL
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
v_diff = (sysH2O.shifts(1)-4.65)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).d{1});
assertElementsAlmostEqual(sig, before_sig);
end

%TEST ONE FULL ROTATION. EXPECTED TO BE SAME AS STARTING SIGNAL
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
v_diff = (sysH2O.shifts(1)-4.65)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).d{1});
assertElementsAlmostEqual(before_sig, sig, 'relative', 0.0001);
end

%TESTING ONE HALF ROTATION. EXPECTED TO BE NEGATIVE STARTING SINGAL
function test_half_rotation()
try
    load H2O.mat
catch
    moxunit_throw_test_skipped_exception('could not load H2O.mat, please add FID-A')
end
met = {sysH2O, sysH2O; sysH2O, sysH2O};
gamma=42577000;  %[Hz/T]
B0 = 3;
phantom = MRSI_build_phantom([0.2, 0.2], met, B0);
phantom = MRSI_excite(phantom, 90, 'y');
v_diff = (sysH2O.shifts(1) - 4.65)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period/2);
sig1 = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
sig2 = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
sig3 = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
sig4 = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
assertElementsAlmostEqual(-1, sig1, 'relative', 0.001);
assertElementsAlmostEqual(-1, sig2, 'relative', 0.001);
assertElementsAlmostEqual(-1, sig3, 'relative', 0.001);
assertElementsAlmostEqual(-1, sig4, 'relative', 0.001);
end

%TESTING QUAERTER ROTATION. SHOULD BE IMAGINARY SIGNAL
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
v_diff = (sysH2O.shifts(1)-4.65)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period/4);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).d{1});
assertElementsAlmostEqual(sig, -1i*before_sig, 'absolute', 0.001);
end

%TESTING THREE QUARTERS ROTATION. SHOULD BE NEGATIVE IMAGINARY
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
v_diff = (sysH2O.shifts(1)-4.65)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, 3*period/4);
sig = trace((Fx + 1i*Fy)*evolved_phan(1,1).d{1});
before_sig = trace((Fx + 1i*Fy)*phantom(1,1).d{1});
assertElementsAlmostEqual(sig, 1i*before_sig, 'absolute', 0.001);
end

%TESTING ONE EIGHTH ROTATION. SHOULD BE cos(pi/4) + 1i*sin(pi/4)
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
v_diff = (sysH2O.shifts(1)-4.65)*gamma*B0/10^6;
period = 1/v_diff;

Fx = phantom(1,1).met(1).Fx;
Fy = phantom(1,1).met(1).Fy;

evolved_phan = MRSI_evolve(phantom, period/8);
sig = trace((Fx + 1i*Fy)*squeeze(evolved_phan.spins{1}(1,1,:,:)));
assertElementsAlmostEqual(sig, sqrt(1/2)+1i*sqrt(1/2), 'absolute', 0.001);
end