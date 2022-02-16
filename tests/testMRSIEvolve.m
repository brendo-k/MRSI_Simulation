classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture({'tools', '../MRSI_Trajectory_Simulation/'})})...
        testMRSIEvolve < matlab.unittest.TestCase
    properties
        metabolite
        phantom
        phantomSize
        resolution
        b0
        relativeTolerance
        period
    end

    methods(TestMethodSetup)
        function createPhantom(testCase)
            h2o = load('H2O.mat');
            testCase.metabolite = h2o.sysH2O;
            testCase.b0 = 3;
            metabolites = {testCase.metabolite};
            testCase.phantomSize = [200, 200];
            testCase.phantom = MRSI_build_phantom(testCase.phantomSize, metabolites, testCase.b0);
            testCase.resolution = size(metabolites);
            testCase.period = testMRSIEvolve.calculatePeriod(testCase.metabolite, testCase.b0);
            
            testCase.relativeTolerance = single(0.000001);
        end
    end
    methods(Test)
        %TESTING IF THERE IS NO ROTATION. EXPECTED TO BE NO SIGNAL
        function test_no_spin(testCase)
            evolved_phan = MRSI_evolve(testCase.phantom, 0);
            testCase.verifyEqual(evolved_phan.d{1}, testCase.phantom.d{1});
        end

        %TESTING FOR ROTATION BUT NO EXCITATION. EXPECTED TO BE NO SIGNAL
        function test_rotation_no_excite(testCase)
            evolved_phan = MRSI_evolve(testCase.phantom, testCase.period);
            sig = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            before_sig = calculateSignal(testCase.phantom.met(1), testCase.phantom.spins{1});
            testCase.verifyEqual(sig, before_sig, 'RelTol', testCase.relativeTolerance);
        end

        %TEST ONE FULL ROTATION. EXPECTED TO BE SAME AS STARTING SIGNAL
        function test_full_rotation(testCase)
            excitedPhantom = MRSI_excite(testCase.phantom, 90, 'y');
            evolved_phan = MRSI_evolve(excitedPhantom, testCase.period);
            sig = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            before_sig = calculateSignal(excitedPhantom.met(1), excitedPhantom.spins{1});

            testCase.verifyEqual(sig, before_sig, 'RelTol', testCase.relativeTolerance);
        end

        %TESTING ONE HALF ROTATION. EXPECTED TO BE NEGATIVE STARTING SINGAL
        function test_half_rotation(testCase)
            excitedPhantom = MRSI_excite(testCase.phantom, 90, 'y');
            evolved_phan = MRSI_evolve(excitedPhantom, testCase.period/2);
            sig = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(-1), 'RelTol', testCase.relativeTolerance);
        end

        %TESTING QUAERTER ROTATION. SHOULD BE IMAGINARY SIGNAL
        function test_quater_rotation(testCase)
            excitedPhantom = MRSI_excite(testCase.phantom, 90, 'y');
            evolved_phan = MRSI_evolve(excitedPhantom, testCase.period/4);
            sig = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(-1i), 'RelTol', testCase.relativeTolerance);
        end

        %TESTING THREE QUARTERS ROTATION. SHOULD BE NEGATIVE IMAGINARY
        function test_three_quater_rotation(testCase)
            excitedPhantom = MRSI_excite(testCase.phantom, 90, 'y');
            evolved_phan = MRSI_evolve(excitedPhantom, 3*testCase.period/4);
            sig = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(1i));
        end

        %TESTING ONE EIGHTH ROTATION. SHOULD BE cos(pi/4) + 1i*sin(pi/4)
        function test_one_eight_rotation(testCase)
            excitedPhantom = MRSI_excite(testCase.phantom, 90, 'y');
            evolved_phan = MRSI_evolve(excitedPhantom, testCase.period/8);
            sig = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(cos(pi/4) - 1i*sin(pi/4)));
        end

        function test_spin_echo(testCase)
            excitedPhantom = MRSI_excite(testCase.phantom, 90, 'y');
            evolved_phan = MRSI_evolve(excitedPhantom, testCase.period/4);
            evolved_phan = MRSI_excite(evolved_phan, 180, 'x');
            evolved_phan = MRSI_evolve(evolved_phan, testCase.period/4);
            spinEchoSignal = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            startingSignal = calculateSignal(excitedPhantom.met(1), excitedPhantom.spins{1});
            testCase.verifyEqual(spinEchoSignal, startingSignal);
        end

        function testSpinEchoWithGradient(testCase)
            evolved_phan = MRSI_evolve(testCase.phantom, testCase.period/4);
            evolved_phan = MRSI_excite(evolved_phan, 180, 'x');
            evolved_phan = MRSI_evolve(evolved_phan, 50);
            evolved_phan = MRSI_evolve(evolved_phan, (testCase.period/4) - 50);
            spinEchoSignal = calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            startingSignal = calculateSignal(testCase.phantom.met(1), testCase.phantom.spins{1});
            testCase.verifyEqual(spinEchoSignal, startingSignal, 'RelTol', testCase.relativeTolerance);

        end

    end
    methods(Static)
        function period = calculatePeriod(met, B0)
            gamma = getGamma('overTwoPi', true) / 1000;
            frequency = ((met.shifts(1) - 4.65) * gamma * B0) / 10^6;
            period = 1/frequency;
        end
    end

end
