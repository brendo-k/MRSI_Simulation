
classdef testMRSIEvolve < matlab.unittest.TestCase
    properties
        metabolite
        phantom
        phantomSize
        resolution
        b0
        relativeTolerance
    end

    methods(TestMethodSetup)
        function createPhantom(testCase)
            h2o = load('H2O.mat');
            testCase.metabolite = h2o.sysH2O;
            testCase.b0 = 3;
            metabolites = {testCase.metabolite};
            testCase.phantomSize = [200, 200];
            testCase.phantom = MRSI_build_phantom(testCase.phantomSize, metabolites, testCase.b0);
            testCase.phantom = MRSI_excite(testCase.phantom, 90, 'y');  
            testCase.resolution = size(metabolites);
            
            testCase.relativeTolerance = single(0.0001);
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

            met = testCase.metabolite;
            B0 = testCase.b0;
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period);
            sig = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            before_sig = testMRSIEvolve.calculateSignal(testCase.phantom.met(1), testCase.phantom.spins{1});
            testCase.verifyEqual(sig, before_sig, 'RelTol', testCase.relativeTolerance);
        end

        %TEST ONE FULL ROTATION. EXPECTED TO BE SAME AS STARTING SIGNAL
        function test_full_rotation(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;     
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period);
            sig = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            before_sig = testMRSIEvolve.calculateSignal(testCase.phantom.met(1), testCase.phantom.spins{1});
            testCase.verifyEqual(sig, before_sig, 'RelTol', testCase.relativeTolerance);
        end

        %TESTING ONE HALF ROTATION. EXPECTED TO BE NEGATIVE STARTING SINGAL
        function test_half_rotation(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;     
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period/2);
            sig = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(-1), 'RelTol', testCase.relativeTolerance);
        end

        %TESTING QUAERTER ROTATION. SHOULD BE IMAGINARY SIGNAL
        function test_quater_rotation(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;     
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period/4);
            sig = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(-1i), 'RelTol', testCase.relativeTolerance);
        end

        %TESTING THREE QUARTERS ROTATION. SHOULD BE NEGATIVE IMAGINARY
        function test_three_quater_rotation(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;     
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, 3*period/4);
            sig = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(1i), 'RelTol', testCase.relativeTolerance);
        end

        %TESTING ONE EIGHTH ROTATION. SHOULD BE cos(pi/4) + 1i*sin(pi/4)
        function test_one_eight_rotation(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;     
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period/8);
            sig = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            testCase.verifyEqual(sig, single(cos(pi/4) - 1i*sin(pi/4)), 'RelTol', testCase.relativeTolerance);
        end

        function test_spin_echo(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period/2);
            evolved_phan = MRSI_excite(evolved_phan, 90, 'x');
            evolved_phan = MRSI_evolve(evolved_phan, period/2);
            spinEchoSignal = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            startingSignal = testMRSIEvolve.calculateSignal(testCase.phantom.met(1), testCase.phantom.spins{1});
            testCase.verifyEqual(spinEchoSignal, startingSignal, 'RelTol', testCase.relativeTolerance);
        end

        function testSpinEchoWithGradient(testCase)
            met = testCase.metabolite;
            B0 = testCase.b0;
            period = testMRSIEvolve.calculatePeriod(met, B0);

            evolved_phan = MRSI_evolve(testCase.phantom, period/2);
            evolved_phan = MRSI_excite(evolved_phan, 90, 'x');
            evolved_phan = MRSI_evolve(evolved_phan, 1e-5);
            evolved_phan = MRSI_evolve(evolved_phan, (period/2)-1e-5);
            spinEchoSignal = testMRSIEvolve.calculateSignal(evolved_phan.met(1), evolved_phan.spins{1});
            startingSignal = testMRSIEvolve.calculateSignal(testCase.phantom.met(1), testCase.phantom.spins{1});
            testCase.verifyEqual(spinEchoSignal, startingSignal, 'RelTol', testCase.relativeTolerance);

        end

    end
    methods(Static)
        function period = calculatePeriod(met, B0)
            gamma = 42577000;
            frequencyDiff = ((met.shifts(1) - 4.65) * gamma * B0) / 10^6;
            period = 1/frequencyDiff;
        end

        function signal = calculateSignal(phantom, spins)
            Fx = phantom.Fx;
            Fy = phantom.Fy;
            signal = trace(pagemtimes((Fx + 1i*Fy), spins));
        end
    end

end
