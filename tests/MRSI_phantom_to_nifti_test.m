classdef testMRSIPhantomtoNifti

    methods(TestMethodSetup)
        function createPhantom(testCase)
            load H2O.mat
            met = cell(64,64);
            met(1,1) = {sysH2O};
            phantom = MRSI_build_phantom([0.2,0.2], met, 3);
            tempFile = 'test.nii';
            MRSI_phantom_to_nifti(phantom, tempFile);
            file = fopen('test.nii', 'r');
            bytes = fread(file);
        end

    end

    function teardown(file)
        fclose(file);
        delete test.nii
    end

    function test_header_size()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(int8(bytes(1:4))', 'int32');
        assertEqual(header, int32(348));
        teardown(file)
    end

    function test_dims()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(int8(bytes(41:56))', 'int16');
        assertEqual(header, int16([2,64,64,0,0,0,0,0]));
        teardown(file)
    end

    function test_data_type()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(int8(bytes(71:72))', 'int16');
        assertEqual(header, int16(512));
        teardown(file);
    end

    function test_bitpix()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(int8(bytes(73:74))', 'int16');
        assertEqual(header, int16(16));
        teardown(file);
    end

    function test_pixdim()
        [tempFile, file, bytes] = setup();
        met = cell(64,64);
        phantom = MRSI_build_phantom([0.2,0.2], met, 3);
        deltaX = phantom.x(2) - phantom.x(1);
        deltaY = phantom.y(2) - phantom.y(1);
        bytes = bytes;
        header = typecast(uint8(bytes(77:108))', 'single');
        assertEqual(header, single([-1, deltaX*1000, deltaY*1000, 0, 0, 0, 0, 0]));
        teardown(file);
    end

    function test_vox_offset()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(uint8(bytes(109:112))', 'single');
        assertEqual(header, single(352));
        teardown(file);
    end

    function test_xyzt_units()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = uint8(bytes(124));
        assertEqual(header, uint8(2));
        teardown(file);
    end

    function test_description()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = nonzeros(uint8(bytes(149:228)));
        assertEqual(char(header'), 'Simulated phantom')
        teardown(file);
    end

    function test_qform_code()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(uint8(bytes(253:254))', 'uint16');
        assertEqual(header, uint16(0))
        teardown(file);
    end

    function test_sform_code()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = typecast(uint8(bytes(255:256))', 'uint16');
        assertEqual(header, uint16(1))
        teardown(file);
    end

    function test_srow()
        [tempFile, file, bytes] = setup();
        met = cell(64,64);
        phantom = MRSI_build_phantom([0.2,0.2], met, 3);
        bytes = bytes;
        header = typecast(uint8(bytes(281:328))', 'single');
        deltaX = phantom.y(2) - phantom.y(1);
        deltaY = phantom.x(2) - phantom.x(1);
        row_x = header(1:4);
        row_y = header(5:8);
        row_z = header(9:12);
        assertEqual(row_x, single([deltaX, 0, 0, phantom.x(1)]))
        assertEqual(row_y, single([0, deltaY, 0, phantom.y(1)]))
        assertEqual(row_z, single([0, 0, 1, 0]))
        teardown(file);
    end

    function test_magic_str()
        [tempFile, file, bytes] = setup();
        bytes = bytes;
        header = char(nonzeros(uint8(bytes(345:348)))');
        assertEqual(header, 'n+1')
        teardown(file);
    end
end
