function V = MRSI_phantom_to_nifti(phantom, filename)
    %open the file
    f  = fopen(filename, 'w');
    
    %sizeof_hdr: Needs to be 348
    fwrite(f, 348, 'int');
    
    %36 next bytes used for ANALYZE compatability
    fwrite(f, zeros(36,1), 'int8');
    
    %get the dimension sizes
    phantom_size = size(phantom);
    phantom_size = [2, phantom_size];
    
    %make sure dim is of length 8
    dim = zeros(8,1);
    dim(1:length(phantom_size)) = phantom_size;
    %dim[8]: data array dimensions (first element is the number of dimensions)
    fwrite(f, dim, 'short');
    
    %intent header fields not needed. Put zeros there
    fwrite(f, zeros(14,1), 'int8');
    
    %datatype: 512 corresponds to short datatype
    fwrite(f, 512, 'short');
    %bitpix: Nuber of bits per pixel. Needs to be 16 bits for the short data type
    fwrite(f, 16, 'short');
    
    %slice_start: First slice index. For now it is zero
    fwrite(f, zeros(2,1), 'int8');
    
    %pixdim: dimensions of a pixel in x,y,z directions. pixdim[0] is used
    %for the direction of the affine transformation
    pixdim = zeros(8,1);
    delta_x = phantom(2,1).x - phantom(1,1).x;
    delta_y = phantom(1,2).y - phantom(1,1).y;
    pixdim(1:3) = [-1, delta_x, delta_y];
    fwrite(f, pixdim, 'single');
    
    %vox_offset: number of bytes before the first signal position. Minimum
    %352
    vox_offset = 352;
    fwrite(f, vox_offset, 'float32');
    
    %float	scl_slope   Data scaling, slope.
    %float	scl_inter	Data scaling, offset.
    %short	slice_end	Last slice index.
    %char	slice_code 	Slice timing order.
    %all can be zero for now since only in 2D
    fwrite(f, zeros(11,1), 'int8');
    
    %units are in mm. No temporal units.
    %xyzt_units: units for xyz and t. first 3 bits corrispond to xyz, and
    %next 3 corrispont to time. 011 -> mm 000 -> no time.
    xyzt_units = 3;
    fwrite(f, xyzt_units, 'int8');
    
    % float	cal_max	124B	4B	Maximum display intensity.
    % float	cal_min	128B	4B	Minimum display intensity.
    % float	slice_duration	132B	4B	Time for one slice.
    % float	toffset	136B	4B	Time axis shift.
    % int	glmax	140B	4B	Not used; compatibility with analyze.
    % int	glmin	144B	4B	Not used; compatibility with analyze.
    % All can be zero
    fwrite(f, zeros(24,1), 'int8');
    
    %descrip: description of the nifti file. 80 bytes long.
    descrip = 'Simulated phantom';
    descrip(80) = 0;
    fwrite(f, descrip, 'char*1');
    
    %char	aux_file[24]	228B	24B	Auxiliary filename. No aux file
    %needed
    fwrite(f, zeros(24,1), 'int8');
    
    %qform_code: Use the quaternion fields. No need to rotate.
    fwrite(f, 0, 'short');
    
    %sform code: Used to scale and translate image. 1 -> scanner coordinate
    %system
    fwrite(f, 1, 'short');
    
    % float	quatern_b	256B	4B	Quaternion b parameter.
    % float	quatern_c	260B	4B	Quaternion c parameter.
    % float	quatern_d	264B	4B	Quaternion d parameter.
    % float	qoffset_x	268B	4B	Quaternion x shift.
    % float	qoffset_y	272B	4B	Quaternion y shift.
    % float	qoffset_z	276B	4B	Quaternion z shift.
    %quaternion paramaters used for rotation
    fwrite(f, zeros(24,1), 'int8');
    
    %srows: rows of the affine matrix used to scale the points
    srow_x = [1, 0, 0, phantom(1,1).x];
    srow_y = [0, 1, 0, phantom(1,1).y];
    srow_z = [0, 0, 1, 0];
    
    fwrite(f, [srow_x, srow_y, srow_z], 'single');
    
    % char	intent_name[16]	328B	16B	Name or meaning of the data.
    fwrite(f, zeros(16,1), 'int8');
    
    %magic[4]: must end with n+1 or or ni1 for the file to be considered
    %nifti format
    magic_str = zeros(4,1);
    magic_str(1:3) = 'n+1';
    fwrite(f, magic_str, 'char*1');
    fwrite(f, [0,0,0,0], 'int8');
    
    signal = [phantom(:).dI];
    for i = 1:length(signal)
        if signal(i) == 1
            signal(i) = 32767;
        end
    end
    fwrite(f, signal, 'int16');
    
    fclose(f);
    
end
