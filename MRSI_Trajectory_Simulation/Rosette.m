% Rosette.m
% Brenden Kadota, McGill University 2019.
% 
% USAGE:
% Rosette(par)
% 
% DESCRIPTION:
% This function takes in a parameter sturucture (par) and outputs the
% corresponding rosette trajectory form the parameters. The formulas to
% create the rosette trajectory comes from Schirda et al. 2009.
%
% INPUTS:
% sw = spectral width
% Fov = Fov array in m, should be 2x1 dimension in mm 
% ImgeSize = Number of points in the x and y dimensions, 2x1 dimension.
%
% OUTPUTS:
% obj = Trajectory obj

function [obj] = Rosette(sw, Fov, imageSize, parameters)
    arguments
        sw (1,1) double = 2000
        Fov (1,2) double = [200,200]
        imageSize (1,3) double = [20, 20, 512]
        parameters.omega1 (1,1) double
        parameters.adc_points (1,1) double
        parameters.N_AngleInts (1,1) double
        parameters.dwellTime (1,1) double
        parameters.kMax (1,1) double
        parameters.spacial_points (1,1) double
        
    end

    %initalize constant
    gamma = 42.577478518e6;    
    
    if(isempty(fieldnames(parameters)))
        %get omega1 and omega2 from sw
        omega1 = sw*pi;
        omega2 = omega1;
        %get kFov
        kFov = imageSize(1)/Fov(1);
        %get KMax
        kMax = kFov/2;
        %calculate gMax
        gMax = kMax * max(omega1, omega2)/(gamma); %T/m

        %dwell time calculation
        dwellTime = 1/(gamma*gMax*Fov(1));
        spectral_points =  (1/sw)/dwellTime;
        dwellTime = (1/sw)/ceil(spectral_points);
        %number of angular interleves (number of turns)
        N_AngleInts = floor((pi*imageSize(1))/(sqrt(1+3*(omega2^2/omega1^2))));
        
        adc_points = (1/sw)*(imageSize(3) - dwellTime);
        %get time vector
        t = 0:dwellTime:adc_points;
    else
        omega1 = parameters.omega1;
        omega2 = parameters.omega1;
        kMax = parameters.kMax;
        dwellTime = parameters.dwellTime;
        N_AngleInts = parameters.N_AngleInts;
        t = 0:dwellTime:dwellTime*(parameters.adc_points-1);
        spectral_points =  (1/sw)/dwellTime;
   
    end
    
    %initalize array size
    traj = complex(zeros(N_AngleInts, length(t)), 0);
    %Rosette trajectory formulat from shirda paper
    traj(1,:) = kMax*sin(omega1*t).*exp(1i*omega2*t);


    %creating rosette trajectory in k space
    for rotation = 1:N_AngleInts
        %get the radians to rotate each trajectory by
        rotationAngle = ((rotation - 1)/N_AngleInts)*2*pi;
        
        %rotate trajecotry in the complex plane using e^(i*theta)
        rotationConverter = exp(1i*rotationAngle); 
        
        %apply to the first trajectory to get the rotated trajectory.
        traj(rotation,:) = rotationConverter*traj(1,:);
    end
    spatial_points = ceil(spectral_points);

    %Create trajectory object
    obj = Trajectory('rosette', traj, imageSize, Fov, dwellTime, sw, t, spatial_points);
end

function mustBeCSV(file_name)
    if isempty(regexp(file_name, '.*\.csv','ignorecase')) 
        eidType = 'mustBeCSV:notCSV';
        msgType = 'Input file must be a csv filename';
        throwAsCaller(MException(eidType,msgType))
    end
end
