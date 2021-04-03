load('H2O.mat');

b0 = 3;
data = [];
%Phantom for now
for i = 1:64
    for j = 1:64
        if(i>=20 && i<=40 && j>=20 && j<=40)
            data(i,j).metabolites = sysH2O;
        else
            data(i,j).metabolites = sysH2O;
        end
    end
end

for i = 1:size(data,1)
    for j = 1:size(data,2)
        data(i,j).b0 = b0;
    end
end

for i = 1:size(data,1)
    for j = 1:size(data,2)
        data(i,j).b0 = b0;
    end
end

phantom.b0 = b0;
phantom.fovX = .2; %in m
phantom.fovZ = .2; %in m
phantom.deltaX = .2/size(data, 1);
phantom.deltaZ = .2/size(data, 2);

%Coordinates start at the center of the voxel
x = -phantom.fovX/2+phantom.deltaX/2:phantom.deltaX:phantom.fovX/2-phantom.deltaX/2;
z = -phantom.fovZ/2 +phantom.deltaZ/2:phantom.deltaZ:phantom.fovZ/2-phantom.deltaZ/2;
[Z,X] = meshgrid(x,z);

%Z dimension 2, X dimension 1
for i = 1:size(data,1)
    for j = 1:size(data,2)
        data(i,j).coordX = X(i,j);
        data(i,j).coordZ = Z(i,j);
    end
end
phantom.data = data;
phantom.data = add_grad(phantom.data, 5, 2);


for i = 1:size(phantom.data, 1)
    [h,d1]=sim_Hamiltonian(phantom.data(1,i).metabolites, phantom.data(1,i).b0);
    d1=sim_excite(d1,h,'x');                            
    [out,dout]=sim_readout(d1,h,64,5000,3,90);
    out_arr(i) = out;
    outd_arr(i) = dout; 
    H(i) = h;
    d(i) = d1;
end

figure
hold on
for i = 1:size(out_arr,2)
    plot(out_arr(i).ppm, real(out_arr(i).specs))
end
hold off

final = zeros(size(out_arr(1).fids));
for i = 1:size(phantom.data, 1)
    final = final + out_arr(i).fids;
end
figure
plot(real(final));

figure
plot(abs(fftshift(fft(final))))

%sim_onepulse(size(phantom.data, 1), 51092, 2.8, 1, phantom.data(1,i).metabolites);

%adding gradient to a single dimension (gradient in mt/m)
function obj = add_grad(obj, gradient, dimension)
    if(length(size(obj)) < dimension)
        error('Pick a valid dimension')
    end
    if (dimension == 1)
        for i = 1:size(obj,1)
            local_grad = gradient*obj(i,1).coordX*10^-3 ;
            for j = 1:size(obj,2)
                obj(i,j).b0 = obj(i,j).b0 + local_grad;
            end
        end

    elseif (dimension == 2)
        for j = 1:size(obj,2)
            local_grad = gradient*obj(1,j).coordZ*10^-3;
            for i = 1:size(obj,1)
                obj(i,j).b0 = obj(i,j).b0 + local_grad;
            end
        end
    end 
end



            
       
