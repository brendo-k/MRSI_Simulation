clearvars;
load('tmp.mat');

D = myexpm_(Dz,[],[],true,false,1);
a=myexpm_(Dz(:,:,:,[1 11]),[],[],true,false,1);
b=myexpm_(Dz(:,:,1000,end),[],[],true,false,1);
a(:,:,1000,end) % correct
b(:,:,1,1)
D(:,:,1000,end) % wrong
expm(Dz(:,:,1000,end)) % the Matlab function, the correct one