% my_script.m
clear variables
tic
result = zeros(64, 64,3000, 'like', 1i);
toc
%%assign real values to a real vector
disp('assign real values to a real vector');
idx = get_diag_index(64, 3000);
tic
exp = complex(ones([64,3000]),1);
toc
result(1) = 0 + 1i;
tic
for i = length(idx):-1:1
result(idx(i)) = exp(i);
end
toc