%load spins from binary file
% 
function spins = load_spins(phantom, m)

%open file 
file_id = fopen(phantom.file{m});
%read file as single data type and cast to single
spins = fread(file_id,inf, 'single=>single');
%close file and delete
%delete(file_name);
fclose(file_id);

%reshape spins back into proper dimensions
spins = reshape(spins, [size(phantom.d{m}, [1,2]), length(phantom.y), length(phantom.x)]);
end