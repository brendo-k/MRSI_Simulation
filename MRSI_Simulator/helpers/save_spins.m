%Take spin matrix from MRSI object and save it to a binary file

function file_name = save_spins(spins, name, file_name)
filePath = mfilename('fullpath');
[pathRoot, ~, ~] = fileparts(filePath);
tempFileName = ['/temp/' char(datetime) '-' name '-' file_name];
file_name = fullfile(pathRoot, '../../', tempFileName);
fileID = fopen(file_name,'w');
fwrite(fileID, spins, 'single');
fclose(fileID);
end