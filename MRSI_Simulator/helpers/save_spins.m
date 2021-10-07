%Take spin matrix from MRSI object and save it to a binary file

function file_name = save_spins(spins, name, file_name)
proj = currentProject;
root = proj.RootFolder;
file_name = append(root, ['/temp/' char(datetime) '-' name '-' file_name]);
fileID = fopen(file_name,'w');
fwrite(fileID, spins, 'single');
fclose(fileID);
end