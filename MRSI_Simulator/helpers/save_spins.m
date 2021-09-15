function file_name = save_spins(spins, name)
tic
    file_name = ['temp/' char(datetime) '-' name];
    fileID = fopen(file_name,'w');
    fwrite(fileID, spins, 'single');
    fclose(fileID);
    toc
end