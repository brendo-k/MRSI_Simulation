function delete_temp
    proj = currentProject;
    root = proj.RootFolder;
    delete(append(root, '/temp/*'));
end