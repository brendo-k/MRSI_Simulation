% takes in FID-A metabolite structure name and extracts the metabolite name
function metaboliteName = extractMetaboliteName(metaboliteLabel)
    underscoreIndex = regexp(metaboliteLabel, '_');
    if(isempty(underscoreIndex))
        underscoreIndex = length(metaboliteLabel) + 1;
    end
    % get the metabolite name
    metaboliteName = metaboliteLabel(1:underscoreIndex - 1);
end