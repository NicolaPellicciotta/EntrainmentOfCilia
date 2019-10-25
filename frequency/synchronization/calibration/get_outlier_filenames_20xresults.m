function outlier_filenames = get_outlier_filenames_20xresults(C)

outlier_filenames = {}
outlier_counter = 1
for ii = 1:numel(C)
    if C{ii}.A1 == 0 && ~contains(C{ii}.filename, "P0")
        outlier_filenames{outlier_counter} = C{ii}.filename;
        outlier_counter = outlier_counter + 1;
    end
end

end