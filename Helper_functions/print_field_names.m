function tbOut = print_field_names(tbIn)

    tbOut = table(tbIn.Properties.VariableNames', 'VariableNames', {'fieldNames'});
    disp(' ');
    disp(tbOut)

end