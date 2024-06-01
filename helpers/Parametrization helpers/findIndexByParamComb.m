function [boolindices] = findIndexByParamComb(allVarSelCell, UniqValsCell, keyvalsCell)
    boolindices = true; 
    for i=1:numel(allVarSelCell)
        boolindices = boolindices & findIndexByParamVal(allVarSelCell{i}, UniqValsCell{i}, keyvalsCell{i});
    end
end
