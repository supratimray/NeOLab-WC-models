function [boolindices] = findIndexByParamVal(varSel, OrderedUniqVals, keyval)
    if isempty(keyval)
        boolindices = any(varSel,2);
    else
        boolindices = any(varSel(:,any(keyval(:)==OrderedUniqVals(:)',1))==1, 2);
        boolindices = boolindices(:)==1;
    end
end
