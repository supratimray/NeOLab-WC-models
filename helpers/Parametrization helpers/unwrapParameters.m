function [varVals, varSels, nsims] = unwrapParameters(varuniq)
nparams=numel(varuniq);
vargrid = cell(1,nparams);
for i=1:nparams
    vargrid{i} = varuniq{i};
end

numvals = zeros(nparams,1);
varSels = cell(nparams,1);
for i=1:nparams
    if i>1
        for j=1:i-1
            [vargrid{j},Z] = meshgrid(vargrid{j}(:)', vargrid{i}(:)');
            vargrid{j} = vargrid{j}(:);
        end
        vargrid{i} = Z(:);
    end
    numvals(i) = numel(vargrid{i});
end
for i=1:nparams
    varSels{i} = vargrid{i}(:)==varuniq{i}(:)';
end

varVals = cell2mat(vargrid);
nsims = numel(varVals)/nparams;
end

