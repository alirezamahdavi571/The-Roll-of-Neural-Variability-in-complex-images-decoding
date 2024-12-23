function [class , Mdl] = gen_fx_MC_LDA(sample, train, grp)

    classNo = double(unique(grp));
    cls_nu = length(classNo);
    r = nchoosek(1:cls_nu, 2);
    vote = [];

    % Remove zero variance features
    featureVariance = var(train);
    nonZeroVarIdx = featureVariance > 0;
    train = train(:, nonZeroVarIdx);
    sample = sample(:, nonZeroVarIdx);

    for i = 1:length(r(:, 1))
        ix = (grp == classNo(r(i, 1))) | (grp == classNo(r(i, 2)));
        % Use pseudoLinear to handle zero variance within class
        svmStruct = fitcdiscr(train(ix, :), grp(ix), 'DiscrimType', 'pseudoLinear');
        vote = [vote predict(svmStruct, sample)];
    end
    class = mode(vote, 2);

    % Train final model
    Mdl = fitcdiscr(train, grp, 'DiscrimType', 'pseudoLinear');
end
