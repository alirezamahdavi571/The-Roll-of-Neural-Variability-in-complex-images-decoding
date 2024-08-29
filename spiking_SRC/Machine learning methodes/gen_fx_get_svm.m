function out = gen_fx_get_svm(grp,I,rate,rep)
Pt = [];
Tu = [];
for cnt = 1:rep
    [test train] = gen_fx_get_equal_part(grp,rate);
%     if sum(test)>5 && sum(train) > 5
        cls = gen_fx_MC_SVM(I(test,:),I(train,:),grp(train));
        Pt = [Pt; sum(cls == grp(test))/sum(test)];
        Ct(:,:,cnt)= confusionmat(cls,grp(test));
        Tu = [Tu diag(Ct(:,:,cnt)) ./ (sum(Ct(:,:,cnt)))'];
%     else
%         Pt = [Pt; nan];
%         Ct(:,:,cnt)= nan;
%         Tu = [Tu nan];
%     end
end
out.C = Ct; out.pt = Pt; out.tu = Tu;
end