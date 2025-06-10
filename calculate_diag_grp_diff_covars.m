function calculate_diag_grp_diff_covars
load CompleteData.mat
covars={'dprime_new','NPT_global_score','age','sex','edyears','EstimatedTotalIntraCranialVol','WM','LAN','MEM','EXEC','VIS','totaltau','ratio_Abeta42_40_prec','phosphotau181'};
statstbl=cell2table(cell(length(covars),1),"RowNames",covars,'VariableNames',{'result'});
for a=1:length(covars)
    [~,b,c]=anova1(CompleteData.(covars{a}),CompleteData.diag,'off');
F=round(b{contains(b(:,1),'Groups'),strcmp(b(1,:),'F')},2);
p=round(b{contains(b(:,1),'Groups'),strcmp(b(1,:),'Prob>F')},2);
diff=round(c.means(4)-c.means(1),2);
    statstbl{covars{a},1}={['F(3,156)=',num2str(F),', p=',num2str(p)]};
end
writetable(statstbl,'diagres_copy_paste.xlsx','WriteRowNames',true);

for a=1:length(covars)
    [~,b,c]=anova1(CompleteData.(covars{a}),CompleteData.diag,'off');
F=round(b{contains(b(:,1),'Groups'),strcmp(b(1,:),'F')},2);
p=round(b{contains(b(:,1),'Groups'),strcmp(b(1,:),'Prob>F')},2);
diff=round(c.means(4)-c.means(1),2);
    statstbl{covars{a},1}={['F(3,156)=',num2str(F),', p=',num2str(p),', AD-HC=',num2str(diff)]};
end

writetable(statstbl,'diagres.xlsx','WriteRowNames',true);

end
