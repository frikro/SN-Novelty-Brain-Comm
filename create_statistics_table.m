function create_statistics_table
%CREATE_STATISTICS_TABLE Summary of this function goes here
%   Detailed explanation goes here
cd(fileparts(which('create_statistics_table.m')))
load CompleteData
CompleteData.SN_Volume_man_ts_mask_adjusted_new=zscore(CompleteData.SN_Vol_bilat_adj,[],'omitnan');
CompleteData.SN_med_ratio_pons_raw_adjusted_new=zscore(CompleteData.SN_Med_ratio_bilat_adj,[],'omitnan');
covars={'age','sex','edyears','EstimatedTotalIntraCranialVol'};
for a=1:length(covars)
    if length(unique(CompleteData.(covars{a})))>5
        CompleteData.(covars{a})=zscore(CompleteData.(covars{a}),[],'omitnan');
    end
    if a==1
        covarsstr=[covars{a},'+'];
    else
        covarsstr=[covarsstr,covars{a},'+'];
    end
end
CompleteData.sex=categorical(CompleteData.sex);
tmpmdl=fitlm(CompleteData,['SN_Volume_man_ts_mask_adjusted_new~',covarsstr,'SN_med_ratio_pons_raw_adjusted_new']);
snvolres=['r^2 = ',num2str(round(tmpmdl.Rsquared.Adjusted,3)),' p=',num2str(tmpmdl.Coefficients{"SN_med_ratio_pons_raw_adjusted_new","pValue"}),' n= ',num2str(tmpmdl.NumObservations)];


cogfactors={'dprimenorm','NPT_global_score','VIS','MEM','WM','EXEC','LAN'};
adfactors={'ratio_Abeta42_40_prec','totaltau','phosphotau181'};
for a=1:length(adfactors)
    CompleteData.(adfactors{a})=zscore(CompleteData.(adfactors{a}),[],'omitnan');
end


for a=1:length(cogfactors)
    CompleteData.(cogfactors{a})=zscore(CompleteData.(cogfactors{a}),[],'omitnan');
end
CompleteData.TotalGrayVol=zscore(CompleteData.TotalGrayVol);
addcovars={'all','diag','HCvol','GM_vol','HC','SCD','MCI','ADD'};
snvars={'SN_Volume_man_ts_mask_adjusted_new','SN_med_ratio_pons_raw_adjusted_new'};
cogtests=cell(length(addcovars).*length(cogfactors),1);
idx=1;
for a=1:length(cogfactors)
    for b=1:length(addcovars)
        cogtests{idx}=[cogfactors{a},'_',addcovars{b}];
        idx=idx+1;
    end
end
assres=cell2table(cell(length(covars),length(snvars)),"RowNames",covars,VariableNames=snvars);
cogres=cell2table(cell(length(cogtests),length(snvars)),"RowNames",cogtests,VariableNames=snvars);
adred=cell2table(cell(length(adfactors),length(snvars)),"RowNames",adfactors,VariableNames=snvars);
format_mdl=@(x,y) {['R^2 = ',num2str(round(x.Rsquared.Adjusted,3)),', p=',num2str(round(x.Coefficients{y,"pValue"},3)),', n= ',num2str(x.NumObservations)]};

for a=1:length(snvars)


    for b=1:length(covars)
        idx=1;

        for c=1:length(covars)
            if strcmp(covars{b},covars{c})
                continue
            end
            if idx==1
                tmpcovarsstr=[covars{c},'+'];
            else
                tmpcovarsstr=[tmpcovarsstr,covars{c},'+'];
            end
            idx=idx+1;
        end
        if ~strcmp(covars{b},'sex')
            tmpmdl=fitlm(CompleteData,[covars{b},'~',tmpcovarsstr,snvars{a}]);
            assres{covars{b},snvars{a}}=format_mdl(tmpmdl,snvars{a});
        else
            [~,p,~,stats]=ttest2(CompleteData{CompleteData.sex=="0",snvars{a}},CompleteData{CompleteData.sex=="1",snvars{a}});
            assres{covars{b},snvars{a}}={['t(',num2str(stats.df),')=',num2str(round(stats.tstat,2)),', p=',num2str(round(p,3))]};
        end
    end


    for b=1:length(adfactors)
        tmpmdl=fitlm(CompleteData,[snvars{a},'~',covarsstr,adfactors{b}]);
        adred{adfactors{b},snvars{a}}=format_mdl(tmpmdl,adfactors{b});
    end


    for b=1:length(cogfactors)
        for c=1:length(addcovars)
            switch addcovars{c}
                case 'all'
                    tmpmdl=fitlm(CompleteData,[snvars{a},'~',covarsstr,cogfactors{b}]);
                case 'diag'
                    tmpmdl=fitlm(CompleteData,[snvars{a},'~',covarsstr,'diag+',cogfactors{b}]);
                case 'HCvol'
                    tmpmdl=fitlm(CompleteData,[snvars{a},'~',covarsstr,'HCnorm+',cogfactors{b}]);
                case 'GM_vol'
                    tmpmdl=fitlm(CompleteData,[snvars{a},'~',covarsstr,'TotalGrayVol+',cogfactors{b}]);
                case 'HC'
                    tmpmdl=fitlm(CompleteData(CompleteData.diag==1,:),[snvars{a},'~',covarsstr,cogfactors{b}]);
                case 'SCD'
                    tmpmdl=fitlm(CompleteData(CompleteData.diag==2,:),[snvars{a},'~',covarsstr,cogfactors{b}]);
                case 'MCI'
                    tmpmdl=fitlm(CompleteData(CompleteData.diag==3,:),[snvars{a},'~',covarsstr,cogfactors{b}]);
                case 'ADD'
                    tmpmdl=fitlm(CompleteData(CompleteData.diag==4,:),[snvars{a},'~',covarsstr,cogfactors{b}]);
            end
            cogres([cogfactors{b},'_',addcovars{c}],snvars{a})=format_mdl(tmpmdl,cogfactors{b});
        end
    end
end
allres=[assres;cogres;adred];
relvals=[covars,cellfun(@(x) [x,'_all'],cogfactors,'UniformOutput',false),adfactors]';

pvals=cellfun(@(x) sscanf(x(strfind(x,'p=')+2:strfind(x,'p=')+6),'%f'),table2cell(allres(relvals,:)),'UniformOutput',true);
p_val_corr=round([fdr_BH(pvals(:,1),0.05)',fdr_BH(pvals(:,2),0.05)'],3);

for a=1:length(snvars)

    for b=1:length(relvals)
        tmp=cell2mat(allres{relvals{b},snvars{a}});
        if contains(tmp,'n=')
                    start=tmp(1:strfind(tmp,', n='));
            tmpnew=[start,' q=',num2str(p_val_corr(b,a)),tmp(strfind(tmp,', n='):end)];
        else
            tmpnew=[tmp,', q= ',num2str(p_val_corr(b,a))];

        end
       allres{relvals{b},snvars{a}}={tmpnew};
    end
end
writetable(allres,'all_results.xlsx','WriteRowNames',true,'WriteVariableNames',true)

end

