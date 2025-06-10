function calculate_manual_auto_ICC
if strcmp(detectOS,'nix')
    matlabdrive='/storage/DataTempVolatile/Friedrich/MATLAB-Drive';
end
relvals={'vol_man','vol_auto','contrast_man','contrast_auto','dice','specificity','sensitivity'};
basedir='/DataTempVolatile/Friedrich/final/Betts82_0.75/';
SN_masks=mydir(fullfile(basedir,'raw'),'_man.nii',2);
img=cellfun(@(x) x(1:end-11),SN_masks,'UniformOutput',false);
stats=array2table(NaN(length(img),length(relvals)),'RowNames',img,VariableNames=relvals);
load CompleteData.mat

for a=1:length(img)
    disp(['processing ',img{a}])
    if exist(fullfile(basedir,[img{a},'_SN_man.nii']),'file')==0
        continue
    end
    
    main=mynifti(fullfile(basedir,[img{a},'.nii']));
    pons=mynifti(fullfile(basedir,[img{a},'_Pons.nii']));
    SN_man=mynifti(fullfile(basedir,[img{a},'_SN_man.nii']));
    SN_man.Values(SN_man.Values>0.2)=1;
    SN=mynifti(fullfile(basedir,[img{a},'_SN.nii']));
    SN.Values(SN.Values>0.2)=1;
    SN_ROI=mynifti(SN_mask{contains(SN_mask,img{a})&contains(SN_mask,'SN_ROI')});
    tp=sum(sum(sum(SN_man.Values>0.2&SN.Values>0.2)));
    tn=sum(sum(sum(SN_man.Values==0&SN.Values==0&SN_ROI.Values>0.2)));
    fp=sum(sum(sum(SN_man.Values==0&SN.Values>0.2)));
    fn=sum(sum(sum(SN_man.Values>0.2&SN.Values==0)));
    stats{img{a},'sensitivity'}=tp./(tp+fn);
    stats{img{a},'specificity'}=tn./(tn+fp);
    stats{img{a},'dice'}=max(dice(SN_man.Values,SN.Values));
    stats{img{a},'vol_auto'}=sum(sum(sum(SN.Values>0.2)))*0.75*0.75*0.75;
    stats{img{a},'vol_man'}=sum(sum(sum(SN_man.Values>0.2)))*0.75*0.75*0.75;
    ponsauto=median(main.Values(pons.Values>0.2));
    medauto=median(main.Values(SN.Values>0.2));
    stats{img{a},'contrast_auto'}=(medauto-ponsauto)./ponsauto;
    medman=median(main.Values(SN_man.Values>0.2));
     stats{img{a},'contrast_man'}=(medman-ponsauto)./ponsauto;
end
stats(isnan(stats.vol_man),:)=[];
volstats = ICC(stats{:,{'vol_auto','vol_man'}},'A-1',0.05);
constats=ICC(stats{:,{'contrast_auto','contrast_man'}},'A-1',0.05);
save(fullfile(fileparts(which('calculate_manual_auto_ICC.m')),'stats_betts_manual_automatic.mat'),'stats','volstats','constats');
end