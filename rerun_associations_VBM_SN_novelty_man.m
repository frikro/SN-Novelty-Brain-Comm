function rerun_associations_VBM_SN_novelty_man
%% prep
clear matlabbatch
load('CompleteData')
switch detectOS
    case 'nix'
        fullscans=cellfun(@(x) ['/DataTempVolatile/Friedrich/DELCODE/novelty/raw_maps/',x,'.nii'],...
            CompleteData.SubjID,'UniformOutput',false);
end
covlistcomplete={'sex','site_id','EstimatedTotalIntraCranialVol','edyears','age'};
nosite={'sex','EstimatedTotalIntraCranialVol','edyears','age'};
%% SN median - abeta pos
diagidx=CompleteData.ratio_Abeta42_40_prec<0.08;
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_abeta_pos';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = Completep.SNmednorm;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,covlistcomplete);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end
%% SN median - CSF data
diagidx=~isnan(CompleteData.ratio_Abeta42_40_prec);
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_CSF';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end

if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = Completep.SNmednorm;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,covlistcomplete);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end


%% SN median - abeta pos - harmonized
diagidx=CompleteData.ratio_Abeta42_40_prec<0.08;
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_abeta_pos_harm';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c =zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,nosite);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end
%% SN median - CSF data adjusted
diagidx=~isnan(CompleteData.ratio_Abeta42_40_prec);
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_CSF_harm';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,nosite);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end



%% SN median - CSF data adjusted HC volume
diagidx=~isnan(CompleteData.ratio_Abeta42_40_prec);
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_CSF_harm_HC';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,[nosite,'HCnorm']);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end

%% SN median - CSF data adjusted HC volume
diagidx=~isnan(CompleteData.ratio_Abeta42_40_prec);
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_CSF_harm_diag';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,[nosite,'diag']);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end

%% SN median - CSF data adjusted grey matter
diagidx=~isnan(CompleteData.ratio_Abeta42_40_prec);
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_CSF_harm_GM_vol';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch=return_covs_lm(matlabbatch,Completep,[nosite,'TotalGrayVol']);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end

%% SN median - impaired adjusted
diagidx=CompleteData.diag>2;
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_impaired';
clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch=return_covs_lm(matlabbatch,Completep,nosite);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end


%% SN median - unimpaired adjusted
diagidx=CompleteData.diag<3;
spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_unimpaired';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans(diagidx);
Completep=CompleteData(diagidx,:);
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch=return_covs_lm(matlabbatch,Completep,nosite);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
end

%% SN median - complete

spmdir='/DataTempVolatile/Friedrich/DELCODE/VBM/SNmed_complete_harm';

clear matlabbatch

if exist(spmdir,'dir')==0
    mkdir(spmdir)
end
if isempty(mydir(spmdir))
scans=fullscans;
Completep=CompleteData;
matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = scans;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = zscore(Completep.SN_Med_ratio_bilat_adj);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'SN median';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch=return_covs_lm(matlabbatch,Completep,nosite);

matlabbatch=add_final_touch_lm(matlabbatch);
spm_jobman('run',matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);

end



end