function sincmeup(fname, stepsize)

Y=spm_vol(fname);
ima=spm_read_vols(Y);
dim=size(ima);
[x, y, z]=ndgrid(1:1:dim(1), 1:1:dim(2), 1:1:dim(3));
[x2, y2, z2]=ndgrid(1:stepsize:dim(1), 1:stepsize:dim(2), 1:stepsize:dim(3));
new_tt=interpn(x, y, z, ima, x2, y2, z2,'sinc',0);
new_dim=size(new_tt);
newY=Y;
newY.dim=[new_dim(1) new_dim(2) new_dim(3)];
newY.fname=([Y.fname(1:end-4) '_interp.nii']);
spm_write_vol(newY, new_tt)