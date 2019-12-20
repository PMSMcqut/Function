clear;clc;
openfemm;
opendocument('SFA_LIM_prototype.fem');
mi_saveas('temp.fem');
for k=1:1:20
    freq=k/4;
    mi_probdef(freq,'millimeters','planar',1e-8,148,20);
    mi_analyze(1);
    mi_loadsolution;
    ll=mo_getcircuitproperties('A');
    flux=ll(3);
    mo_groupselectblock(4);
    Trustforce=mo_blockintegral(18);
    fid = fopen(['fluxlinkage','.txt'],'a');
    fprintf(fid, '%4.8f\t',freq,flux,Trustforce);
    fprintf(fid, '\n');
    fclose(fid);
end
closefemm;
delete('temp.fem');
delete('temp.ans');

