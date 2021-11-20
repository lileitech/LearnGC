clc;clear all;
foldtrain='C:\LiLei\2018MedAI_ScarSeg\Data\repeat_exp\TestData_exp3\';
a=dir([foldtrain '*_*']); 
cd (foldtrain); 

P_xy=6;
P_z=8;
rada=0.4;
xy=num2str(2*P_xy+1);
z=num2str(2*P_z+1);
r=num2str(rada);


for i=1:numel(a)
	CaseName=a(i).name; 
	PatchFold=[foldtrain  CaseName '\Patch_' xy '_' xy '_' z];%śÔÓŚ¸ü¸Ä

	%if exist(PatchFold,'dir')

	%------------------------Achieve the scar seg using Graph cut-------------------------
	command=['GraphCut_auto ' foldtrain ' ' CaseName ' ' num2str(P_z)  ' '  num2str(P_xy) ' ' r ];system(command);  

	%load Ground Truth and calculate Dice
    Groundlabel = strcat(foldtrain,CaseName,'\LA_MeshWallLabel_fixed.nii.gz');  %Ground truth
	nii = load_untouch_nii(Groundlabel);
	GLabel = uint16(nii.img);


 	oldSegScarlabel = strcat(PatchFold,'\ScarSegGraphCut.nii.gz'); 
	%SegScarlabel=strrep(oldSegScarlabel, 'GraphCut', 'GraphCut_shift6_reda04');
	%movefile(oldSegScarlabel,SegScarlabel); 
    %SegScarlabel=strcat(foldtrain,CaseName,'\en_seg_msp_la.nii.gz');
	nii = load_untouch_nii(oldSegScarlabel); 
	Seg_Scar_GraphCut = uint16(nii.img);
	Dice_GraphCut(i) = CalculateDice(Seg_Scar_GraphCut,GLabel);
	%end       
end  

xlswrite('C:\LiLei\2018MedAI_ScarSeg\Data\GDice_exp3',Dice_GraphCut(:));






