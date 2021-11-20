clc;clear all;
foldtrain='C:\LeiLi\2020_MICCAI\LA2020\Data_60\test_data_result\';
a=dir([foldtrain 'p*']); 
cd (foldtrain); 

P_xy=6;
P_z=8;
rada=0;
xy=num2str(2*P_xy+1);
z=num2str(2*P_z+1);
r=num2str(rada);

for i=1:numel(a)
    CaseName=a(i).name; 
    PatchFold=[foldtrain  CaseName '\Patch_' xy '_' xy '_' z '_MS'];

%     if exist(PatchFold,'dir')
    %load Ground Truth and calculate Dice
    Groundlabel = strcat(foldtrain,CaseName,'\scarSegImgM_surface.nii.gz');  %Ground truth
    nii = load_untouch_nii(Groundlabel);
    GLabel = int16(nii.img);
    
    SegScarlabel=strcat(foldtrain,CaseName,'\scar_predict_dice_base_surface.nii.gz');
%     SegScarlabel_new = strcat(foldtrain,CaseName,'\test.nii.gz');
%     command=['zxhimageop -int ' SegScarlabel ' -o ' SegScarlabel_new ' -vr 1 1 -vs 420 ']; system(command);
    nii = load_untouch_nii(SegScarlabel); 
    Seg_Scar_GraphCut = nii.img;
    %Dice_GraphCut(i) = CalculateDice(Seg_Scar_GraphCut,GLabel);
    %delete(SegScarlabel_new)
    [Accracy_GraphCut(i),Sensitivity_GraphCut(i),Selectivity_GraphCut(i),Dice_GraphCut(i),GDice_GraphCut(i)] = CalculateROC(Seg_Scar_GraphCut,GLabel);
    
%     end       
end     

%xlswrite('C:\LeiLi\2020_MICCAI\LA2020\Dice_2d',Dice_GraphCut(:));

SaveResult=[Accracy_GraphCut(:),Sensitivity_GraphCut(:),Selectivity_GraphCut(:),Dice_GraphCut(:),GDice_GraphCut(:)];       
xlswrite('C:\LeiLi\2020_MICCAI\LA2020\GMM2',SaveResult);








