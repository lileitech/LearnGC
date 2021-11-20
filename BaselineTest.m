clc;clear all;
foldtrain='C:\LiLei\2020_MICCAI\LA2020\Data\fold_result\';
a=dir([foldtrain 'p*']);  
verbose='-v 0';

for i=1:numel(a)
    CaseName=a(i).name; 
    TargetImage = [foldtrain  CaseName '\enhanced.nii.gz'];
    %TargetLabel = [foldtrain  CaseName '\LAwall_gd.nii.gz'];
    TargetLabel = [foldtrain  CaseName '\LAwall_Unet.nii.gz'];
    
%     LAU_Net = [foldtrain  CaseName '\LA_predict_Unet_BCE.nii.gz'];
%     img2 = [foldtrain  CaseName '\img2.nii.gz'];
%     command=['zxhimageop -int ' LAU_Net ' -o ' img2 ' -DIs 2.5 0 1 ' verbose]; system(command); 
% 	command=['zxhimageop -int ' LAU_Net ' -o ' TargetLabel ' -ERs 1.5 0 1 ' verbose]; system(command);   
% 	command=['zxhimageop -int ' img2 ' -int ' TargetLabel ' -o ' TargetLabel ' -sub ' verbose]; system(command);
%     delete(img2)
    
    im_nii=load_untouch_nii(TargetImage);
    img=double(im_nii.img);   
    lab_nii=load_untouch_nii(TargetLabel);
    lab=double(lab_nii.img); 
    mseg = Ostu(img.*logical(lab));
    im_nii.img=int16(mseg);
    
%     foldtrain_new = strrep(foldtrain, 'fold_3', 'fold_result');
%     newdir = [foldtrain_new CaseName];
%     if ~exist(newdir,'dir')
%         b=['mkdir ' ' ' newdir];
%         system(b); 
%     end
    SaveName=[foldtrain CaseName '\OstuSeg_Unet.nii.gz'];
    
    save_untouch_nii(im_nii,SaveName);
    

%     [Percentage_Ostu(i),PercentageTrue(i)]=CalculatePercentage(Seg,GLabel);
%     [Accracy_GraphCut(i),Sensitivity_GraphCut(i),Selectivity_GraphCut(i),Dice_GraphCut(i)] = CalculateROC(Seg,GLabel);

    
end
 

% SaveResult=[Accracy_GraphCut(:),Sensitivity_GraphCut(:),Selectivity_GraphCut(:),Dice_GraphCut(:),PercentageTrue(:),Percentage_Ostu(:)];       
% xlswrite([foldtrain '\EvaluatedResult_GMM'],SaveResult);

 
