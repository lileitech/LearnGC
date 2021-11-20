
%   实验步骤
% 1.	Smooth LA
% 2.	生成LA_Mesh_M.stl
% 3.	生成LA_label_GauiisanBlur_M.nii.gz
% 4.	生成en_seg_msp_M.nii.gz
% 5.	生成LA_MeshWall_M.nii.gz
% 6.	生成LA_MeshWallLabel_M.nii.gz
% 7.	生成LA_PatchLabel_M.nii.gz
% 8.	生成dist image,从而生成LA_MeshWall_422_prob_M.nii.gz

clc;clear all;
verbose=' -v 0 ' ;
fstdForGenerateProb=2;

foldtrain='C:\LeiLi\2020_MICCAI\LA2020\Data\train_data\';
a=dir([foldtrain '*']); 
cd (foldtrain); 

for i=3:numel(a)
    CaseName=a(i).name; 
    AutoWHS=[foldtrain  CaseName '\en_seg_msp.nii.gz'];

    ManualLASeg=[foldtrain  CaseName '\atriumSegImgMO.nii.gz'];
    SmoothManualLASeg=[foldtrain  CaseName '\atriumSegImgMO_smoooth.nii.gz'];
    
    %------------------smooth the manual LASeg ------------------  
    command=['zxhvolumelabelop  ' ManualLASeg ' ' SmoothManualLASeg '   -sc -sg 0 0 3 '];  system(command); 
    command=['zxhimageop -int ' SmoothManualLASeg  ' zxhimageop -vr 1 100 -vs 420 ' verbose]; system(command);
    
    %------------------Generate the LA_Mesh_M.stl------------------
    command=['GenerateModelsFromLabels ' foldtrain ' ' CaseName]; system(command); 
    
    %------------------Generate LA_label_GauiisanBlur_M.nii.gz------------------
	LABlur_M=[foldtrain  CaseName  '\LA_label_GauiisanBlur_M.nii.gz'];  
    system(['zxhimageop -int ' SmoothManualLASeg ' -o ' LABlur_M ' -gau 4 -v 0 ']);
    
    % ------------------Generate manual WHS: en_seg_msp_M.nii.gz-----------------
    ManualLAAndAutoWHS=[foldtrain  CaseName '\en_seg_msp_M.nii.gz'];      
	command=['zxhimageop -int ' AutoWHS ' -o ' ManualLAAndAutoWHS ' zxhimageop -vr 420 420 -vs 0 ' verbose]; system(command); %remove from WHS?auto?
    command=['zxhimageop -int ' ManualLAAndAutoWHS ' -int ' SmoothManualLASeg ' -o ' ManualLAAndAutoWHS ' -sum ' verbose]; system(command);  %add Manual LA(smooth)   
    % 在mitral valve或其他位置上，如果出现空洞，用closing来填充上空洞的部分，注意选择一个合适的参数，不能让这个closing的范围“没有必要的”过大;  
	command=['zxhimageop -int ' ManualLAAndAutoWHS ' -CL 1.5 ' verbose]; system(command);   
    % 把叠加上Manual LA(smooth后的)后的，叠加label都归为420 （LA收纳后宫啦哈哈）
    command=['zxhimageop -int ' ManualLAAndAutoWHS ' -o ' ManualLAAndAutoWHS ' -vr 625 625 -vr 920 920 -vr 970 970 -vr 1240 1240 -vr 1270 1270 -vs 420 ' verbose]; system(command);
    
    %------------------生成LA_MeshWall_M.nii.gz-----------------
    command=['GenMeshWall ' foldtrain ' ' CaseName];   system(command); 
    
    %------------------生成LA_MeshWallLabel_M.nii.gz-----------------
    command=['GenMeshWallLabel ' foldtrain ' ' CaseName];system(command);  
    
%      %------------------生成LA_PatchLabel_M.nii.gz(找出交叉区域)-----------------
%     command=['GenPatchLabel ' foldtrain ' ' CaseName];system(command);  
  
     %------------------Generate dist file------------------
    enseg=[foldtrain a(i).name '\LA_MeshWallLabel_M.nii.gz'];%changed into fixed;
    verbose=' -v 0 ' ;
	img1=strrep(enseg, 'LA_MeshWallLabel_M', 'LA_MeshWall_test');
    WallImag=strrep(enseg, 'LA_MeshWallLabel_M', 'LA_MeshWall_421_M');
    ScarImag=strrep(enseg, 'LA_MeshWallLabel_M', 'LA_MeshWall_422_M');
    WallImagDismap=strrep(enseg, 'LA_MeshWallLabel_M', 'LA_MeshWall_421_dismap_M');
    ScarImagDismap=strrep(enseg, 'LA_MeshWallLabel_M', 'LA_MeshWall_422_dismap_M');
    
	command=['zxhimageop -int ' enseg ' -o ' img1 ' -vr 1 1000 -VS 1 ' verbose]; system(command);
	command=['zxhimageop -int ' enseg  ' -o ' WallImag ' -vr 421 421 -VS 421 ' verbose]; system(command); 
    command=['zxhimageop -int ' enseg  ' -o ' ScarImag ' -vr 422 422 -VS 422 ' verbose]; system(command);    
    command=['zxhvolumelabelop ' WallImag ' ' WallImagDismap ' -genmap 0 0 ' verbose]; system(command);      
    command=['zxhimageop -float ' WallImagDismap  ' -float ' img1 ' -o ' WallImagDismap  ' -mul' verbose]; system(command);    
    command=['zxhvolumelabelop ' ScarImag  '  '  ScarImagDismap ' -genmap 0 0 ' verbose]; system(command); 
    command=['zxhimageop -float ' ScarImagDismap  ' -float ' img1 ' -o ' ScarImagDismap  ' -mul' verbose]; system(command); 
    
    delete(img1,WallImag,ScarImag);
   
    %------------------------Generate MeshProb-------------------------
    command=['GenMeshProb ' foldtrain ' ' CaseName];system(command); 
    
%     % ------------------------generate prob images------------------------
%     imgLab = strcat(foldtrain,CaseName,'\en_seg_msp_M_test.nii.gz');
%     Prob = strcat(foldtrain,CaseName,'\Prob');
% 	command=['zxhimageop -int ' ManualLAAndAutoWHS ' -o ' imgLab  ' -vr 550 550 -vr 820 820 -vs 421 ' verbose]; system(command);%LA and Aorta are set 421
% 	command=['zxhimageop -int ' imgLab ' -vr 0 400 -vr 430 1000 -vs 0 ' verbose]; system(command); %others are set 0
% 	command=['zxhvolumelabelop ' imgLab ' ' Prob ' -genprobf ' num2str(fstdForGenerateProb) ' ' imgLab verbose]; system(command) ;
    
%     %----------------试图使用closing的办法，检测交界区域----------------
%     ClosingWHS=[foldtrain  CaseName '\en_seg_msp_Closing.nii.gz'];
%     JointRegionLA=[foldtrain  CaseName '\JointRegionLA.nii.gz'];   
%     command=['zxhimageop -int ' ManualLAAndAutoWHS ' -o ' ClosingWHS ' -CL 2 ' verbose]; system(command);
%     command=['zxhimageop -int ' ClosingWHS ' -int ' ManualLAAndAutoWHS ' -o ' JointRegionLA ' -sub ' verbose]; system(command);
%     command=['zxhimageop -int ' JointRegionLA ' -vr 0 400 -vr 551 800 -vr 821 1000 -vs 0 ' verbose]; system(command);

end       
     








