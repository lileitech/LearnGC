clc;clear all;

foldtrain='C:\LiLei\2020_MICCAI\LA2020\Data_patch\test_data\';
a=dir([foldtrain '*_*']); 

P_xy=6;
P_z=8;
xy=num2str(2*P_xy+1);
z=num2str(2*P_z+1);
weight = num2str(0);

% for i=1:numel(a)
%     CaseName=a(i).name;
%     %PatchFold=[foldtrain CaseName '\Patch_' xy '_' xy '_' z '_MS'];
%     PatchFold=[foldtrain CaseName '\Patch_' xy '_' xy '_' z];
%     %if ~exist(PatchFold,'dir')
%     command=['PatchExtractForTest ' foldtrain ' ' CaseName ' ' num2str(P_z)  ' '  num2str(P_xy)]; %_MS_M
%     system(command);%_MS_M        
%     %end
%     delete([PatchFold '\Patch_N_info.txt_____temp.txt'])
%     delete([PatchFold '\Patch_T_info.txt_____temp.txt'])
%     delete([PatchFold '\Patch_T_info_norm.txt_____temp.txt'])
% end

for i=1:numel(a)
    CaseName=a(i).name;  
    PatchFold=[foldtrain CaseName '\Patch_' xy '_' xy '_' z '_MS'];
    %PatchFold=[foldtrain CaseName '\Patch_' xy '_' xy '_' z];
    command=['GraphCut ' foldtrain ' ' CaseName ' ' num2str(P_z)  ' '  num2str(P_xy) ' ' weight]; 
    system(command);        
end


