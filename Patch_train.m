clc;clear all;

foldtrain='C:\LiLei\2020_MICCAI\LA2020\Data_patch\train_data\';
a=dir([foldtrain 'p*']); 

P_xy=6;
P_z=8;
RandShift=P_z;
xy=num2str(2*P_xy+1);
z=num2str(2*P_z+1);
h=num2str(RandShift);

for i=1:numel(a)
    CaseName=a(i).name; 
    
%     b=dir([foldtrain CaseName '\*predict*']); 
%     for j=1:numel(b)
%         delfile = [foldtrain CaseName '\' b(j).name]; 
%         delete(delfile)
%     end
  
    %PatchFold=[foldtrain CaseName '\Patch_' xy '_' xy '_' z '_MS'];
    PatchFold=[foldtrain CaseName '\Patch_' xy '_' xy '_' z];
    if ~exist(PatchFold,'dir')
        command=['PatchExtractForTraining ' foldtrain ' ' CaseName ' ' num2str(P_z)  ' '  num2str(P_xy) ' ' h];system(command); %_MS_M
    end %_MS_M
end
