clc;clear all;

foldtrain='C:\LiLei\2018MedAI_ScarSeg\Data\repeat_exp\data\'; 
a=dir([foldtrain '*_*']); 

trainingdir=strrep(foldtrain, 'data', 'TrainingData_exp5');
testdir=strrep(foldtrain, 'data', 'TestData_exp5');

p = randperm(numel(a));
p1 = p(1:27);
p2 = p(27+1:end);

for i=1:27   
    CaseName=a(p1(i)).name; 
    oldname=[foldtrain CaseName];
    newname=[testdir CaseName];
  	copyfile(oldname,newname);  
end

for i=1:31   
    CaseName=a(p2(i)).name; 
    oldname=[foldtrain CaseName];
    newname=[trainingdir CaseName];
  	copyfile(oldname,newname);  
end



