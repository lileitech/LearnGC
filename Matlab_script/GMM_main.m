clc;clear all;
foldtrain='C:\LiLei\2020_MICCAI\LA2020\Data\fold_result\';
a=dir([foldtrain 'p*']); 
cd (foldtrain); 

nc =[1;2];
stopexp=2;
max_iter=100;
xmin=0;xmax=1000;
fstdForGenerateProb=1;
verbose=' -v 0 ' ;

%============================= read files =================================
for i=1:numel(a)
    CaseName=a(i).name; 
    TargetImage = strcat(foldtrain,CaseName,'\enhanced.nii.gz'); 
    im_nii=load_untouch_nii(TargetImage);% nii headers
    temp=double(im_nii.img); 
    tmax=max(temp(:));tmin=min(temp(:));
    temp=(temp-tmin)*(xmax-xmin)/(tmax-tmin)+xmin;
    originalsize=size(temp);
    x=temp;

    % generate prob images
    Prob = [foldtrain  CaseName '\Prob'];
    imgLab = [foldtrain,CaseName,'\OstuSeg_Unet.nii.gz']; 
    command=['zxhvolumelabelop ' imgLab ' ' Prob ' -genprobf ' num2str(fstdForGenerateProb) ' ' imgLab verbose]; 
    system(command) ;

    probatlasmin=0;
    probatlasmax=1 ;

    ProbImg=[foldtrain,CaseName,'\Prob_Label1.nii.gz']; 
    pa_nii=load_untouch_nii(ProbImg);% nii headers
    probatlas=double(pa_nii.img); 
    pmin=min(probatlas(:)); pmax=max(probatlas(:)) ;
    if pmax-pmin>0
        probatlas=(probatlas-pmin)*(probatlasmax-probatlasmin)/(pmax-pmin)+probatlasmin;
        prob=zeros([originalsize,2]);
        prob(:,:,:,1)=probatlas;
        prob(:,:,:,2)=1-probatlas(:,:,:,1);
    end


    mseg=SVGMMCardSegEM2D_new_lj(x,prob,nc,stopexp,max_iter);
    im_nii.img=int16(mseg);
    SaveName=strcat(foldtrain,CaseName,'\GMMSeg_Unet.nii.gz');
    save_untouch_nii(im_nii,SaveName); 
    command=['zxhimageop -int ' SaveName ' -vr 2 2 -vs 0 ' verbose]; system(command); 
    delete([foldtrain,CaseName,'\Prob_Label0.nii.gz'])
    delete([foldtrain,CaseName,'\Prob_Label1.nii.gz'])
end




 


