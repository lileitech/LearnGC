function results=SVGMMCardSegEM2D_new_lj(x,priorp_cla,nc,stopexp,iter_max)

nclass=length(nc);
ncomp=sum(nc);
xsize=size(x);

%=========initialization==============
classflag=zeros(sum(nc),1);%±êÖ¾Ã¿¸öcomponentËùÊôµÄÀà
%--devide the probability of each class into corresponding conponents
counter=1;
for i=1:nclass
    for j=1:nc(i)
         classflag(counter)=i;
         counter=counter+1;
    end
end
 %--Initializition on sigma and mu of each class
muclass=sum(sum(sum(bsxfun(@times,x,priorp_cla),1),2),3)./sum(sum(sum(priorp_cla,1),2),3);
sigmaclass= sqrt(sum(sum(sum((bsxfun(@minus,x,muclass).^2.*priorp_cla),1),2),3)./sum(sum(sum(priorp_cla,1),2),3));
muclass=reshape(muclass,nclass,1);
sigmaclass=reshape(sigmaclass,nclass,1);
%--End initializition on sigma and mu of each class

%--Initializition on sigma and mu of each component
mu_init=zeros(ncomp,1);
sigma_init=zeros(ncomp,1); 
counter=1;
for i=1:nclass
    if nc(i)==1
        mu_init(counter)=muclass(i);
    else
        mu_init(counter:counter+nc(i)-1)=-sigmaclass(i):(sigmaclass(i)*2)/(nc(i)-1):sigmaclass(i);
        mu_init(counter:counter+nc(i)-1)=mu_init(counter:counter+nc(i)-1)+muclass(i);
    end
    sigma_init(counter:counter+nc(i)-1)=sigmaclass(i)/sqrt(nc(i));
    counter=counter+nc(i);
end

%--Prapare parameters in EM iterations
% priorp_cla=prob;
postp_cla=zeros([xsize,nclass]);
delta=zeros(ncomp,1);
for i=1:ncomp
    delta(i)=1/nc(classflag(i));
end
delta=reshape(delta,1,1,1,[]);
mu1=reshape(mu_init,1,1,1,[]);
sigma1=reshape(sigma_init,1,1,1,[]);
% probatlas parameters - modify when needed

priorp_com=zeros([xsize,ncomp]);
Lx1=0;
dLx=10;

%========EM iterations================

iter_num=1;
while((dLx>=10^(-stopexp))&&(iter_num<=iter_max))    
    %%%%%%%%%%%%%%%%%%%%%%  E-Step  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     phi=gaussian1(x,mu1,sigma1);    
    for i=1:ncomp
        priorp_com(:,:,:,i)=priorp_cla(:,:,:,classflag(i));
    end
    priorp_com=bsxfun(@times,priorp_com,delta);
    postp_com=phi.*priorp_com;
    postp_com=bsxfun(@rdivide,postp_com,sum(postp_com,4));
    for i=1:nclass
        postp_cla(:,:,:,i)=sum(postp_com(:,:,:,classflag==i),4);
    end
    [~,I]=max(postp_cla,[],4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  M-step  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu2=sum(sum(sum(bsxfun(@times,postp_com,x),1),2),3)...
        ./sum(sum(sum(postp_com(:,:,:,1:end),1),2),3);
    sigma2=sqrt(sum(sum(sum(bsxfun(@times,postp_com,bsxfun(@minus,x,mu2).^2),1),2),3)...
        ./sum(sum(sum(postp_com(:,:,:,1:end),1),2),3));
    
    counter=1;
    for i=1:nclass
        temp=sum(sum(sum(postp_com(:,:,:,classflag==i),1),2),3);
        temp=temp/sum(temp);
        delta(counter:counter+nc(i)-1)=temp;
        counter=counter++nc(i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  log-likehood  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi=bsxfun(@rdivide,...
        exp(-0.5*bsxfun(@rdivide,bsxfun(@minus,x,mu2),sigma2).^2),...
        sqrt(2*pi)*sigma2)+eps;    
    
    Lx2=log(sum(postp_com.*phi,4));
    Lx2=sum(sum(sum(Lx2,1),2),3);
    dLx=abs((Lx2-Lx1)/Lx2);
    disp(dLx);

    Lx1=Lx2;
     mu1=mu2;
     sigma1=sigma2;

    iter_num=iter_num+1;
end
results=I;


function phi=gaussian1(x,mu,sigma)
    phi=bsxfun(@rdivide,...
        exp(-0.5*bsxfun(@rdivide,bsxfun(@minus,x,mu),sigma).^2),...
        sqrt(2*pi)*sigma)+eps;
return
