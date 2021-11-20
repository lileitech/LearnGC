
function a=Ostu(a)
   [m,n,s]=size(a); 
   N=m*n*s; 
   L=256; 
   
for i=1:L
    count(i)=length(find(a==(i-1)));
    f(i)=count(i)/(N);
end
 
for i=1:L 
    if count(i)~=0 
        st=i-1; 
        break; 
    end 
end 
for i=L:-1:1 
    if count(i)~=0 
        nd=i-1; 
        break; 
    end 
end 

p=st;   q=nd-st; 
u=0; 
for i=1:q 
    u=u+f(i)*(p+i-1);  
    ua(i)=u;           
end
 
for i=1:q 
    w(i)=sum(f(1:i));  
end
 
w=w+0.0001;
 
d=(u*w-ua).^2./(w.*(1-w)); 
[y,tp]=max(d); 
th=tp+p; 
 
for i=1:m 
    for j=1:n 
        for k=1:s
            if a(i,j,k)>th 
                a(i,j,k)=1; 
            else 
                a(i,j,k)=0; 
            end 
        end
    end 
end  
