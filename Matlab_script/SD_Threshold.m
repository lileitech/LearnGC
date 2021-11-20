function a=SD_Threshold(a, sd_value)   

    Mean_T=sum(a(:))/nnz(a);

    [m,n,s]=size(a); 
    N=m*n*s; 
    sum_P=0;
    for i=1:m 
        for j=1:n 
            for k=1:s
                if a(i,j,k)>0 
                    sum_P=sum_P+double((a(i,j,k)-Mean_T)^2);
                end 
            end
         end 
    end
   
   std_T=sqrt(sum_P/nnz(a));
   
   th=Mean_T+sd_value*std_T;
   
   
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
