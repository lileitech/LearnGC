function Dice = CalculateDice_G( QuerySegResult, QuerySegTrue)

TPnumber=0; %+ + 
FNnumber=0; %+ -
FPnumber=0; %- +
TNnumber=0; %+ -

QueryImageSize = size(QuerySegResult);

for k=1:QueryImageSize(3)
    for j=1:QueryImageSize(2)
        for i=1:QueryImageSize(1)
            
		   if(QuerySegTrue(i,j,k)>421)&&(QuerySegResult(i,j,k)>421)
               TPnumber = TPnumber+1;  
           end
		   
		   if(QuerySegTrue(i,j,k)>421)&&(QuerySegResult(i,j,k)==421)
               FNnumber = FNnumber+1;
           end
		   
		   if(QuerySegTrue(i,j,k)==421)&&(QuerySegResult(i,j,k)>421)
               FPnumber = FPnumber+1;
           end
		   
           if(QuerySegTrue(i,j,k)==421)&&(QuerySegResult(i,j,k)==421)
               TNnumber = TNnumber+1;
           end
       
        end
    end
end

N=TPnumber+FNnumber+FPnumber+TNnumber;
Dice_scar=2*TPnumber/(TPnumber*2+FNnumber+FPnumber);
Dice_myo=2*TNnumber/(TNnumber*2+FNnumber+FPnumber);
Dice=Dice_scar*(TPnumber+FNnumber)/N+Dice_myo*(FPnumber+TNnumber)/N;
end

