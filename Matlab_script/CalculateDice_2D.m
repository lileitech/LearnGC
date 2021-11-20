function Dice = CalculateDice_2D( QuerySegResult, QuerySegTrue)

TPnumber=0; %+ + 
FNnumber=0; %+ -
FPnumber=0; %- +
TNnumber=0; %+ -

T=1;
QueryImageSize = size(QuerySegResult);
Dice_sum = 0;
for i=1:QueryImageSize(3)

    number_True_scar = sum(QuerySegTrue(:, :, i)==T);
    number_Result_scar = sum(QuerySegResult(:, :, i)==T);
    TPnumber=length(find((QuerySegTrue(:, :, i)+QuerySegResult(:, :, i))==2*T));
    FNnumber=number_True_scar-TPnumber;
    FPnumber=number_Result_scar-TPnumber; 
    
    Dice_slice=2*TPnumber/(TPnumber*2+FNnumber+FPnumber);
    Dice_sum = Dice_sum + Dice_slice;   
end

Dice = Dice_sum/QueryImageSize(3);

end
