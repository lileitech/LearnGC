function [Accracy,Sensitivity,Selectivity,Dice_scar,GDice] = CalculateROC( QuerySegResult, QuerySegTrue)

T=421;

number_True_scar = sum(QuerySegTrue(:)>T);
number_Result_scar = sum(QuerySegResult(:)>T);
TPnumber=length(find((QuerySegTrue(:)+QuerySegResult(:))==2*(T+1)));
FNnumber=number_True_scar-TPnumber;
%FNnumber=length(find((QuerySegTrue(:)-QuerySegResult(:))==1)) + length(find((QuerySegTrue(:)-QuerySegResult(:))==T+1));
% FPnumber=length(find((QuerySegResult(:)-QuerySegTrue(:))==1)) + length(find((QuerySegResult(:)-QuerySegTrue(:))==T+1));
FPnumber=number_Result_scar-TPnumber;
TNnumber=length(find((QuerySegResult(:)+QuerySegTrue(:))==2*T)) + length(find(abs(QuerySegResult(:)-QuerySegTrue(:))==T));

N=TPnumber+FNnumber+FPnumber+TNnumber;
Accracy = (TPnumber+TNnumber)/N;
Sensitivity=TPnumber/(TPnumber+FNnumber);
Selectivity=TNnumber/(TNnumber+FPnumber);
Dice_scar=2*TPnumber/(TPnumber*2+FNnumber+FPnumber);

number_True_myo = sum(QuerySegTrue(:)==T);
number_Result_myo = sum(QuerySegResult(:)==T);
TPnumber_myo=length(find((QuerySegTrue(:)+QuerySegResult(:))==2*T));
Dice_myo=2*TPnumber_myo/(number_Result_myo+number_True_myo);

GDice=Dice_scar*number_True_scar/N+Dice_myo*number_True_myo/N;

% TPnumber=0; %+ + 
% FNnumber=0; %+ -
% FPnumber=0; %- +
% TNnumber=0; %+ -
%QueryImageSize = size(QuerySegResult);
% for k=1:QueryImageSize(3)
%     for j=1:QueryImageSize(2)
%         for i=1:QueryImageSize(1)
%             
%                if(QuerySegTrue(i,j,k)>T)&&(QuerySegResult(i,j,k)>T)
%                    TPnumber = TPnumber+1;  
%                end
% 
%                if(QuerySegTrue(i,j,k)>T)&&(QuerySegResult(i,j,k)<=T)
%                    FNnumber = FNnumber+1;
%                end
% 
%                if(QuerySegTrue(i,j,k)<=T)&&(QuerySegResult(i,j,k)>T)
%                    FPnumber = FPnumber+1;
%                end
% 
%                if(QuerySegTrue(i,j,k)<=T)&&(QuerySegResult(i,j,k)<=T)&&((QuerySegTrue(i,j,k)+QuerySegResult(i,j,k))>0)
%                    TNnumber = TNnumber+1;
%                end
%        
%         end
%     end
% end

% N=TPnumber+FNnumber+FPnumber+TNnumber;
% Accracy = (TPnumber+TNnumber)/N;
% Sensitivity=TPnumber/(TPnumber+FNnumber);
% Selectivity=TNnumber/(TNnumber+FPnumber);
% Dice=2*TPnumber/(TPnumber*2+FNnumber+FPnumber);
% 
% Dice_scar=2*TPnumber/(TPnumber*2+FNnumber+FPnumber);
% Dice_myo=2*TNnumber/(TNnumber*2+FNnumber+FPnumber);
% GDice=Dice_scar*(TPnumber+FNnumber)/N+Dice_myo*(FPnumber+TNnumber)/N;

end

