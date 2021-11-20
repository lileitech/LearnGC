function Dice = CalculateDice( QuerySegResult, QuerySegTrue)

T=421;

%TPnumber=length(find((QuerySegTrue(:)+QuerySegResult(:))>=2));
TPnumber=length(find((QuerySegTrue(:)+QuerySegResult(:))==2*T+2));
number_True = sum(QuerySegTrue(:)>T);
number_Result = sum(QuerySegResult(:)>T);

%Dice=2*TPnumber/(number_Result+number_True);

Dice = number_Result/sum(QuerySegResult(:)==T);
end

