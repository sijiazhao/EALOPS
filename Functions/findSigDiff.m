function sigSignal=findSigDiff(data_B, pval)
% a function that looks at the bootstrapped data and for each point in time
% compuptes whether a certain percentage (X) of signals is above or below
% 0. the out put is a time series where all the significant points are set
% to 1 and the non significant points are set to 0.

[t NS]=size(data_B);
% t is length of each signal, NS is the number of bootstrap iterations.

sigSignal=zeros(1,t);
for (i=1:t)
    AllIt=data_B(i,:);
    numPos=length(find(AllIt>0));
    numNeg=length(find(AllIt<0));
    
    percntPos=numPos/NS;
    percntNeg=numNeg/NS;
    
    if((percntPos<pval)&& (abs(mean(AllIt))>0))
          sigSignal(i)=-1;
    end
    if ((percntNeg<pval)&& (abs(mean(AllIt))>0))
        sigSignal(i)=1;
    end
    
end

sigSignal(find(sigSignal==0))=NaN;
end