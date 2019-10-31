function meanB=bootstrap(data)
[NS t]=size(data);
meandata=mean(data);
%bootstrap context


B=10000; %number of iterations
p=repmat(1:NS,B,1);
p=p(reshape(randperm(B*NS),B,NS));


meandataB=zeros(length(data),B);

for(b=1:B)
    meandata_temp=  mean(data(p(b,:),:));
    meandataB(:,b)=meandata_temp;
end


meandataBmedian=mean(meandataB,2);
meandataBstd=std(meandataB,1,2);


meanB=meandataB;
end
