function matrix = genmatrix(fid)

%fid='rb_bravo.txt';
data=dlmread(fid);
%zerolocations=find(data==0);
%data(zerolocations)=[];

matrix=zeros(length(data),max(data(:)));
for n=1:length(data)
    for m=1:max(data(:))
        if any(data(n,:)== m) == 1
            matrix(n,m)=1;
        else 
            matrix(n,m)=0;
        end
    end
end

