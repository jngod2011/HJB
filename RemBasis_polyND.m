function Pf=RemBasis_polyND(Pf,N)


PP=[];
for i =1:1:size(Pf,1)
    p=sum(Pf{i}(1,2:end));
    
    if min(abs(p-N))==0
      PP=[PP,i];  
        
    end
    
    
end

Pf(PP)=[];









