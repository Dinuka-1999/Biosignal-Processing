function peaks=findPeaks(data,Limit)
    locs=find(data>Limit);
    %values = data(locs);
    
    actual_locs=[];
    
    temp_locs=[];
    
    max_loc=locs(1);
    j=1;
    i=1;
    for r=1:length(locs)
        if locs(r)-max_loc>100
            if isempty(temp_locs)
                temp_locs(j)=locs(r);
                j=j+1;
                max_loc=locs(r);
            else
                temp_values=data(temp_locs);
                [x,y]=max(temp_values);
                actual_locs(i)=temp_locs(y);
                i=i+1;
                temp_locs=[];
                j=1;
                temp_locs(j)=locs(r);
                j=j+1;
                max_loc=locs(r);
            end
            
        elseif ~isempty(temp_locs)
            temp_locs(j)=locs(r);
            j=j+1;
                
        end
    end
    peaks=actual_locs;
end