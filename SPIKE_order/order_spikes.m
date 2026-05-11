function results=order_spikes(tmin,tmax,spikes)
% Aggregate pairwise order vectors for all spike trains.
    n=length(spikes);
    results=cell(n,1); 

    for i=1:n
        aggregated=zeros(1,length(spikes{i}));
        
        % Sum pairwise orderings with all other spike trains
        for j=1:n
            if i~=j
                pairwise=pairwise_order(tmin,tmax,spikes,i,j);
                aggregated=aggregated+pairwise; 
            end
        end
        results{i}=aggregated/(n-1);
    end
end
