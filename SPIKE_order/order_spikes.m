function [results, SO_matrix]=order_spikes(tmin,tmax,spikes)
% Aggregate pairwise order vectors for all spike trains.
    n=length(spikes);
    results=cell(n,1);
    SO_matrix=zeros(n,n);

    if n==0
        return;
    end

    if n==1
        results{1}=zeros(1,length(spikes{1}));
        return;
    end

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

    if nargout>1
        for i=1:n-1
            for j=i+1:n
                so_profile_i=pairwise_order(tmin,tmax,spikes,i,j);
                so_profile_j=pairwise_order(tmin,tmax,spikes,j,i);
                num_pair_spikes=length(spikes{i})+length(spikes{j});

                if num_pair_spikes>0
                    so_sum=sum(so_profile_i)+sum(so_profile_j);
                    so_value=so_sum/num_pair_spikes;
                    SO_matrix(i,j)=so_value;
                    SO_matrix(j,i)=so_value;
                end
            end
        end
    end
end
