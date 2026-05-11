function res=pairwise_order(tmin,tmax,spikes,spike_ind1,spike_ind2)
    % Compute pairwise spike ordering between two spike trains

    n=length(spikes);
    if spike_ind1>n || spike_ind2>n || spike_ind1<1 || spike_ind2<1
        error('Index out of bounds');
    end

    s1=spikes{spike_ind1};
    s2=spikes{spike_ind2};
    res=zeros(1,length(s1));

    for i=1:length(s1)
        for j=1:length(s2)
            if abs(s1(i)-s2(j))<coincidence_window(tmin,tmax,spikes,spike_ind1,spike_ind2,i,j)
                if s2(j)>s1(i)
                    res(i)=1;
                elseif s2(j)<s1(i)
                    res(i)=-1;
                else
                    res(i)=0;
                end
                break;  % Move to the next spike in s1 after finding a coincident spike in s2
            end
        end
    end
end


