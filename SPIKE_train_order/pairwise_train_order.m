function [res1,res2]=pairwise_train_order(tmin,tmax,spikes,spike_ind1,spike_ind2)
    % Compute bidirectional pairwise spike-train ordering
    % For each coincidence, both trains receive the same ordering value

    n=length(spikes);
    if spike_ind1>n || spike_ind2>n || spike_ind1<1 || spike_ind2<1
        error('Index out of bounds');
    end

    s1=spikes{spike_ind1};
    s2=spikes{spike_ind2};
    res1=zeros(1,length(s1));
    res2=zeros(1,length(s2));

    for i=1:length(s1)
        for j=1:length(s2)
            if abs(s1(i)-s2(j))<coincidence_window(tmin,tmax,spikes,spike_ind1,spike_ind2,i,j)
                % Both trains get the same ordering sign based on relative timing
                if s1(i)<s2(j)
                    signValue=1;
                elseif s1(i)>s2(j)
                    signValue=-1;
                else
                    signValue=0;
                end
                res1(i)=signValue;
                res2(j)=signValue;
                break;
            end
        end
    end
end