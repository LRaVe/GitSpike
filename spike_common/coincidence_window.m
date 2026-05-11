function window=coincidence_window(tmin,tmax,spikes,spike_ind1,spike_ind2,ind1,ind2)
    % Calculate the coincidence window for two specific spikes
    
    n=length(spikes);
    if spike_ind1>n || spike_ind2>n || spike_ind1<1 || spike_ind2<1
        error('Index out of bounds');
    end

    s1=spikes{spike_ind1};
    s2=spikes{spike_ind2};

    if ind1>length(s1) || ind1<1
        error('Index out of bounds');
    end
    if ind2>length(s2) || ind2<1
        error('Index out of bounds');
    end

    if ind1==1
        prev1=s1(ind1)-tmin;
    else
        prev1=s1(ind1)-s1(ind1-1);
    end

    if ind1==length(s1)
        next1=tmax-s1(ind1);
    else
        next1=s1(ind1+1)-s1(ind1);
    end

    if ind2==1
        prev2=s2(ind2)-tmin;
    else
        prev2=s2(ind2)-s2(ind2-1);
    end

    if ind2==length(s2)
        next2=tmax-s2(ind2);
    else
        next2=s2(ind2+1)-s2(ind2);
    end

    % Coincidence window is half the minimum of all inter-spike distances
    window=min([prev1,next1,prev2,next2])/2;
end