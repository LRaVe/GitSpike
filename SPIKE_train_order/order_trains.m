function [results,order_matrix]=order_trains(tmin,tmax,spikes)
    % Compute spike-train-order values for all spikes and pairwise relationships
    % Aggregates pairwise train-order vectors for all spike trains

    n=length(spikes);
    results=cell(n,1);
    order_matrix=zeros(n,n);

    % Initialize result vectors for each train
    for i=1:n
        results{i}=zeros(1,length(spikes{i}));
        order_matrix(i,i)=0;  % Diagonal is always 0 (train compared to itself)
    end

    % Compute pairwise orderings between all distinct train pairs
    for i=1:n-1
        for j=i+1:n
            [res_i,res_j]=pairwise_train_order(tmin,tmax,spikes,i,j);
            results{i}=results{i}+res_i;
            results{j}=results{j}+res_j;
            order_matrix(i,j)=mean(res_i);
            order_matrix(j,i)=-mean(res_i);
        end
    end

    % Normalize results by number of other trains (n-1)
    if n>1
        for i=1:n
            results{i}=results{i}/(n-1);
        end
    end
end