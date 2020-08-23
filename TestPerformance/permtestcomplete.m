function p_value=permtestcomplete(data1,data2,type,tail, nPerm)
tf = isequal(data1,data2);
if tf==1
    p_value=1;
else 
% disp('Calculation of real and null data values...');
data=cat(1,data1,data2); nsub1=size(data1,1); nsub2=size(data2,1);
switch type
    case 'Paired' % Hypothesis: [data1data2] comes from a distribution with zero median
        if nsub1~=nsub2; error('Number of subjects must be the same for a paired test!'); end;
        real_res=(nansum(data1-data2,1))';
    case 'Unpaired' % Hypothesis: [data1] and [data2] contain samples from continuous distributions with equal medians
        real_res=(nanmean(data1,1)-nanmean(data2,1))';
end
dummy_res=NaN(size(real_res,1),nPerm);
for i=1:nPerm  
%     if mod(i, 100)==0; disp(['  ',num2str(i),' out of ',num2str(nPerm),' permutations completed']); end
    switch type
        case 'Paired'
            signs=randi(2,nsub1,1); signs(signs==2)=-1; % Swapping conditions in paired test simply reverses the sign of the difference
            dummy_res(:,i)=nansum((data1-data2).*repmat(signs,[1 size(data1,2)]),1); 
        case 'Unpaired'
            rand_data=data(randperm(nsub1+nsub2),:); % Rearrangement of subjects in the two examined groups
            rand_data1=rand_data(1:nsub1,:);
            rand_data2=rand_data(nsub1+1:end,:);
            dummy_res(:,i)=nanmean(rand_data1,1)-nanmean(rand_data2,1);
    end
end; clear i;

[p_value, ~, ~] = pval_unt_woverbose(real_res,dummy_res,tail);
end
end
