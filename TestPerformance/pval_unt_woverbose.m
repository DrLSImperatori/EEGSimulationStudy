function [p_real, s_real, z_real] = pval_unt_woverbose(real_res,dummy_res,tail)
%% [p_real, s_real, z_real] = pval_unt(real_res,dummy_res,tail)
% Giulio Bernardi [giulioberna@gmail.com], 2018.05.20
% 2017.11.30 - Added z-score calculation
% 2018.05.20 - Added options for p-value calculation

% disp('Calculation of significance levels...');

null_mean=nanmean(dummy_res,2); % mean of null distribution
null_stdv=nanstd(dummy_res,[],2); % sample standard deviation of null distribution
p_real=NaN(size(real_res)); s_real=NaN(size(real_res)); z_real=NaN(size(real_res));
for i=1:size(real_res,1)
    z_real(i)=(real_res(i,:)-null_mean(i,:))./null_stdv(i,:); % z-score calculation
    switch tail
        case 'right'
%             p_real(i) = (length(find(dummy_res(i,:) >= real_res(i)))) / (size(dummy_res,2));
%             p_real(i) = (length(find(dummy_res(i,:) >= real_res(i)))+1) / (size(dummy_res,2)+1);
            p_real(i) = length(find(dummy_res(i,:) >= real_res(i))) / size(dummy_res,2); 
            s_real(i) = 1 ; 
        case 'left' 
%             p_real(i) = (length(find(dummy_res(i,:) <= real_res(i)))) / (size(dummy_res,2));
%             p_real(i) = (length(find(dummy_res(i,:) <= real_res(i)))+1) / (size(dummy_res,2)+1);
            p_real(i) = length(find(dummy_res(i,:) <= real_res(i))) / size(dummy_res,2); 
            s_real(i) = -1; 
        case 'both'
            if real_res(i) >= median(dummy_res(i,:))
%                 p_real(i) = (length(find(abs(dummy_res(i,:)) >= abs(real_res(i))))) / (size(dummy_res,2));
                 p_real(i) = (length(find(abs(dummy_res(i,:)) >= abs(real_res(i))))+1) / (size(dummy_res,2)+1);
%                 p_real(i) = (length(find(dummy_res(i,:) >= real_res(i))) / size(dummy_res,2)).*2;
                s_real(i) = 1;
            elseif real_res(i) <= median(dummy_res(i,:))
%                 p_real(i) = (length(find(abs(dummy_res(i,:)) >= abs(real_res(i))))) / (size(dummy_res,2));
                 p_real(i) = (length(find(abs(dummy_res(i,:)) >= abs(real_res(i))))+1) / (size(dummy_res,2)+1);
%                 p_real(i) = (length(find(dummy_res(i,:) <= real_res(i))) / size(dummy_res,2)).*2;
                s_real(i) = -1;
            end;
    end;
    
%     if p_real(i)==0 % adjust in case p-value is zero (fix p just below threshold)
%         p_real(i)=1/(size(dummy_res,2)+1); 
%         warning(['P-value smaller than perm-threshold']);
%     end; 
    
%     if p_real(i)>1; disp(num2str(p_real(i))); p_real(i)=1; % adjust in case p-value is larger than one
%         warning(['P-value larger than 1: adjusted']);
%     end; 
    
end; clear i;

end