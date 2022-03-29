function [mcost,pr_w] = TSS_print_elasticities_markups(theta,inp,R2,steplength)

warning('off','MATLAB:xlswrite:AddSheet'); % drop matlab warning

frmlabels = inp.frmcl; catlabels = inp.catcl;
K = inp.K;
S2 = 16; % number of different firms for profit max foc

FileName = ['Output/' 'Elasticities_markups','_',datestr(now, 'ddmm_yyyy_HHMM'),'.xlsx'];

% sheets 1-3
elast_labels = {'q,d,c all free.', 'q,d free; c fixed.', 'q free; d,c fixed.' };
fix_entry = {'empty', 'c', 'dc'};
%disp('Warning: standard elasticities not computed! only decomposed ones...')
median_elast = zeros(K,K,3,S2);
for r=1:3  % 'decomposed' elasticities, where some parts of consumer choice are held fixed
    xlswrite(FileName,{['number of draws for nu: ' num2str(R2)]},r,'M1')
    xlswrite(FileName,{['step length for demand derivatives: ' num2str(steplength)]},r,'M2')
    inp.fix = fix_entry{r};
    if r==1
        [markupsetc,elast] = TSS_multipledraws_markups(theta,inp,S2,R2,steplength);
        save eletc markupsetc elast
        load eletc
        mcost = markupsetc(:,2);
        pr_w = markupsetc(:,1); % weighted price
    elseif r>1
        [~,elast] = TSS_multipledraws_markups(theta,inp,S2,R2,steplength);
    end
    for s=1:S2
       median_elast(:,:,r,s) = elast((s-1)*K+1:s*K,(s-1)*K+1:s*K);
    end
    
    xlswrite(FileName,elast,r,'C5:DZ132')
    frmcat = cell(K*S2,1);
    for s=1:S2
        frmcat{(s-1)*K+1,1} = frmlabels{s,1};
        for k=1:K
            frmcat{(s-1)*K+k,2} = catlabels{k,1};
        end
    end
    xlswrite(FileName,frmcat',r,'C3:DZ4')
    xlswrite(FileName,frmcat,r,'A5:B132')
    heading = {'Price elasticities of category demands (quantities)';
        '(elasticity of row quantity wrt. column price)'; elast_labels{r}};
    xlswrite(FileName,heading,r,'A1:A2')
    heading = {elast_labels{r}};
    xlswrite(FileName,heading,r,'F2')
end
save median_elast median_elast
% sheet 4
xlswrite(FileName,{['number of draws for nu: ' num2str(R2)]},4,'M1')
xlswrite(FileName,{['steplength for demand derivatives: ' num2str(steplength)]},4,'M2')

heading = {'Implied marginal cost and markups, under supermarket pricing and delegation'};
xlswrite(FileName,heading,4,'A2')
heading = {'price' 'marg.cost' 'markup' '% markup' 'marg.cost deleg.' 'markup deleg.' '% markup deleg.'};
xlswrite(FileName,heading,4,'C4:I4')
xlswrite(FileName,frmcat,4,'A5:B132')
xlswrite(FileName,markupsetc,4,'C5:I132')

