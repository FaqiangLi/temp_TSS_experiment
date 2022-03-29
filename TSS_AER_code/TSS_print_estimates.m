function TSS_print_estimates(theta,se,inp,inp2)

frmlabels = inp2.frmcl9; catlabels = inp2.catcl;
stchlabels = inp2.stchcl; L1 = size(stchlabels,1);
hhchlabels = inp2.hhchcl; L2 = size(hhchlabels,1);
rnlabels = inp2.rncl; Lr = size(rnlabels,1);
splabels = inp2.spcl; Lsp = size(splabels,1);
cat_comb = inp.cat_comb; 
D = length(theta); 
K = inp.K; S = inp.S;

FileName = ['Output/' 'ParameterEstimates','_',datestr(now, 'ddmm_yyyy_HHMM'),'.xlsx'];
xlswrite(FileName,[theta se],1,['E5:F5' num2str(D+4)]) 
headings = {'Estimates' 'Standard errors'};
xlswrite(FileName,headings,1,'E4:F4') 
L0 = 5; % startline for writing
heading = {'Scaling parameter'};
xlswrite(FileName,heading,1,['B' num2str(L0)])
xlswrite(FileName,catlabels(1:K-1,1),1,['C' num2str(L0) ':C' num2str(L0+K-2)])
heading = {'Constant'};
xlswrite(FileName,heading,1,['B' num2str(L0+K-1)])
L0 = L0 + K-1;
xlswrite(FileName,stchlabels,1,['B' num2str(L0+1) ':C' num2str(L0+L1)])
L0 = L0 + L1;
xlswrite(FileName,hhchlabels,1,['B' num2str(L0+1) ':C' num2str(L0+L2)])
L0 = L0 + L2;
xlswrite(FileName,rnlabels,1,['B' num2str(L0+1) ':C' num2str(L0+Lr)])
L0 = L0 + Lr;
heading = {'Second-order terms'};
xlswrite(FileName,heading,1,['B' num2str(L0+1)])
xlswrite(FileName,catlabels(1:K,1),1,['C' num2str(L0+1) ':C' num2str(L0+K)])

heading = {'Second-order interaction terms'};
xlswrite(FileName,heading,1,['B' num2str(L0+K+1)])
l = 1;
for t=1:length(cat_comb)
    c_cat_t = cat_comb{t};
    no_cat = size(c_cat_t,2);
    if no_cat>1
        combs = combnk(c_cat_t,2);
        CC = size(combs,1);
        for cc=1:CC
            heading = [catlabels(combs(cc,1),1) catlabels(combs(cc,2),1) ];
            xlswrite(FileName,heading,1,['C' num2str(L0+K+l) ':D' num2str(L0+K+l) ])
            l = l + 1;
        end
    end
end
L0 = L0 + K + l - 1;
heading = {'Price coefficients'};
xlswrite(FileName,heading,1,['B' num2str(L0+1)])
heading = {'Price'};
xlswrite(FileName,heading,1,['C' num2str(L0+1)])
heading = {'Price/per capita income'};
xlswrite(FileName,heading,1,['C' num2str(L0+2)])
L0 = L0+2;
xlswrite(FileName,splabels,1,['B' num2str(L0+1) ':C' num2str(L0+Lsp)])

L0 = 5 + D - K*(S-1); 
heading = {'Firm-category effects'};
xlswrite(FileName,heading,1,[ 'B' num2str(L0)])
for k=1:K
    xlswrite(FileName,catlabels(k,1),1,['C' num2str(L0+(k-1)*(S-1))]) 
    xlswrite(FileName,frmlabels(2:end,1),1,...
        ['D' num2str(L0+(k-1)*(S-1)) ':D' num2str(L0+k*(S-1)-1)]) 
end

