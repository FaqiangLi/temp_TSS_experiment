function TSS_print_predictions(theta,inp,inp2)

frmlabels = inp2.frmcl; catlabels = inp2.catcl;
frmlabels9 = inp2.frmcl9;

warning('off','MATLAB:xlswrite:AddSheet'); % drop matlab warning
predictions = TSS_predictions(theta,inp);
[obsvis,vis,obsrev,revenue,table] = TSS_pairs(theta,inp);

TSS_distance_histogram(theta,inp);


% sheet 1
K = inp.K; S = inp.S;

FileName = ['Output/' 'Predictions','_',datestr(now, 'ddmm_yyyy_HHMM'),'.xlsx'];
xlswrite(FileName,predictions,1,'D6:U77') 
headings = {'Quantities' '	'	'Visitors'	' '	'Quantities, 1SS' ' ' 'Visitors, 1SS'	' '	'Quantities, 2SS'	' '	'Visitors, 2SS'	' '	'Average income' ' ' 'Avg. Hhsize'	' '	'Avg. dist. Travelled' ' ';	
'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred'	'Obs' 	'Pred'	'Obs' 	'Pred'	'Obs' 	'Pred'};
xlswrite(FileName,headings,1,'D4:U5') 
for k=1:K
    xlswrite(FileName,catlabels(k,1),1,['B' num2str(6+(k-1)*S)]) 
    xlswrite(FileName,frmlabels9,1,...
        ['C' num2str(6+(k-1)*S) ':C' num2str(6+k*S-1)]) 
end

% sheet 2
heading = {' ' ' ' 'Revenue, % of total' ' ' 'visits, % of total no. of households*' ' ' 'Average household size' ' ' 'Average household income' ' ' 'Average distance travelled to stores' ' ' ' '	'*) adds to more than 100 percent since 2-firm visits count doubly.';
' ' ' '		'1-firm sh.'	'2-firm sh.'	'1-firm sh.'	'2-firm sh.'	'1-firm sh.'	'2-firm sh.'	'1-firm sh.'	'2-firm sh.'	'1-firm sh.'	'2-firm sh.'	    ' ' 'Note: here the distinction is 1-firm/2-firm shoppers, not 1-store/2-store shoppers.'};
xlswrite(FileName,heading,2,'B5:O6')
xlswrite(FileName,table,2,'D7:M38')
heading = {'observed'; 'predicted'};
for s = 1:16
    xlswrite(FileName,frmlabels(s,1),2,['B' num2str(7+2*(s-1)) ])
    xlswrite(FileName,heading,2,['C' num2str(7+2*(s-1))...
        ':C' num2str(7+2*(s-1)+1)])
end

% sheet 3
xlswrite(FileName,frmlabels',3,'E2:T2')
xlswrite(FileName,frmlabels,3,'D3:D18')
heading = {'Visitors by firm combination (observed)'};
xlswrite(FileName,heading,3,'C1')
xlswrite(FileName,obsvis,3,'E3:T18')
xlswrite(FileName,frmlabels',3,'E21:T21')
xlswrite(FileName,frmlabels,3,'D22:D37')
heading = {'Visitors by firm combination (predicted)'};
xlswrite(FileName,heading,3,'C20')
xlswrite(FileName,vis,3,'E22:T37')

% sheet 4
xlswrite(FileName,frmlabels',4,'E2:T2')
xlswrite(FileName,frmlabels,4,'D3:D18')
heading = {'Revenue by firm combination (observed)'};
xlswrite(FileName,heading,4,'C1')
xlswrite(FileName,obsrev,4,'E3:T18')
xlswrite(FileName,frmlabels',4,'E21:T21')
xlswrite(FileName,frmlabels,4,'D22:D37')
heading = {'Revenue by firm combination (predicted)'};
xlswrite(FileName,heading,4,'C20')
xlswrite(FileName,revenue,4,'E22:T37')


% sheet 5
outofsample = true;
N = 4000;
[inp,~] = TSS_input(N,inp.cat_comb,inp.IV_comb,outofsample);
predictions = TSS_predictions(theta,inp);
xlswrite(FileName,predictions,5,'D6:U77') 
headings = {'Quantities' '	'	'Visitors'	' '	'Quantities, 1SS' ' ' 'Visitors, 1SS'	' '	'Quantities, 2SS'	' '	'Visitors, 2SS'	' '	'Average income' ' ' 'Avg. Hhsize'	' '	'Avg. dist. Travelled' ' ';	
'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred' 	'Obs' 	'Pred'	'Obs' 	'Pred'	'Obs' 	'Pred'	'Obs' 	'Pred'};
xlswrite(FileName,headings,5,'D4:U5') 
for k=1:K
    xlswrite(FileName,catlabels(k,1),5,['B' num2str(6+(k-1)*S)]) 
    xlswrite(FileName,frmlabels9,5,...
        ['C' num2str(6+(k-1)*S) ':C' num2str(6+k*S-1)]) 
end