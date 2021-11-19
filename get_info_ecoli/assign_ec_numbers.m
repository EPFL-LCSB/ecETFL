cd /Users/omid/Desktop/GITO/GECKO/geckomat/get_enzyme_data
% model_data = getEnzymeCodes(modelr,'ignore'); % model is the GEM of the organism

% iJO1366ecnumbers is the table of gene names and their associated EC
% numbers
EC_nembers = table2cell(iJO1366ecnumbers);
ind = ones(length(modelr.rxns),1);
gather = cell(length(modelr.rxns),20);
for i=1:size(EC_nembers,1)
    this_gene = EC_nembers{i,1};
    ind_rxns = find(contains(modelr.grRules,this_gene));
    for j=1:length(ind_rxns)
        gather{ind_rxns(j),ind(ind_rxns(j))} = char(EC_nembers{i,2});
    end
    ind(ind_rxns) = ind(ind_rxns) +1;
end

kcats = matchKcats(model_data, 'escherichia coli');

kcat_fwd = max(kcats.forw.kcats,[],2);
kcat_bwd = max(kcats.back.kcats,[],2);

