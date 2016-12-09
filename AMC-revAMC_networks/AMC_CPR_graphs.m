% Generate network diagrams for AMC and CPR landscapes 
% Requires Matlab's Bioinformatics Toolbox

fig = figure;
maxGrowth = 1.914; % for standardizing subplots
maxFreq = 0.2258; 

drugs = {'AMC','CPR','X1'};
genotypes = {'0000','1000','0100','0010','0001','1100','1010','1001','0110','0101','0011','1110','1101','1011','0111','1111'};

d=3; % select drug
%subplot(1,2,2)

% Import transition probabilities and make adjacency matrices
drugTrans = zeros(length(genotypes)+1,length(genotypes)+1);
drugAs = zeros(length(genotypes),length(genotypes));
filename=strcat(drugs{d},'_transition_probs.csv');
drugTrans=importTrans(filename);
drugAs = drugTrans(2:length(genotypes)+1,2:length(genotypes)+1);

% Import genotype frequencies
filename=strcat(drugs{d},'_adaptive_seq_geno_frequencies.csv');
freqs = importFreqs(filename,2,2);

% Import growth rates
filename=strcat(drugs{d},'_growth_rate_table.csv');
growth_rates = importGR(filename);

% Plot graph 
G = digraph(drugAs,genotypes,'OmitSelfLoops');
p=plot(G,'Layout','Force');
%p=plot(G);
axis off
axis square
% title(drugs(d))

% Adjust overlapping edges
% with bespoke coordinates
p.XData=[-2.00037818576126,-0.0172192295146225,-1.45263135870231,-2.53758312824889,-0.390000000000000,0.776463233211535,-0.370000000000000,1.76975269845423,-1.76975269852666,0.650000000000000,-0.776463233181368,0.274498759539414,2.53758312836726,1.45263135879342,0.180000000000000,2.00037818558011];
p.YData=[-1.91982690277580,-2.53263360374787,-0.239087533806790,0.160000000000000,-1.43406574276901,-0.819207999598380,-0.400000000000000,-1.82774612849780,1.82774612820811,0.600000000000000,0.819207999419269,1.43406574268422,0.0700000000000000,0.239087533708053,2.50000000000000,1.91982690269715];

% Weight edges
% via darkness
%p.EdgeCData = G.Edges.Weight;
%colormap(gray) % requires two colormaps, crap

% via thickness
G.Edges.LWidths=7*G.Edges.Weight/max(G.Edges.Weight);
p.LineWidth = G.Edges.LWidths;
p.EdgeColor = 'black';

% Color nodes by fitness
p.NodeCData = growth_rates;
colormap(viridis())
caxis([0 maxGrowth])
%if d==2
%    subplot(1,3,3)
    axis off
    cb = colorbar('southoutside');
    cmap = colormap(viridis());
    caxis([0 maxGrowth])
    cpos = cb.Position;
    currentHeight = cpos(4);
    currentY = cpos(2);
    cpos(4)= cpos(4)/3; % changes height
    newY = currentY + currentHeight/2 - cpos(4)/2;
    cpos(2)=newY;
    cb.Position=cpos;
    %ylabel(cb,'Growth rate (','Font','Times','FontSize',9);
    cb.Label.String='Growth rate (\times 10^{-3})';
    cb.Label.FontSize=9;
    cb.Label.FontName='Times';
    cb.FontSize=8;
    cb.FontName='Times';
    %ax.Position = axpos;
%end

%Change node size to reflect frequency
for n = 1:length(freqs)
    highlight(p,n,'MarkerSize',3+20*freqs(n)/maxFreq);
end

% Highlight edge between 1101 and 0101
% if d==1
%     highlight(p,[10 13],'EdgeColor','r');
% elseif d==2
%     highlight(p,[13 10],'EdgeColor','r');
% end

saveas(gcf,'AMC_demo.fig','fig')
fig.PaperUnits = 'inches';
fig.PaperSize = 1.2*[3.75 5];
fig.PaperPosition = 1.2*[0,0,3.75,5];
set(gcf,'renderer','painters')
print(gcf,'AMC_demo','-dpdf')

saveas(gcf,'AMC_unforced.fig','fig')
fig.PaperUnits = 'inches';
fig.PaperSize = 1.5*[4 5];
fig.PaperPosition = 1.5*[0,0,4,5];
print(gcf,'AMC_unforced','-depsc')

saveas(gcf,'rAMC_freq.fig','fig')
fig.PaperUnits = 'inches';
fig.PaperSize = 1.2*[3.75 5];
fig.PaperPosition = 1.2*[0,0,3.75,5];
set(gcf,'renderer','painters')
print(gcf,'rAMC_freq','-dpdf')
