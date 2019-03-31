    %%  ploting Probe
    % memory accuracy
    y1=squeeze(mean(MeanATLPower_Probe(AccuracyMask,:,tasktime)));
    y2 = squeeze(mean(MeanATLPower_Probe(~AccuracyMask,:,tasktime)));   
    y3=y1-y2;
    
    SubjResults(isub).Probe_AccurateMeanPower=y1;
    SubjResults(isub).Probe_InAccurateMeanPower=y2;
    SubjResults(isub).Probe_AccuracyMainEffect=y3;
    
    allmin= -0.4; %round(min(min(y3)),1)*0.5;
    allmax=-allmin;
    
    h=figure(101);clf
    set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
    
    subplot(2,3,1);
    contourf(timePercent(tasktime),waveletFreqs,y1,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax]);    xlabel('Probe Onset');ylabel('freq');
    
    title('Accurate Recall')
    
    subplot(2,3,2);
    contourf(timePercent(tasktime),waveletFreqs,y2,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax]);    xlabel('Probe Onset');ylabel('freq');
    
    title('Inaccurate Recall')
    
    subplot(2,3,3);
    contourf(timePercent(tasktime),waveletFreqs,y3,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax]);    xlabel('Probe Onset');ylabel('freq');
    
    title('Diff, Accurate-Inaccurate Recall')
    
    % memorability
    y4=squeeze(mean(MeanATLPower_Probe(HighMemMask,:,tasktime)));
    y5=squeeze(mean(MeanATLPower_Probe(~HighMemMask,:,tasktime)));
    
    y6=y4-y5;
    
    SubjResults(isub).Probe_HighMemoMeanPower=y4;
    SubjResults(isub).Probe_LowMemoMeanPower=y5;
    SubjResults(isub).Probe_MemorabiltyMainEffect=y6;
    
    subplot(2,3,4);
    contourf(timePercent(tasktime),waveletFreqs,y4,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    xlabel('Probe Onset');ylabel('freq');
    
    title('Higher Memorablilty')
    
    subplot(2,3,5);
    contourf(timePercent(tasktime),waveletFreqs,y5,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    xlabel('Probe Onset');ylabel('freq');
    title('Lower Memorablilty')
    
    subplot(2,3,6);
    contourf(timePercent(tasktime),waveletFreqs,y6,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    title('Diff Higher-Lower Memorablilty');
    xlabel('Probe Onset');ylabel('freq');
    printfilename=fullfile(outpath,['Probe_Subj_' num2str(isub) '_MainEffects.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
    
    % Interaction Effect
    y7=squeeze(mean(MeanATLPower_Probe(HighMemMask & AccuracyMask,:,tasktime),1));
    y8=squeeze(mean(MeanATLPower_Probe(HighMemMask & ~AccuracyMask,:,tasktime),1));
    y9=squeeze(mean(MeanATLPower_Probe(~HighMemMask & AccuracyMask,:,tasktime),1));
    y10=squeeze(mean(MeanATLPower_Probe(~HighMemMask & ~AccuracyMask,:,tasktime),1));
    y11= (y7-y8)-(y9-y10);
    
    SubjResults(isub).Probe_InteractionEff=y11;
    
    h=figure(102);clf
    contourf(timePercent(tasktime),waveletFreqs,y11,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    title('Interaction Higher-Lower Memorablilty (Accurate-Inaccurate)')
    xlabel('Probe Onset');ylabel('freq');
    printfilename=fullfile(outpath,['Probe_Subj_' num2str(isub) '_InteractionEffect.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
    
    
    %%  ploting Recall
    % memory accuracy
    y1=squeeze(mean(MeanATLPower_Recall(AccuracyMask,:,tasktime)));
    y2=squeeze(mean(MeanATLPower_Recall(~AccuracyMask,:,tasktime)));   
    y3=y1-y2;
    
    SubjResults(isub).Recall_AccurateMeanPower=y1;
    SubjResults(isub).Recall_InAccurateMeanPower=y2;
    SubjResults(isub).Recall_AccuracyMainEffect=y3;
    
    allmin= -0.4; %round(min(min(y3)),1)*0.5;
    allmax=-allmin;
    
    h=figure(101);clf
    set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
    
    subplot(2,3,1);
    contourf(timePercent(tasktime),waveletFreqs,y1,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax]);    xlabel('Recall Onset');ylabel('freq');
    
    title('Accurate Recall')
    
    subplot(2,3,2);
    contourf(timePercent(tasktime),waveletFreqs,y2,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax]);    xlabel('Recall Onset');ylabel('freq');
    
    title('Inaccurate Recall')
    
    subplot(2,3,3);
    contourf(timePercent(tasktime),waveletFreqs,y3,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax]);    xlabel('Recall Onset');ylabel('freq');
    
    title('Diff, Accurate-Inaccurate Recall')
    
    % memorability
    y4=squeeze(mean(MeanATLPower_Recall(HighMemMask,:,tasktime)));
    y5=squeeze(mean(MeanATLPower_Recall(~HighMemMask,:,tasktime)));
    
    y6=y4-y5;
    
    SubjResults(isub).Recall_HighMemoMeanPower=y4;
    SubjResults(isub).Recall_LowMemoMeanPower=y5;
    SubjResults(isub).Recall_MemorabiltyMainEffect=y6;
    
    subplot(2,3,4);
    contourf(timePercent(tasktime),waveletFreqs,y4,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    xlabel('Recall Onset');ylabel('freq');
    
    title('Higher Memorablilty')
    
    subplot(2,3,5);
    contourf(timePercent(tasktime),waveletFreqs,y5,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    xlabel('Recall Onset');ylabel('freq');
    title('Lower Memorablilty')
    
    subplot(2,3,6);
    contourf(timePercent(tasktime),waveletFreqs,y6,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    title('Diff Higher-Lower Memorablilty');
    xlabel('Recall Onset');ylabel('freq');
    printfilename=fullfile(outpath,['Recall_Subj_' num2str(isub) '_MainEffects.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
    
    % Interaction Effect
    y7=squeeze(mean(MeanATLPower_Recall(HighMemMask & AccuracyMask,:,tasktime),1));
    y8=squeeze(mean(MeanATLPower_Recall(HighMemMask & ~AccuracyMask,:,tasktime),1));
    y9=squeeze(mean(MeanATLPower_Recall(~HighMemMask & AccuracyMask,:,tasktime),1));
    y10=squeeze(mean(MeanATLPower_Recall(~HighMemMask & ~AccuracyMask,:,tasktime),1));
    y11= (y7-y8)-(y9-y10);
    
    SubjResults(isub).Recall_InteractionEff=y11;
    
    h=figure(102);clf
    contourf(timePercent(tasktime),waveletFreqs,y11,200,'linecolor','none');
    set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
    colorbar;caxis([allmin allmax])
    title('Interaction Higher-Lower Memorablilty (Accurate-Inaccurate)')
    xlabel('Recall Onset');ylabel('freq');
    printfilename=fullfile(outpath,['Recall_Subj_' num2str(isub) '_InteractionEffect.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
    
