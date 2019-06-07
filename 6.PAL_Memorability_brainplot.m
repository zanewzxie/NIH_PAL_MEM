% print(gcf,'scaled_decodingacc','-dpdf','-bestfit');
%% PLOT ALL Electrodes
close all

bp = brainplotter()
bd = braindata2()
bd.loadAverage()
bp =bd.ezplot(bp)
bp.view('ventral')

datapath='/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/Scripts/3_Mono/data/ATL_05_19';

for isub = allsub
    isub
     savefilename=fullfile(datapath,['SubjResults_ATL' num2str(isub) '.mat']);
    
    
    if exist(savefilename)>0
        
        load(savefilename);

        % identify and limit the electrodes for analysis, only look at the surface ALT electrodes.
        subj     = ['NIH0' num2str(isub)];
        
        % load subject braindata
        bd = braindata2(subj,rootEEGdir)
        
        % load each electrodes mesh info on the standardized brain
        t  = bd.roi.lead_mesh;
        
        clear elecs
        % get example electrodes
        elecs=SubjResults.allcleanChans;
        
        
        % Electrodes to exclude for misalignment
        % NIH026 G27
        % NIH029 G17
        % NIH032 ROF3
        % NIH036 OF3, OF4
        % NIH062 OF4
        % NIH066 TG84, TG123
        
        %         if strcmp(subj,'NIH026')
        %             elecs=elecs(~ismember(elecs,{'G27'}));
        %         elseif strcmp(subj,'NIH029')
        %             elecs=elecs(~ismember(elecs,{'G17'}));
        %         elseif strcmp(subj,'NIH032')
        %             elecs=elecs(~ismember(elecs,{'ROF3'}));
        %         elseif strcmp(subj,'NIH036')
        %             elecs=elecs(~ismember(elecs,{'OF3','OF4'}));
        %         elseif strcmp(subj,'NIH062')
        %             elecs=elecs(~ismember(elecs,{'OF4'}));
        %         elseif strcmp(subj,'NIH066')
        %             elecs=elecs(~ismember(elecs,{'TG84', 'TG123'}));
        %         end
        
        
        for j = 1:length(elecs)
            % get channels index within lead mesh table
            chan = strcmp(t.chanName,elecs{j});
            
            % get surface to plot on (left or right hemisphere)
            s = ['pial_' t.whichHemi{chan}];
            bp.plotPoint(t(chan,:),...
                'radius',1,...
                'mesh_clust_d',0,...
                'alpha',0.5,...
                'color',[1 0 0 ],...
                'surf',s) %'label', elecs{j}
            
            %   mesh_clust_d - If you are passing mesh index lists, by default a point is placed at the first mesh index in the list. However, for sulcus-
            %             spanning electrodes, we may want to plot 2+ points at each group of mesh indices. Pass a distance here if you want to clusterize
            %             such that mesh_indices >=~ mesh_clust_d apart will be plotted as separate points. Default is 3. 0 to just use the first index.
            % ,...
            
            
        end
        
        %keyboard
    end
    
end



% Electrodes to exclude for misalignment
% NIH026 G27
% NIH029 G17
% NIH032 ROF3
% NIH036 OF3, OF4
% NIH062 OF4
% NIH066 TG84, TG123


