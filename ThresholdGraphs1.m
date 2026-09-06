%% STEP ONE %%%%%%%%%%%%%%% IMPORT DATA %%%%%%%%%%%%%%%%%%%%
% This script is the statistics-and-visualization stage of a larger
% pipeline. It expects a handful of preprocessed structures to already
% exist in the MATLAB workspace before it is run:
%   AW      - one struct per subject/scan, raw computed network metrics
%             and demographic fields (age, sex, diagnosis group, etc.)
%   As, Es  - AW reorganized into struct-of-arrays form, one field per
%             metric/subnetwork combination (e.g. "IFunctionalThresh...")
%   Efields - fieldnames(Es), the list of metric/subnetwork combinations
%             the main loop below iterates over
% The earlier data-preparation stage that builds these from SAOM results 
% and network analyses, and the underlying participant dataset itself, are not
% included in this repository -- the dataset involves clinical/
% neuroimaging data from human subjects and is not publicly shareable.
% This repo exists to show the statistical-analysis and figure-generation
% approach, not to be run end-to-end by someone without access to that
% data.

%% STEP TWO %%%%%%%%%%%%%% USER SETTINGS %%%%%%%%%%%%%%%%%%%%
% The toggles reflect the analysis setting that can be chosen. 
% Turning them to false omits those particular analyses or graphs.

% Saving Locations %
% Where graphs and result tables get saved. Actual plot files land under
% one of three subfolders (see resolveSavePath): "Individual Plots" and
% "Grid Plots and Heat Maps" each further split into
% <Functional|Structural>\<Subnet, blank for whole-brain>\; "Aggregate Plots"

Settings.outputDir = "./Output/"; % relative to the working directory MATLAB is launched from

% DATA TOGGLES %
Settings.runFunctionalWholeBrain = true;
Settings.runFunctionalSubnets    = true;
Settings.runStructuralWholeBrain = true;
Settings.runStructuralSubnets    = true;

% OUTLIER SETTING %
% isoutlier ThresholdFactor used when removing outliers from each subnetwork's data.
Settings.outlierThreshold = 3;

% ANALYSIS TOGGLES %

% Disorder Comparisons 
Settings.tests.Disorder.includeQuadratic  = false;
Settings.tests.Disorder.doAgeTrendPlot    = true;
Settings.tests.Disorder.doViolinPlot      = true;
Settings.tests.Disorder.doGridPlot        = true;

% MCI/Control Comparisons
Settings.tests.MCI.includeQuadratic  = false;
Settings.tests.MCI.doAgeTrendPlot    = true;
Settings.tests.MCI.doViolinPlot      = true;
Settings.tests.MCI.doGridPlot        = true; 

%AD/Control Comparisons
Settings.tests.AD.includeQuadratic  = false;
Settings.tests.AD.doAgeTrendPlot    = true;
Settings.tests.AD.doViolinPlot      = true;
Settings.tests.AD.doGridPlot        = true; % added 2026-09-04 -- was missing while Disorder/Sex/SexxDisorder already had it

% Sex Comparisons
Settings.tests.Sex.includeQuadratic  = false;
Settings.tests.Sex.doAgeTrendPlot    = true;
Settings.tests.Sex.doViolinPlot      = true;
Settings.tests.Sex.doGridPlot        = true;

%Sex and Disorder Anovas
Settings.tests.SexxDisorder.includeQuadratic  = false;
Settings.tests.SexxDisorder.doAgeTrendPlot    = true;
Settings.tests.SexxDisorder.doViolinPlot      = true;
Settings.tests.SexxDisorder.doGridPlot        = true;

%Sex and MCI Anovas
Settings.tests.SexxMCI.includeQuadratic  = false;
Settings.tests.SexxMCI.doAgeTrendPlot    = true;
Settings.tests.SexxMCI.doViolinPlot      = true;
Settings.tests.SexxMCI.doGridPlot        = false; 

%Sex and AD Anovas
Settings.tests.SexxAD.includeQuadratic  = false;
Settings.tests.SexxAD.doAgeTrendPlot    = true;
Settings.tests.SexxAD.doViolinPlot      = true;
Settings.tests.SexxAD.doGridPlot        = false; 

% Impairment Severity Anova
Settings.tests.Severity.doViolinPlot    = true;
Settings.tests.Severity.doGridPlot      = true; 

% Overall Age Analysis
Settings.tests.Overall.doAgeTrendPlot   = true; 

% Aggregate Figures 
Settings.aggregateFigures.mainFigures   = true;
Settings.aggregateFigures.subnetFigures = true;
Settings.aggregateFigures.numHeadlineScatters = 3; % how many Figure 3b candidates to flag

% TOGGLE OVERRIDE %

% Toggle overrides for resolveToggle(): exceptions only, most-specific
% wins. Empty until Phase 1's toggle migration actually wires call sites
% through resolveToggle() -- see resolveToggle's own comments for the
% dotted-path/override convention. Field names use "_" in place of "."
% (e.g. Settings.toggles.Individual_TTest_Impairment_Functional = true).
Settings.toggles = struct();

if ~isfolder(Settings.outputDir)
    mkdir(Settings.outputDir);
end


% GRAPHICAL SETTINGS %

% Establishing some settings for the graphing including the names of the factors. 
% Future work may expand or change this list of factors. 
load("Graphing.mat") % color palettes/markers/group-code lookup arrays designed by me
longnames2=[
            "Degree"   
            "Transitive Triads"
            "Betweenness"
            "Distance 2" 
            "4 Cycles"
            "Assortativity"
            "Euclidean Dist."
            "Euclidean Neighbor Dist."
            "Same Subnet"
            "Same Subnet Trans. Triplets"
            "Jump Subnet Trans. Triplets"
            "Same Subnet 4-Cycles"
            "Same Subnet Distance 2"
        ];

formulanames=[
            "_{D}"   
            "_{T}"
            "_{B}"
            "_{D2}" 
            "_{C4}"
            "_{A}"
            "_{ED}"
            "_{ND}"
            "_{S}"
            "_{ST}"
            "_{JT}"
            "_{S4}"
            "_{S2}"
        ];

MetricNames=["Rate";longnames2]; % one name per column of V/VSAS/etc., in the same order

% Subnetwork naming conventions
SubnetLongNames=["VN","SMN","DAN","SN","LN","FPN","DMN"];

% PULLING CORRECT DATA %

% This script comes at the end of serveral analyses here called the thresholded data

indeces=find(contains(Efields,'ThresholdedVS'));

indeces1=[find(contains(Efields,'fMRIThresholdedVS')),find(contains(Efields,'DTIThresholdedVS'))];

indeces=[indeces(ismember(indeces,indeces1));indeces(~ismember(indeces,indeces1))];

num=1;

% CREATING RESULTS FILES %

% Accumulates one row per (subnetwork, comparison, group, metric) with its
% weighted mean/SD, saved to a table at the end of the run.

MuSigmaResults=table('Size',[0,7],'VariableTypes',{'string','string','string','string','double','double','double'}, ...
    'VariableNames',{'Subnetwork','Comparison','Group','Metric','Mu','Sigma','N'});

% Accumulates one row per standardized effect (a main effect, an interaction,
% or a per-group age slope) -- the data Figures 1/2/3a read from. 'Group' is
% blank for a whole-comparison effect (e.g. the overall Impairment effect) and
% named for a per-group effect (e.g. the CN group's own age slope).

EffectsResults=table('Size',[0,12],'VariableTypes',{'string','string','string','string','string','double','double','double','double','double','double','string'}, ...
    'VariableNames',{'Subnetwork','Comparison','Term','Group','Metric','TStat','DF','Effect','CILow','CIHigh','PValue','Type'});

% CREATING RESULTS STRUCTS %

Datasets=struct();
Results=struct();

%% STEP THREE %%%%%%%%%%%%%%% DATA CLEANING %%%%%%%%%%%%%%%%%%%%
% For each threshold/subnetwork entry (one iteration per row of `indeces`),
% resolve scope, load raw values, filter/align them across every downstream
% array, and drop NaN/Inf/outlier/zero-error rows -- the shared numeric
% inputs (VSAS/ESAS/Diso/Diag/Sex/Age/Years/ID) that every group-subsetting,
% Mu/Sigma, and test step below is built from. The %% headers below mark
% each sub-step for the editor's next/previous-section navigation.

for i=indeces'

    %% Scope setup: Functional vs Structural, subnetwork name, output folder

    % GRAPH NAMES %

    F=strfind(Efields{i},"fMRI");
    D=strfind(Efields{i},"DTI");
    notskip=true;
    if ~isempty(F)
        set1="fMRI";
        SetName="Functional";
    elseif ~isempty(D)
        set1="DTI";
        SetName="Structural";
    else
        notskip=false;
    end


    if notskip

        % SUBNET SETTINGS %
        Subnetname=Efields{i}((strfind(Efields{i}(1:(strfind(Efields{i},"Thresh")-1)),"I")+1):(strfind(Efields{i},"Thresh")-1));
        

        if isempty(Subnetname)
            SubnetName="";
            TitleName1=SetName;
        else
            SubnetName=SubnetLongNames(find(FSnames==Subnetname));
            TitleName1=strcat(SetName," ",SubnetName);
        end

         % NODE SETTINGS %
        
        if ~isempty(Subnetname)
            Nodes=FSsizes(find(FSnames==Subnetname));
        else
            Nodes=100;
        end

        % FACTORS INFORMATION %

        eval(strcat("AllFactors=AW(1).",Efields{i}(1:(strfind(Efields{i},"Thresh")-1)),"AllFactors;"))
        

        if contains(TitleName1,"Structural")
            coe="S";
        else
            coe="F";
        end
        
        % DETERMINING TESTS TO RUN %

        % Skip this entry entirely if its category is switched off in Settings
        IsWholeBrain=isempty(Subnetname);
        IsFunctional=strcmp(SetName,"Functional");
        if IsFunctional && IsWholeBrain && ~Settings.runFunctionalWholeBrain
            continue
        end
        if IsFunctional && ~IsWholeBrain && ~Settings.runFunctionalSubnets
            continue
        end
        if ~IsFunctional && IsWholeBrain && ~Settings.runStructuralWholeBrain
            continue
        end
        if ~IsFunctional && ~IsWholeBrain && ~Settings.runStructuralSubnets
            continue
        end

       
        % TEST RESULT FILE %

        fid = fopen( strcat(Settings.outputDir,TitleName1,".txt"), 'wt' );
           
        %% Load raw factor values (rate + theta) for this threshold entry

        % JUST THRESHOLDED DATA %

        eval(strcat("ThisRate1=As(1).",Efields{i}(1:(strfind(Efields{i},"Thresh")-1)),set1,"Thresh","Rate';"))
        eval(strcat("ThisTheta1=Es(1).",Efields{i},"';"));
        eval(strcat("ThisRateES1=As(1).",Efields{i+1}(1:(strfind(Efields{i+1},"Thresh")-1)),set1,"Thresh","RateES';"))
        eval(strcat("ThisThetaES1=Es(1).",Efields{i+1},"';"));

        % Row count before any cleaning step -- the "before" half of the
        % new sample-size diagnostic (Datasets.(scope).nBefore/.nAfter).
        nBeforeCleaning=length(ThisRateES1);

        
        % DEMOGRAPHIC VARIABLES %
        
        DisoO=nanmean([[AW.CDRSB1];[AW.CDRSB2]])>=0.5;
        DiagO=((nanmean([[AW.CDRSB1];[AW.CDRSB2]])>=0.5)+(nanmean([[AW.CDRSB1];[AW.CDRSB2]])>=4));
        eval(strcat("SexO=Bs(1).Sex;"));
        eval(strcat("AgeO=As(1).Age2;"));
        eval(strcat("YearsO=As(1).Years;"));
        eval(strcat("IDO={AW.ID_subject};"))
        
        % REMOVE NA %

        ThisRateS=ThisRate1(~isnan(ThisRateES1));
        ThisThetaS=ThisTheta1(~isnan(ThisRateES1),:);
        ThisThetaESS=ThisThetaES1(~isnan(ThisRateES1),:);

        Diso=DisoO(~isnan(ThisRateES1));
        Diag=DiagO(~isnan(ThisRateES1));
        Sex=SexO(~isnan(ThisRateES1));
        Age=AgeO(~isnan(ThisRateES1));
        Years=YearsO(~isnan(ThisRateES1));
        ID=IDO(~isnan(ThisRateES1));
        ThisRateESS=ThisRateES1(~isnan(ThisRateES1));

        % CREATE VALUE AND ERROR ARRAYS %

        VSS=[ThisRateS,ThisThetaS];
        ESS=[ThisRateESS,ThisThetaESS];

        % CORRECT RATE FOR YEARS %

        VSS(:,1)=VSS(:,1).*Years';
        ESS(:,1)=ESS(:,1).*Years';
        
        % REMOVE ANY NA/INF AFTER YEARS CORRECTION %

        %If the years between is 0 then we will end up with infinitities
        
        These=(any(isnan(VSS'))+any(isnan(ESS')))>0;
        VSS(These,:)=[];
        ESS(These,:)=[];
        Diso(These)=[];
        Diag(These)=[];
        Sex(These)=[];
        Age(These)=[];
        Years(These)=[];
        ID(These)=[];

        Diso(find(any(isinf(VSS),2))) = [];
        Diag(find(any(isinf(VSS),2))) = [];
        Sex(find(any(isinf(VSS),2))) = [];
        Age(find(any(isinf(VSS),2))) = [];
        Years(find(any(isinf(VSS),2))) = [];
        ID(find(any(isinf(VSS),2))) = [];

        ESS(any(isinf(VSS),2),:) = [];
        VSS(any(isinf(VSS),2),:) = [];
        

        % OUTLIER REMOVAL %

        NotOutlier=~any(isoutlier([VSS(:,1),log(VSS(:,1)),VSS(:,(2:end)),ESS(:,1),ESS(:,(2:end))],ThresholdFactor=Settings.outlierThreshold),2);
        DisoA=Diso(NotOutlier);
        DiagA=Diag(NotOutlier);
        SexA=Sex(NotOutlier);
        AgeA=Age(NotOutlier);
        YearsA=Years(NotOutlier);
        IDA=ID(NotOutlier);
        ESAS=ESS(NotOutlier,:);
        VSAS=VSS(NotOutlier,:);


        % REMOVE ANYTHING WITH NO ERROR %

        ZeroErr=all(ESAS==0,2);
        DisoA(ZeroErr)=[];
        DiagA(ZeroErr)=[];
        SexA(ZeroErr)=[];
        AgeA(ZeroErr)=[];
        YearsA(ZeroErr)=[];
        IDA(ZeroErr)=[];
        VSAS(ZeroErr,:)=[];
        ESAS(ZeroErr,:)=[];

        % DATA TABLE %
        
        T=table(IDA(:),AgeA(:),YearsA(:),SexA(:),DisoA(:),DiagA(:), ...
            'VariableNames',{'ID','Age','Years','Sex','Disorder','Diagnosis'});
        T.Value=VSAS;
        T.Error=ESAS;

        % CALCULATE MU AND SIGMA
        [MuOverall,SigmaOverall]=weightedMuSigma(VSAS,ESAS);

        scopeField=matlab.lang.makeValidName(TitleName1);
        Datasets.(scopeField)=struct('Label',TitleName1,'T',T, ...
            'MuOverall',MuOverall,'SigmaOverall',SigmaOverall, ...
            'Nodes',Nodes,'SubnetName',SubnetName,'SetName',SetName, ...
            'nBefore',nBeforeCleaning,'nAfter',size(VSAS,1));

        % HISTOGRAMS %
        
        % raw vs. outlier-removed distribution per metric
        
        for f=1:size(VSAS,2)
            if ~all(VSAS(:,f)==0)
                fig=figure;
                fig.Visible='off';
                set(gca, 'box', 'off');
                subplot(1,2,1)
                hold on
                histfit(VSS(:,f))
                histogram(ESS(:,f))
                title("Original Distribution")
                subplot(1,2,2)
                hold on
                histfit(VSAS(:,f))
                histogram(ESAS(:,f))
                title("Original Distribution Remove Outliers")

                if f==1
                    callit=strcat(TitleName1,"/","Rate","_Histograms");
                    sgtitle(strcat(coe,"_\rho"))
                else
                    callit=strcat(TitleName1,"/",longnames2{f-1},"_Histograms");
                    sgtitle(strcat(coe,formulanames(f-1)))
                end
                set(gca, 'box', 'off');
                saveas(fig,strcat(Settings.outputDir,callit,".png"))
            end
        end

        % Side-by-side Functional-vs-Structural histogram (2026-09-06,
        % deferred): commented out until the Tests/Graphing reorg is
        % actually finished. This loop only ever runs with ONE modality's
        % data in scope (VSAS is whichever of Functional/Structural
        % Efields{i} currently is) -- VSAF/VSF/TitleName below were never
        % defined at this point, so as written this can't run. Once the
        % reorg gives us a point in the script where both modalities'
        % cleaned data exist together (the originally-planned Graphing
        % section), this belongs there instead, comparing that scope's
        % Datasets.(Functional...) and Datasets.(Structural...) tables.
        % for f=1:size(VSAS,2)
        %     if ~all(VSAS(:,f)==0) || ~all(VSAF(:,f)==0)
        %         fig=figure;
        %         fig.Visible='on';
        %         subplot(2,2,1)
        %         histfit(VSS(:,f))
        %         title("Structural Distribution")
        %         subplot(2,2,2)
        %         histfit(VSAS(:,f))
        %         title("Structural Distribution Remove Outliers")
        %         subplot(2,2,3)
        %         histfit(VSF(:,f))
        %         title("Functional Distribution")
        %         subplot(2,2,4)
        %         histfit(VSAF(:,f))
        %         title("Functional Distribution Remove Outliers")
        %
        %         if f==1
        %             callit=strcat(TitleName,"/",Subnetname,"/","Rate","_Histograms");
        %             sgtitle(strcat(TitleName," ","Rate"))
        %         else
        %             callit=strcat(TitleName,"/",Subnetname,"/",longnames2{f-1},"_Histograms");
        %             sgtitle(strcat(TitleName," ",longnames2(f-1)))
        %         end
        %         saveas(fig,strcat(Settings.outputDir,callit,".png"))
        %     end
        % end

        

        %%%%%%%%%%%%%%%%% ESTABLISH SUBGROUPS %%%%%%%%%%%%%%%%%%%%
        % Every comparison group used below, sliced from the same 7 cleaned
        % per-subject arrays (VSAS/ESAS/YearsA/IDA/AgeA/SexA/DisoA) by a
        % logical mask via subsetGroup -- see that function for why this
        % used to be 7 hand-written lines per group.

        
        %Control
        [VSS0,ESS0,Years0,ID0,Age0,Sex0,Diso0]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,DisoA==0);
        %Diagnosis
        [VSS1,ESS1,Years1,ID1,Age1,Sex1,Diso1]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,DisoA==1);
        %Male
        [VSS2,ESS2,Years2,ID2,Age2,Sex2,Diso2]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,SexA==0);
        %Female
        [VSS3,ESS3,Years3,ID3,Age3,Sex3,Diso3]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,SexA==1);
        %Male Control
        [VSS4,ESS4,Years4,ID4,Age4,Sex4,Diso4]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==0).*(DisoA==0))==1);
        %Male Diagnosed
        [VSS5,ESS5,Years5,ID5,Age5,Sex5,Diso5]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==0).*(DisoA==1))==1);
        %Female Control
        [VSS6,ESS6,Years6,ID6,Age6,Sex6,Diso6]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==1).*(DisoA==0))==1);
        %Female Diagnosed
        [VSS7,ESS7,Years7,ID7,Age7,Sex7,Diso7]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==1).*(DisoA==1))==1);
        %Control (shared MCI/AD)
        [VSS8,ESS8,Years8,ID8,Age8,Sex8,Diso8]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,DiagA==0);
        %MCI
        [VSS9,ESS9,Years9,ID9,Age9,Sex9,Diso9]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,DiagA==1);
        %AD
        [VSS10,ESS10,Years10,ID10,Age10,Sex10,Diso10]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,DiagA==2);
        %Male Control (shared MCI/AD)
        [VSS11,ESS11,Years11,ID11,Age11,Sex11,Diso11]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==0).*(DiagA==0))==1);
        %Male MCI
        [VSS12,ESS12,Years12,ID12,Age12,Sex12,Diso12]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==0).*(DiagA==1))==1);
        %Female Control (shared MCI/AD)
        [VSS13,ESS13,Years13,ID13,Age13,Sex13,Diso13]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==1).*(DiagA==0))==1);
        %Female MCI
        [VSS14,ESS14,Years14,ID14,Age14,Sex14,Diso14]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==1).*(DiagA==1))==1);
        %Male AD (Control group is shared with the MCI comparison: VSS11/ESS11)
        [VSS15,ESS15,Years15,ID15,Age15,Sex15,Diso15]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==0).*(DiagA==2))==1);
        %Female AD (Control group is shared with the MCI comparison: VSS13/ESS13)
        [VSS16,ESS16,Years16,ID16,Age16,Sex16,Diso16]=subsetGroup(VSAS,ESAS,YearsA,IDA,AgeA,SexA,DisoA,((SexA==1).*(DiagA==2))==1);

        %%%%%%%%%%%%%%%%% MU AND SIGMA CALCULATION %%%%%%%%%%%%%%%%%%%%

        [MuSigmaResults,MuA,SigmaA]=computeGroupMuSigma(VSAS,ESAS,TitleName1,"Overall","All",MetricNames,MuSigmaResults,"A",Efields,i,fid,"Overall Mu and Sigma");

        [MuSigmaResults,Mu0,Sigma0]=computeGroupMuSigma(VSS0,ESS0,TitleName1,"Disorder","Control",MetricNames,MuSigmaResults,"0",Efields,i,fid,"Control Mu and Sigma");
        [MuSigmaResults,Mu1,Sigma1]=computeGroupMuSigma(VSS1,ESS1,TitleName1,"Disorder","Diagnosed",MetricNames,MuSigmaResults,"1",Efields,i,fid,"Impaired Mu and Sigma");
        [MuSigmaResults,Mu2,Sigma2]=computeGroupMuSigma(VSS2,ESS2,TitleName1,"Sex","Male",MetricNames,MuSigmaResults,"2",Efields,i,fid,"Male Mu and Sigma");
        [MuSigmaResults,Mu3,Sigma3]=computeGroupMuSigma(VSS3,ESS3,TitleName1,"Sex","Female",MetricNames,MuSigmaResults,"3",Efields,i,fid,"Female Mu and Sigma");
        [MuSigmaResults,Mu4,Sigma4]=computeGroupMuSigma(VSS4,ESS4,TitleName1,"SexxDisorder","MaleControl",MetricNames,MuSigmaResults,"4",Efields,i,fid,"Male Control Mu and Sigma");
        [MuSigmaResults,Mu5,Sigma5]=computeGroupMuSigma(VSS5,ESS5,TitleName1,"SexxDisorder","MaleDiagnosed",MetricNames,MuSigmaResults,"5",Efields,i,fid,"Male Impaired Mu and Sigma");
        [MuSigmaResults,Mu6,Sigma6]=computeGroupMuSigma(VSS6,ESS6,TitleName1,"SexxDisorder","FemaleControl",MetricNames,MuSigmaResults,"6",Efields,i,fid,"Female Control Mu and Sigma");
        [MuSigmaResults,Mu7,Sigma7]=computeGroupMuSigma(VSS7,ESS7,TitleName1,"SexxDisorder","FemaleDiagnosed",MetricNames,MuSigmaResults,"7",Efields,i,fid,"Female Impaired Mu and Sigma");
        [MuSigmaResults,Mu8,Sigma8]=computeGroupMuSigma(VSS8,ESS8,TitleName1,"Diagnosis","Control",MetricNames,MuSigmaResults,"8",Efields,i,fid,""); % shared Control for MCI/AD -- original never logged this one
        [MuSigmaResults,Mu9,Sigma9]=computeGroupMuSigma(VSS9,ESS9,TitleName1,"Diagnosis","MCI",MetricNames,MuSigmaResults,"9",Efields,i,fid,"Mild Impairment Mu and Sigma");
        [MuSigmaResults,Mu10,Sigma10]=computeGroupMuSigma(VSS10,ESS10,TitleName1,"Diagnosis","AD",MetricNames,MuSigmaResults,"10",Efields,i,fid,"Significant Impairment Mu and Sigma");
        [MuSigmaResults,Mu11,Sigma11]=computeGroupMuSigma(VSS11,ESS11,TitleName1,"SexxMCI","MaleControl",MetricNames,MuSigmaResults,"11",Efields,i,fid,""); % shared Control for MCI/AD -- original never logged this one
        [MuSigmaResults,Mu12,Sigma12]=computeGroupMuSigma(VSS12,ESS12,TitleName1,"SexxMCI","MaleMCI",MetricNames,MuSigmaResults,"12",Efields,i,fid,"Male Mild Impairment Mu and Sigma");
        [MuSigmaResults,Mu13,Sigma13]=computeGroupMuSigma(VSS13,ESS13,TitleName1,"SexxMCI","FemaleControl",MetricNames,MuSigmaResults,"13",Efields,i,fid,""); % shared Control for MCI/AD -- original never logged this one
        [MuSigmaResults,Mu14,Sigma14]=computeGroupMuSigma(VSS14,ESS14,TitleName1,"SexxMCI","FemaleMCI",MetricNames,MuSigmaResults,"14",Efields,i,fid,"Female Mild Impairment Mu and Sigma");
        [MuSigmaResults,Mu15,Sigma15]=computeGroupMuSigma(VSS15,ESS15,TitleName1,"SexxAD","MaleAD",MetricNames,MuSigmaResults,"15",Efields,i,fid,"");
        [MuSigmaResults,Mu16,Sigma16]=computeGroupMuSigma(VSS16,ESS16,TitleName1,"SexxAD","FemaleAD",MetricNames,MuSigmaResults,"16",Efields,i,fid,"");


        %OVERALL WITH AGE
        if Settings.tests.Overall.doAgeTrendPlot
        for q=1
            if q==1
                V=VSAS;
                R=ESAS;
                A=AgeA;
                I=IDA;
                TitleName=TitleName1;
            else
                V=VSAF;
                R=ESAF;
                A=AgeA;
                I=IDA;
                TitleName=TitleName2;
            end
            for fac=1:size(V,2)
                if ~all(V(:,fac)==0)

                    V0=A(:);
                    V1=V(:,fac);
                    E1=R(:,fac);

                    
                    plus=8;
                    % if fac==1
                    %     Weights= ones(size(E1));
                    % else
                        Weights= inverseVarianceWeights(E1);
                    % end

                    data.ID=categorical(I(:));
                    data.AgeC=V0-mean(V0);
                    data.Age=V0;
                    data.Weights=Weights;
                    data.Score=V1;

                    myTable = struct2table(data);

                    try
                        lme=fitlme(myTable, 'Score ~ AgeC + (1|ID)','Weights',data.Weights);
                    catch ME
                        fprintf("Skipping Overall Age fac=%d (not enough usable data): %s\n", fac, ME.message);
                        continue
                    end

                    fixed_effects = lme.Coefficients;

                    Intercept = fixed_effects.Estimate(1);
                    Slope = fixed_effects.Estimate(2);
                    Se = fixed_effects.SE(2);
                    Tstat = fixed_effects.tStat(2);
                    p = fixed_effects.pValue(2);

                    % WM0=sum(Weights .* V0) / sum(Weights);
                    % WM1=sum(Weights .* V1) / sum(Weights);
                    % 
                    % WV0=sum(Weights .* (V0-WM0).^2)/sum(Weights);
                    % WV1=sum(Weights .* (V1-WM1).^2)/sum(Weights);
                    % 
                    % WCV = sum(Weights .* (V0 - WM0) .* (V1 - WM1)) / sum(Weights);
                    % 
                    % 
                    % %Linear Fit
                    % 
                    % XL = [ones(size(V0)), V0];  % Design matrix: intercept + slope
                    % 
                    % W = diag(Weights);           % Weight matrix
                    % 
                    % % Weighted least squares solution: (X'WX)^(-1) X'Wy
                    % betaL = (XL' * W * XL) \ (XL' * W * V1);
                    % 
                    % % Trendline values
                    % Age_line = linspace(min(V0), max(V0), 100)';
                    % VS_fit = betaL(1) + betaL(2) * Age_line;
                    % 
                    % RValL= WCV/sqrt(WV0*WV1);
                    % 
                    % ZL = 0.5 * log((1 + RValL) / (1 - RValL));
                    % 
                    % N_eff = (sum(Weights))^2 / sum(Weights.^2);
                    % 
                    % SE_z = 1 / sqrt(N_eff - 3);
                    % 
                    % ZStatL=ZL/SE_z;
                    % 
                    % pL=2*(1 - normcdf(abs(ZStatL)));
                    % 
                    % %Quadratic Fit
                    % 
                    % XQ = [V0.^2, V0, ones(size(V0))];
                    % 
                    % % Weighted least squares solution: (X'WX)^(-1) X'Wy
                    % betaQ = (XQ' * W * XQ) \ (XQ' * W * V1);
                    % 
                    % x_fit = linspace(min(V0), max(V0), 100);
                    % y_fit = betaQ(1)*x_fit.^2 + betaQ(2)*x_fit + betaQ(3);
                    % 
                    % y_pred = XQ * betaQ;
                    % residuals = V1 - y_pred;
                    % % Residual and total sum of squares (weighted)
                    % ybar = sum(Weights .* V1) / sum(Weights);  % weighted mean
                    % SS_res = sum(Weights .* residuals.^2);
                    % SS_tot = sum(Weights .* (V1 - ybar).^2);
                    % 
                    % % R-squared and adjusted R-squared
                    % R_squared = 1 - SS_res / SS_tot;
                    % adj_R_squared = 1 - (1 - R_squared) * (length(V1) - 1) / (length(V1) - size(XQ,2));
                    % 
                    % % Variance estimate (weighted MSE)
                    % sigma2 = SS_res / (length(V1) - size(XQ,2));
                    % 
                    % % Covariance matrix of beta estimates
                    % cov_beta = sigma2 * inv(XQ' * W * XQ);
                    % SE_beta = sqrt(diag(cov_beta));        % standard errors
                    % t_stats = betaQ ./ SE_beta;
                    % p_valsQ = 2 * (1 - tcdf(abs(t_stats), length(V1) - size(XQ,2)));  % two-sided p-values
                    % 
                    % % Overall model F-test
                    % F = (R_squared / (size(XQ,2) - 1)) / ((1 - R_squared) / (length(V1) - size(XQ,2)));
                    % pQ = 1 - fcdf(F, size(XQ,2) - 1, length(V1) - size(XQ,2));
                    % 
                    % [p,which]=min([pL,pQ]);

                    Age_line = linspace(min(data.Age), max(data.Age), 100)';
                    VS_fit = Intercept + Slope * Age_line;

                    fig=figure('Visible','off');
                    ax=gca; % captured once so the sig-box-only bonus save below can find every data child of THIS axes specifically
                    set(ax, 'box', 'off');
                    set(fig,'units','inches','outerposition',[0 0 11 6.4],'windowstyle','normal') % item 3, 2026-09-05: matches plotAgeTrendByGroup's height (room for the lifted title)
                    hold on
                    scatter(V0, V1, 40 * normalize(Weights, 'range')+plus, 'filled','black','MarkerFaceAlpha',0.9)
                    % scatter(V0, V1, 'filled','black','MarkerFaceAlpha',0.65)
                    plot(Age_line, VS_fit, 'LineWidth', 2,'Color','black');
                    % if (which==1)
                    %
                    % end
                    % if (which==2)
                    %     plot(x_fit, y_fit, 'LineWidth', 2,'Color','black');
                    % end
                    yline(0,'--','Color','k')

                    % eval(strcat("xticklabels({num2str(",Efields{i}(1:(strfind(Efields{i},"Thresh")-1)),"Mu0(pass)),num2str(",Efields{i}(1:(strfind(Efields{i},"Thresh")-1)),"Mu1(pass))});"))
                    if fac==1
                        callit=strcat(TitleName,"/","Rate_Age_Original");
                        logLabel="Rate Age Fit: p= ";
                        titleHandle=title(strcat(coe,"_\rho"));
                    else
                        callit=strcat(TitleName,"/",longnames2{fac-1},"_Age_Original");
                        logLabel=strcat(longnames2{fac-1}," Age Fit: p= ");
                        titleHandle=title(strcat(coe,formulanames(fac-1)));
                    end
                    titleHandle.Units='normalized';
                    titleHandle.Position(2)=titleHandle.Position(2)+0.02; % item 3, 2026-09-05: same lift as the other scatter/violin titles
                    % One companion text file per plot (same name, same
                    % folder as the PNG) instead of writing into the big
                    % shared per-scope log.
                    savePath=resolveSavePath(Settings.outputDir,"Individual Plots",callit);
                    txtFid=fopen(strcat(savePath,".txt"),'wt');
                    fprintf(txtFid,'%s',logLabel);
                    fprintf(txtFid,strcat("%5e","\n"),p);
                    fclose(txtFid);

                    fontsize(14,'points')

                    limx=xlim;
                    limy=ylim;
                    fracx=(limx(2)-limx(1))./20;
                    fracy=(limy(2)-limy(1))./20;
                    % Item 3, 2026-09-05: this box had its own hand-rolled
                    % star logic and plain (unbordered, unstyled) text,
                    % inconsistent with every other box in the file --
                    % brought in line with plotAgeTrendByGroup's box:
                    % boxSigLabel for the same soft-font stars/"n.s."
                    % styling, placeFitTextBox for the same corner
                    % placement, and the same real-extent-based axis
                    % clearance.
                    mylabel=strcat("Overall Fit ",boxSigLabel(p,12));
                    boxText=placeFitTextBox(limx(2)-0.2*fracx,limy(2)+3*fracy,mylabel,12,'right','top');
                    ext=boxText.Extent; % read before hiding below, in case Visible='off' ever affects Extent computation
                    boxText.Visible='off'; % 2026-09-05, user request: the box now lives only on the sig-box-only bonus image, turned back on right before that save
                    xlim([limx(1),max(limx(2)+4*fracx,ext(1)+ext(3)+0.5*fracx)]);
                    ylim([limy(1),max(limy(2)+4*fracy,ext(2)+ext(4)+0.5*fracy)]);
                    % Item 5 (prior round): same "nice ticks through zero"
                    % treatment as the other scatter/violin plots. Exponent
                    % forced off (item 2, prior round) -- MATLAB's own
                    % auto power-of-ten scaling applies OUR decimal format
                    % to the SCALED mantissa, not the raw tick value,
                    % producing artifacts like "5.0000" on a x10^4 scale.
                    % Follow-up (2026-09-05): ytickformat's own interaction
                    % with a forced Exponent turned out to be unreliable
                    % the OTHER way too (rounding small-scale values to
                    % "0.00" in some plots) -- explicit tick label STRINGS
                    % (computed ourselves, via yticklabels) bypass
                    % ytickformat/Exponent entirely, so nothing is left to
                    % second-guess our own decimal count.
                    [niceY,stepY]=niceTicksThroughZero(limy(1),limy(2),5);
                    decY=decimalsForStep(stepY);
                    yticks(niceY)
                    yticklabels(arrayfun(@(v) sprintf(strcat('%.',num2str(decY),'f'),v), niceY, 'UniformOutput', false))
                    xlabel("Age")
                    ylabel("Factor Weight")

                    % annotation('textbox',[.2 .74 .25 .15],'String',strcat("p-value=",num2str(round(p,4))),'FitBoxToText','on')

                    set(gca, 'box', 'off');
                    saveas(fig,strcat(savePath,".png"))

                    % "Bonus" clean sig-box-only version (2026-09-05, user
                    % request, same mechanism as plotGroupViolin/
                    % plotAgeTrendByGroup): same figure, same axes/
                    % scaling/ticks, but the scatter/fitted line/zero-line
                    % hidden, leaving only the "Overall Fit" box against
                    % the correctly-scaled, otherwise-empty axes. Border
                    % framing just the box text.
                    dataChildren=ax.Children(ax.Children~=boxText);
                    set(dataChildren,'Visible','off')
                    set(boxText,'EdgeColor',[0.3,0.3,0.3],'Margin',4,'Visible','on') % turned back on here -- hidden on the main save above, per user request
                    saveas(fig,strcat(savePath,"_SigBoxOnly.png"))

                end
            end
            close all
        end
        end

        %DIAGNOSIS WITH AGE (Phase B, 2026-09-05: sourced from T instead of
        % VSAS/AgeA/SexA/DisoA/IDA)
        for fac=1:size(T.Value,2)
            if ~all(T.Value(:,fac)==0)
                FactorWeight=T.Value(:,fac);
                Weights=inverseVarianceWeights(T.Error(:,fac));
                Age=T.Age;
                Sex=T.Sex;
                Impairment=T.Disorder;
                ID=T.ID;
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                GroupLabel=strings(height(tbl),1);
                GroupLabel(~Impairment)="CN";
                GroupLabel(Impairment)="DS";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Impairment + AgeCentered + Impairment:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Impairment + AgeCentered + AgeCentered^2 + Impairment:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.Disorder.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["CN","DS"];
                cfg.mainEffects=struct('pattern',{'^Impairment_[^:]+$','^AgeCentered(\^2)?$','^Impairment_[^:]+:AgeCentered$'}, ...
                    'label',{"Impairment","Age","Interaction"}, ...
                    'type',{"diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=8;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Diagnosis_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Diagnosis Age Fit (Mixed Model)";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Diagnosis_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Diagnosis Age Fit (Mixed Model)");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.Disorder.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"Disorder",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end

            end
        end
        close all

        %DIAGNOSIS (Phase B, 2026-09-05: sourced from T instead of
        % VSS0/VSS1/subsetGroup; test math moved into runWeightedTTest,
        % now mixed-effects like every other T-Test per user direction)
        maskControl=T.Disorder==0;
        maskDiagnosed=T.Disorder==1;
        V=T.Value(maskControl,:); Vv=T.Value(maskDiagnosed,:);
        Rv=T.Error(maskControl,:); Rr=T.Error(maskDiagnosed,:);
        I0=T.ID(maskControl); I1=T.ID(maskDiagnosed);

        DisorderGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0=V(:,fac); V1=Vv(:,fac); E0=Rv(:,fac); E1=Rr(:,fac);
            [p,tStat,dfStat]=runWeightedTTest(V0,E0,I0,V1,E1,I1);

            if ~isnan(p)
                cfg=struct();
                cfg.groupLevels=["CN","DS"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=2.5;
                comparisons=struct('leftIdx',1,'rightIdx',2,'pValue',p,'label',"Diagnosis Difference (Mixed Model)",'tStat',tStat,'df',dfStat);
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Diagnosis_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Diagnosis Difference (Mixed Model)";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Diagnosis_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Diagnosis Difference (Mixed Model)");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                DisorderGridData{fac}=struct('values',{{V0,V1}},'mu',[Mu0(fac),Mu1(fac)],'sigma',[Sigma0(fac),Sigma1(fac)], ...
                    'comparisons',struct('leftIdx',1,'rightIdx',2,'pValue',p));

                cfg.drawPlot=Settings.tests.Disorder.doViolinPlot;
                Effects=plotGroupViolin({V0,V1},[Mu0(fac),Mu1(fac)],[Sigma0(fac),Sigma1(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"Disorder",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
        end
        % Stashed for the Graphing section (below the main loop), which
        % draws this comparison's per-scope grid, the combined
        % Functional+Structural grid, and the subnetwork-aggregate grid --
        % all from this one Results entry, regardless of which toggles
        % are on.
        Results.Disorder.(scopeField)=DisorderGridData;

        %SEX WITH AGE (Phase B, 2026-09-05: sourced from T instead of
        % VSAS/AgeA/SexA/DisoA. Converted to mixed-effects 2026-09-05 per
        % user direction, joining Disorder's AgeTrend as ID-clustered.)
        for fac=1:size(T.Value,2)
            if ~all(T.Value(:,fac)==0)
                FactorWeight=T.Value(:,fac);
                Weights=inverseVarianceWeights(T.Error(:,fac));
                Age=T.Age;
                Sex=T.Sex;
                Impairment=T.Disorder;
                ID=T.ID;
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                GroupLabel=strings(height(tbl),1);
                GroupLabel(T.Sex==0)="M";
                GroupLabel(T.Sex==1)="F";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Sex+AgeCentered+Sex:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Sex+AgeCentered^2+Sex:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.Sex.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["M","F"];
                cfg.mainEffects=struct('pattern',{'^Sex_[^:]+$','^AgeCentered(\^2)?$','^Sex_[^:]+:AgeCentered$'}, ...
                    'label',{"Sex","Age","Interaction"}, ...
                    'type',{"diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=6;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Sex Age Fit";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Sex Age Fit");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.Sex.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"Sex",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end

            end
        end
        close all

        %SEX (Phase B, 2026-09-05: sourced from T instead of VSS2/VSS3/
        % subsetGroup; test math moved into runWeightedTTest, mixed-effects
        % like every other T-Test per user direction, and now ID-clustered
        % same as Disorder always was)
        maskMale=T.Sex==0;
        maskFemale=T.Sex==1;
        V=T.Value(maskMale,:); Vv=T.Value(maskFemale,:);
        Rv=T.Error(maskMale,:); Rr=T.Error(maskFemale,:);
        I0=T.ID(maskMale); I1=T.ID(maskFemale);

        SexGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0=V(:,fac); V1=Vv(:,fac); E0=Rv(:,fac); E1=Rr(:,fac);
            [p,tStat,dfStat]=runWeightedTTest(V0,E0,I0,V1,E1,I1);

            if ~isnan(p)
                cfg=struct();
                cfg.groupLevels=["M","F"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=2.5;
                comparisons=struct('leftIdx',1,'rightIdx',2,'pValue',p,'label',"Sex Difference (Mixed Model)",'tStat',tStat,'df',dfStat);
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Sex_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Sex Difference (Mixed Model)";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Sex_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Sex Difference (Mixed Model)");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                SexGridData{fac}=struct('values',{{V0,V1}},'mu',[Mu2(fac),Mu3(fac)],'sigma',[Sigma2(fac),Sigma3(fac)],'comparisons',comparisons);

                cfg.drawPlot=Settings.tests.Sex.doViolinPlot;
                Effects=plotGroupViolin({V0,V1},[Mu2(fac),Mu3(fac)],[Sigma2(fac),Sigma3(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"Sex",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
        end
        % Stashed for the Graphing section (below the main loop).
        Results.Sex.(scopeField)=SexGridData;



        %DIAGNOSIS AND SEX WITH AGE (Phase B, 2026-09-05: sourced from T
        % instead of VSAS/AgeA/SexA/DisoA. Converted to mixed-effects
        % 2026-09-05 per user direction.)
        for fac=1:size(T.Value,2)
            if ~all(T.Value(:,fac)==0)
                FactorWeight=T.Value(:,fac);
                Weights=inverseVarianceWeights(T.Error(:,fac));
                Age=T.Age;
                Sex=T.Sex;
                Impairment=T.Disorder;
                ID=T.ID;
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                GroupLabel=strings(height(tbl),1);
                GroupLabel((T.Sex==0)&(~T.Disorder))="MCN";
                GroupLabel((T.Sex==0)&(T.Disorder==1))="MDS";
                GroupLabel((T.Sex==1)&(~T.Disorder))="FCN";
                GroupLabel((T.Sex==1)&(T.Disorder==1))="FDS";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Impairment+Sex+AgeCentered+Sex:AgeCentered+Impairment:AgeCentered+Impairment:Sex:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Impairment+Sex+AgeCentered^2+Sex:AgeCentered+Impairment:AgeCentered+Impairment:Sex:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.SexxDisorder.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["MCN","MDS","FCN","FDS"];
                cfg.mainEffects=struct('pattern',{'^Impairment_[^:]+$','^Sex_[^:]+$','^AgeCentered(\^2)?$','^Sex_[^:]+:AgeCentered$'}, ...
                    'label',{"Impairment","Sex","Age","Interaction"}, ...
                    'type',{"diff","diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=6;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Diagnosis_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Diagnosis Sex Age Fit";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Diagnosis_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Diagnosis Sex Age Fit");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.SexxDisorder.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"SexxDisorder",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end

            end
        end
        close all

        %DIAGNOSIS AND SEX (Phase B, 2026-09-05: sourced from T instead of
        % VSS4-7/subsetGroup; F-test + pairwise math moved into
        % runWeightedTwoWayAnova. Converted to mixed-effects 2026-09-05
        % per user direction (fitlme+(1|ID), pooled-variance contrasts,
        % Bonferroni correction).)
        maskMC=(T.Sex==0)&(T.Disorder==0); maskMD=(T.Sex==0)&(T.Disorder==1);
        maskFC=(T.Sex==1)&(T.Disorder==0); maskFD=(T.Sex==1)&(T.Disorder==1);
        V=T.Value(maskMC,:); Vv=T.Value(maskMD,:); Vvv=T.Value(maskFC,:); Vvvv=T.Value(maskFD,:);
        R=T.Error(maskMC,:); Rr=T.Error(maskMD,:); Rrr=T.Error(maskFC,:); Rrrr=T.Error(maskFD,:);

        maskAllSxD=maskMC|maskMD|maskFC|maskFD;
        Factor1AllSxD=T.Sex(maskAllSxD);
        Factor2AllSxD=T.Disorder(maskAllSxD);
        IDAllSxD=T.ID(maskAllSxD);

        SexxDisorderGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0 = V(:,fac);
            V1 = Vv(:,fac);
            V2 = Vvv(:,fac);
            V3 = Vvvv(:,fac);
            E0 = R(:,fac);
            E1 = Rr(:,fac);
            E2 = Rrr(:,fac);
            E3 = Rrrr(:,fac);

            W0=1./(E0.^2);
            W1=1./(E1.^2);
            W2=1./(E2.^2);
            W3=1./(E3.^2);

            FactorWeightAllSxD=T.Value(maskAllSxD,fac);
            % Skip metrics that don't apply to this scope at all (value=0
            % for every subject, a structural placeholder rather than a
            % real measurement) instead of handing fitlme an all-NaN-weight
            % table -- same guard the T-Test/AgeTrend blocks already use
            % (2026-09-05).
            if ~all(FactorWeightAllSxD==0)
            WeightsAllSxD=inverseVarianceWeights(T.Error(maskAllSxD,fac));
            aov=runWeightedTwoWayAnova(FactorWeightAllSxD,WeightsAllSxD,Factor1AllSxD,Factor2AllSxD,IDAllSxD,["M","F"],["CN","DS"]);

            if ~isnan(aov.pI)
                cfg=struct();
                cfg.groupLevels=["MCN","MDS","FCN","FDS"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=4.7; % item 7, 2026-09-05: extra room for the box
                comparisons=struct( ...
                    'leftIdx',{1,3,1,2}, ...
                    'rightIdx',{2,4,3,4}, ...
                    'pValue',{aov.pA,aov.pB,aov.pC,aov.pD}, ...
                    'label',{"Male Impairment","Female Impairment","Control Sex","Impaired Sex"}, ...
                    'tStat',{aov.tStatA,aov.tStatB,aov.tStatC,aov.tStatD}, ...
                    'df',{aov.dfA,aov.dfB,aov.dfC,aov.dfD});
                cfg.overallEffects=struct('pValue',{aov.pCD,aov.pAB,aov.pI},'label',{"Impairment","Sex","Interaction"}, ...
                    'tStat',{sqrt(aov.FStatCD),sqrt(aov.FStatAB),sqrt(aov.FStatI)},'df',{aov.DFR,aov.DFR,aov.DFR});
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Diagnosis_Sex_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Diagnosis and Sex Differences";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Diagnosis_Sex_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Diagnosis and Sex Differences");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                SexxDisorderGridData{fac}=struct('values',{{V0,V1,V2,V3}},'mu',[Mu4(fac),Mu5(fac),Mu6(fac),Mu7(fac)], ...
                    'sigma',[Sigma4(fac),Sigma5(fac),Sigma6(fac),Sigma7(fac)],'comparisons',comparisons,'overallEffects',cfg.overallEffects); % 2026-09-05, user request: carries the sig box through to the grid/aggregate plots

                cfg.drawPlot=Settings.tests.SexxDisorder.doViolinPlot;
                Effects=plotGroupViolin({V0,V1,V2,V3},[Mu4(fac),Mu5(fac),Mu6(fac),Mu7(fac)],[Sigma4(fac),Sigma5(fac),Sigma6(fac),Sigma7(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"SexxDisorder",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
            end
        end
        % Stashed for the Graphing section (below the main loop).
        Results.SexxDisorder.(scopeField)=SexxDisorderGridData;

        %DIAGNOSIS WITH AGE (MCI; Phase B, 2026-09-05: sourced from T
        % instead of VSAS/AgeA/SexA/DisoA/DiagA. Converted to mixed-effects
        % 2026-09-05 per user direction.)
        maskCNorMCI=T.Diagnosis<2;
        for fac=1:size(T.Value,2)
            if ~all(T.Value(maskCNorMCI,fac)==0)
                FactorWeight=T.Value(maskCNorMCI,fac);
                Weights=inverseVarianceWeights(T.Error(maskCNorMCI,fac));
                Age=T.Age(maskCNorMCI);
                Sex=T.Sex(maskCNorMCI);
                Impairment=T.Disorder(maskCNorMCI);
                ID=T.ID(maskCNorMCI);
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                GroupLabel=strings(height(tbl),1);
                GroupLabel(~Impairment)="CN";
                GroupLabel(Impairment)="MCI";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Impairment+AgeCentered+Impairment:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Impairment+AgeCentered^2+Impairment:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.MCI.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["CN","MCI"];
                cfg.mainEffects=struct('pattern',{'^Impairment_[^:]+$','^AgeCentered(\^2)?$','^Impairment_[^:]+:AgeCentered$'}, ...
                    'label',{"Mild Impairment","Age","Interaction"}, ...
                    'type',{"diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=8;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_MCI_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate MCI Age Fit";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_MCI_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," MCI Age Fit");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.MCI.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"MCI",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end

            end
        end
        close all

        %DIAGNOSIS (MCI; Phase B, 2026-09-05: sourced from T instead of
        % VSS8/VSS9/subsetGroup; test math moved into runWeightedTTest,
        % mixed-effects like every other T-Test per user direction)
        maskControlMCI=T.Diagnosis==0;
        maskMCI=T.Diagnosis==1;
        V=T.Value(maskControlMCI,:); Vv=T.Value(maskMCI,:);
        Rv=T.Error(maskControlMCI,:); Rr=T.Error(maskMCI,:);
        I0=T.ID(maskControlMCI); I1=T.ID(maskMCI);

        MCIGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0=V(:,fac); V1=Vv(:,fac); E0=Rv(:,fac); E1=Rr(:,fac);
            [p,tStat,dfStat]=runWeightedTTest(V0,E0,I0,V1,E1,I1);

            if ~isnan(p)
                cfg=struct();
                cfg.groupLevels=["CN","MCI"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=2.5;
                comparisons=struct('leftIdx',1,'rightIdx',2,'pValue',p,'label',"MCI Difference (Mixed Model)",'tStat',tStat,'df',dfStat);
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_MCI_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate MCI Difference (Mixed Model)";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_MCI_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," MCI Difference (Mixed Model)");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                MCIGridData{fac}=struct('values',{{V0,V1}},'mu',[Mu8(fac),Mu9(fac)],'sigma',[Sigma8(fac),Sigma9(fac)],'comparisons',comparisons);

                cfg.drawPlot=Settings.tests.MCI.doViolinPlot;
                Effects=plotGroupViolin({V0,V1},[Mu8(fac),Mu9(fac)],[Sigma8(fac),Sigma9(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"MCI",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
        end
        % Stashed for the Graphing section (below the main loop).
        Results.MCI.(scopeField)=MCIGridData;

        %DIAGNOSIS WITH AGE (AD; Phase B, 2026-09-05: sourced from T
        % instead of VSAS/AgeA/SexA/DisoA/DiagA. Converted to mixed-effects
        % 2026-09-05 per user direction.)
        maskCNorAD=T.Diagnosis~=1;
        for fac=1:size(T.Value,2)
            if ~all(T.Value(maskCNorAD,fac)==0)
                FactorWeight=T.Value(maskCNorAD,fac);
                Weights=inverseVarianceWeights(T.Error(maskCNorAD,fac));
                Age=T.Age(maskCNorAD);
                Sex=T.Sex(maskCNorAD);
                Impairment=T.Disorder(maskCNorAD);
                ID=T.ID(maskCNorAD);
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                GroupLabel=strings(height(tbl),1);
                GroupLabel(~Impairment)="CN";
                GroupLabel(Impairment)="AD";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Impairment+AgeCentered+Impairment:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Impairment+AgeCentered^2+Impairment:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.AD.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["CN","AD"];
                cfg.mainEffects=struct('pattern',{'^Impairment_[^:]+$','^AgeCentered(\^2)?$','^Impairment_[^:]+:AgeCentered$'}, ...
                    'label',{"Significant Impairment","Age","Interaction"}, ...
                    'type',{"diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=8;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_AD_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate AD Age Fit";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_AD_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," AD Age Fit");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.AD.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"AD",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end

            end
        end
        close all

        %DIAGNOSIS (AD; Phase B, 2026-09-05: sourced from T instead of
        % VSS8/VSS10/subsetGroup; test math moved into runWeightedTTest,
        % mixed-effects like every other T-Test per user direction)
        maskControlAD=T.Diagnosis==0;
        maskAD=T.Diagnosis==2;
        V=T.Value(maskControlAD,:); Vv=T.Value(maskAD,:);
        Rv=T.Error(maskControlAD,:); Rr=T.Error(maskAD,:);
        I0=T.ID(maskControlAD); I1=T.ID(maskAD);

        ADGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0=V(:,fac); V1=Vv(:,fac); E0=Rv(:,fac); E1=Rr(:,fac);
            [p,tStat,dfStat]=runWeightedTTest(V0,E0,I0,V1,E1,I1);

            if ~isnan(p)
                cfg=struct();
                cfg.groupLevels=["CN","AD"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=2.5;
                comparisons=struct('leftIdx',1,'rightIdx',2,'pValue',p,'label',"AD Difference (Mixed Model)",'tStat',tStat,'df',dfStat);
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_AD_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate AD Difference (Mixed Model)";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_AD_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," AD Difference (Mixed Model)");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                ADGridData{fac}=struct('values',{{V0,V1}},'mu',[Mu8(fac),Mu10(fac)],'sigma',[Sigma8(fac),Sigma10(fac)],'comparisons',comparisons);

                cfg.drawPlot=Settings.tests.AD.doViolinPlot;
                Effects=plotGroupViolin({V0,V1},[Mu8(fac),Mu10(fac)],[Sigma8(fac),Sigma10(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"AD",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
        end
        % Stashed for the Graphing section (below the main loop).
        Results.AD.(scopeField)=ADGridData;

        %SEVERITY: CONTROL VS MCI VS AD (weighted one-way ANOVA + pairwise;
        % Phase B, 2026-09-05: sourced from T instead of VSS8/VSS9/VSS10/
        % subsetGroup; F-test moved into runWeightedOneWayAnova. Converted
        % to mixed-effects 2026-09-05 per user direction (fitlme+(1|ID),
        % pooled-variance contrasts, Bonferroni correction) -- retires the
        % separate weightedPairwiseTTest calls below in favor of the
        % contrasts returned alongside the omnibus test.)
        maskSevControl=T.Diagnosis==0; maskSevMCI=T.Diagnosis==1; maskSevAD=T.Diagnosis==2;
        VSev0=T.Value(maskSevControl,:); ESev0=T.Error(maskSevControl,:);
        VSev1=T.Value(maskSevMCI,:); ESev1=T.Error(maskSevMCI,:);
        VSev2=T.Value(maskSevAD,:); ESev2=T.Error(maskSevAD,:);

        maskAllSev=maskSevControl|maskSevMCI|maskSevAD;
        GroupAllSev=T.Diagnosis(maskAllSev);
        IDAllSev=T.ID(maskAllSev);

        SeverityGridData=cell(1,size(VSev0,2));
        for fac=1:size(VSev0,2)
            V0=VSev0(:,fac); E0=ESev0(:,fac);
            V1=VSev1(:,fac); E1=ESev1(:,fac);
            V2=VSev2(:,fac); E2=ESev2(:,fac);

            W0=1./(E0.^2);
            W1=1./(E1.^2);
            W2=1./(E2.^2);

            FactorWeightAllSev=T.Value(maskAllSev,fac);
            % Skip metrics that don't apply to this scope at all (value=0
            % for every subject, a structural placeholder rather than a
            % real measurement) instead of handing fitlme an all-NaN-weight
            % table -- same guard the T-Test/AgeTrend blocks already use
            % (2026-09-05).
            if ~all(FactorWeightAllSev==0)
            WeightsAllSev=inverseVarianceWeights(T.Error(maskAllSev,fac));
            aov=runWeightedOneWayAnova(FactorWeightAllSev,WeightsAllSev,GroupAllSev,IDAllSev,["CN","MCI","AD"]);
            pANOVA=aov.pOmnibus; DFW=aov.DF2Omnibus;
            pCM=aov.pCM; tCM=aov.tStatCM; dfCM=aov.dfCM;
            pCA=aov.pCA; tCA=aov.tStatCA; dfCA=aov.dfCA;
            pMA=aov.pMA; tMA=aov.tStatMA; dfMA=aov.dfMA;

            if ~isnan(pANOVA)
                comparisons=struct( ...
                    'leftIdx',{1,2,1}, ...
                    'rightIdx',{2,3,3}, ...
                    'pValue',{pCM,pMA,pCA}, ...
                    'label',{"Control vs MCI","MCI vs AD","Control vs AD"}, ...
                    'tStat',{tCM,tMA,tCA}, ...
                    'df',{dfCM,dfMA,dfCA});
                cfg=struct();
                cfg.groupLevels=["CN","MCI","AD"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=3.7; % item 7, 2026-09-05: extra room for the box
                cfg.overallEffects=struct('pValue',pANOVA,'label',"Severity",'tStat',NaN,'df',DFW); % item 5, 2026-09-05: dropped "(Mixed Model ANOVA)"
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_Severity_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate Severity (Control/MCI/AD) ANOVA";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_Severity_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," Severity (Control/MCI/AD) ANOVA");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                SeverityGridData{fac}=struct('values',{{V0,V1,V2}},'mu',[Mu8(fac),Mu9(fac),Mu10(fac)], ...
                    'sigma',[Sigma8(fac),Sigma9(fac),Sigma10(fac)],'comparisons',comparisons,'overallEffects',cfg.overallEffects); % 2026-09-05, user request: carries the sig box through to the grid/aggregate plots

                cfg.drawPlot=Settings.tests.Severity.doViolinPlot;
                Effects=plotGroupViolin({V0,V1,V2},[Mu8(fac),Mu9(fac),Mu10(fac)],[Sigma8(fac),Sigma9(fac),Sigma10(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"Severity",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
            end
        end
        close all

        % Stashed for the Graphing section (below the main loop).
        Results.Severity.(scopeField)=SeverityGridData;

        %DIAGNOSIS AND SEX WITH AGE (MCI; Phase B, 2026-09-05: sourced
        % from T instead of VSAS/AgeA/SexA/DisoA/DiagA. Converted to
        % mixed-effects 2026-09-05 per user direction.)
        maskCNorMCI2=T.Diagnosis<2;
        for fac=1:size(T.Value,2)
            if ~all(T.Value(maskCNorMCI2,fac)==0)
                FactorWeight=T.Value(maskCNorMCI2,fac);
                Weights=inverseVarianceWeights(T.Error(maskCNorMCI2,fac));
                Age=T.Age(maskCNorMCI2);
                Sex=T.Sex(maskCNorMCI2);
                Impairment=T.Disorder(maskCNorMCI2);
                ID=T.ID(maskCNorMCI2);
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                SexSub=T.Sex(maskCNorMCI2);
                DisoSub=T.Disorder(maskCNorMCI2);
                GroupLabel=strings(height(tbl),1);
                GroupLabel((SexSub==0)&(~DisoSub))="MCN";
                GroupLabel((SexSub==0)&(DisoSub==1))="MMCI";
                GroupLabel((SexSub==1)&(~DisoSub))="FCN";
                GroupLabel((SexSub==1)&(DisoSub==1))="FMCI";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Impairment+Sex+AgeCentered+Sex:AgeCentered+Impairment:AgeCentered+Impairment:Sex:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Impairment+Sex+AgeCentered^2+Sex:AgeCentered+Impairment:AgeCentered+Impairment:Sex:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.SexxMCI.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["MCN","MMCI","FCN","FMCI"];
                cfg.mainEffects=struct('pattern',{'^Impairment_[^:]+$','^Sex_[^:]+$','^AgeCentered(\^2)?$','^Sex_[^:]+:AgeCentered$'}, ...
                    'label',{"Mild Impairment","Sex","Age","Interaction"}, ...
                    'type',{"diff","diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=10;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_MCI_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate MCI Sex Age Fit";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_MCI_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," MCI Sex Age Fit");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.SexxMCI.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"SexxMCI",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end

            end
        end
        close all

        %DIAGNOSIS AND SEX (MCI; Phase B, 2026-09-05: sourced from T
        % instead of VSS11-14/subsetGroup; F-test + pairwise math moved
        % into runWeightedTwoWayAnova. Converted to mixed-effects
        % 2026-09-05 per user direction (fitlme+(1|ID), pooled-variance
        % contrasts, Bonferroni correction).)
        maskMCn=(T.Sex==0)&(T.Diagnosis==0); maskMMCI=(T.Sex==0)&(T.Diagnosis==1);
        maskFCn=(T.Sex==1)&(T.Diagnosis==0); maskFMCI=(T.Sex==1)&(T.Diagnosis==1);
        V=T.Value(maskMCn,:); Vv=T.Value(maskMMCI,:); Vvv=T.Value(maskFCn,:); Vvvv=T.Value(maskFMCI,:);
        R=T.Error(maskMCn,:); Rr=T.Error(maskMMCI,:); Rrr=T.Error(maskFCn,:); Rrrr=T.Error(maskFMCI,:);

        maskAllSxDg=maskMCn|maskMMCI|maskFCn|maskFMCI;
        Factor1AllSxDg=T.Sex(maskAllSxDg);
        Factor2AllSxDg=T.Diagnosis(maskAllSxDg);
        IDAllSxDg=T.ID(maskAllSxDg);

        SexxMCIGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0 = V(:,fac);
            V1 = Vv(:,fac);
            V2 = Vvv(:,fac);
            V3 = Vvvv(:,fac);
            E0 = R(:,fac);
            E1 = Rr(:,fac);
            E2 = Rrr(:,fac);
            E3 = Rrrr(:,fac);

            W0=1./(E0.^2);
            W1=1./(E1.^2);
            W2=1./(E2.^2);
            W3=1./(E3.^2);

            FactorWeightAllSxDg=T.Value(maskAllSxDg,fac);
            % Skip metrics that don't apply to this scope at all (value=0
            % for every subject, a structural placeholder rather than a
            % real measurement) instead of handing fitlme an all-NaN-weight
            % table -- same guard the T-Test/AgeTrend blocks already use
            % (2026-09-05).
            if ~all(FactorWeightAllSxDg==0)
            WeightsAllSxDg=inverseVarianceWeights(T.Error(maskAllSxDg,fac));
            aov=runWeightedTwoWayAnova(FactorWeightAllSxDg,WeightsAllSxDg,Factor1AllSxDg,Factor2AllSxDg,IDAllSxDg,["M","F"],["CN","MCI"]);

            if ~isnan(aov.pI)
                cfg=struct();
                cfg.groupLevels=["MCN","MMCI","FCN","FMCI"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=4.7; % item 7, 2026-09-05: extra room for the box
                comparisons=struct( ...
                    'leftIdx',{1,3,1,2}, ...
                    'rightIdx',{2,4,3,4}, ...
                    'pValue',{aov.pA,aov.pB,aov.pC,aov.pD}, ...
                    'label',{"Male Mild Impairment","Female Mild Impairment","Control Sex","Mildly Impaired Sex"}, ...
                    'tStat',{aov.tStatA,aov.tStatB,aov.tStatC,aov.tStatD}, ...
                    'df',{aov.dfA,aov.dfB,aov.dfC,aov.dfD});
                cfg.overallEffects=struct('pValue',{aov.pCD,aov.pAB,aov.pI},'label',{"Mild Impairment","Sex","Interaction"}, ...
                    'tStat',{sqrt(aov.FStatCD),sqrt(aov.FStatAB),sqrt(aov.FStatI)},'df',{aov.DFR,aov.DFR,aov.DFR});
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_MCI_Sex_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate MCI and Sex Differences";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_MCI_Sex_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," MCI and Sex Differences");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                SexxMCIGridData{fac}=struct('values',{{V0,V1,V2,V3}},'mu',[Mu11(fac),Mu12(fac),Mu13(fac),Mu14(fac)], ...
                    'sigma',[Sigma11(fac),Sigma12(fac),Sigma13(fac),Sigma14(fac)],'comparisons',comparisons,'overallEffects',cfg.overallEffects); % 2026-09-05, user request: carries the sig box through to the grid/aggregate plots

                cfg.drawPlot=Settings.tests.SexxMCI.doViolinPlot;
                Effects=plotGroupViolin({V0,V1,V2,V3},[Mu11(fac),Mu12(fac),Mu13(fac),Mu14(fac)],[Sigma11(fac),Sigma12(fac),Sigma13(fac),Sigma14(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"SexxMCI",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
            end
        end
        % Stashed for the Graphing section (below the main loop).
        Results.SexxMCI.(scopeField)=SexxMCIGridData;

        %DIAGNOSIS AND SEX WITH AGE (AD; Phase B, 2026-09-05: sourced
        % from T instead of VSAS/AgeA/SexA/DisoA/DiagA. Converted to
        % mixed-effects 2026-09-05 per user direction.)
        maskCNorAD2=T.Diagnosis~=1;
        for fac=1:size(T.Value,2)
            if ~all(T.Value(maskCNorAD2,fac)==0)
                FactorWeight=T.Value(maskCNorAD2,fac);
                Weights=inverseVarianceWeights(T.Error(maskCNorAD2,fac));
                Age=T.Age(maskCNorAD2);
                Sex=T.Sex(maskCNorAD2);
                Impairment=T.Disorder(maskCNorAD2);
                ID=T.ID(maskCNorAD2);
                tbl = table(FactorWeight, Impairment, Sex, Age, Weights, ID);
                tbl.Impairment = categorical(tbl.Impairment);
                tbl.Sex = categorical(tbl.Sex);
                tbl.ID = categorical(tbl.ID);
                tbl.AgeCentered = tbl.Age - mean(tbl.Age);
                SexSub=T.Sex(maskCNorAD2);
                DisoSub=T.Disorder(maskCNorAD2);
                GroupLabel=strings(height(tbl),1);
                GroupLabel((SexSub==0)&(~DisoSub))="MCN";
                GroupLabel((SexSub==0)&(DisoSub==1))="MAD";
                GroupLabel((SexSub==1)&(~DisoSub))="FCN";
                GroupLabel((SexSub==1)&(DisoSub==1))="FAD";
                tbl.GroupLabel=categorical(GroupLabel);

                cfg=struct();
                cfg.formulaLinear="FactorWeight ~ Impairment+Sex+AgeCentered+Sex:AgeCentered+Impairment:AgeCentered+Impairment:Sex:AgeCentered + (1|ID)";
                cfg.formulaQuadratic="FactorWeight ~ Impairment+Sex+AgeCentered^2+Sex:AgeCentered+Impairment:AgeCentered+Impairment:Sex:AgeCentered + (1|ID)";
                cfg.includeQuadratic=Settings.tests.SexxAD.includeQuadratic;
                cfg.useMixedEffects=true;
                cfg.groupLevels=["MCN","MAD","FCN","FAD"];
                cfg.mainEffects=struct('pattern',{'^Impairment_[^:]+$','^Sex_[^:]+$','^AgeCentered(\^2)?$','^Sex_[^:]+:AgeCentered$'}, ...
                    'label',{"Significant Impairment","Sex","Age","Interaction"}, ...
                    'type',{"diff","diff","slope","diff"});
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.democheat=democheat;
                cfg.markerSizeOffset=10;
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_AD_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate AD Sex Age Fit";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_AD_Sex_Age_Fit");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," AD Sex Age Fit");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                cfg.drawPlot=Settings.tests.SexxAD.doAgeTrendPlot;
                Effects=plotAgeTrendByGroup(tbl,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"SexxAD",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
        end
        close all

        %DIAGNOSIS AND SEX (AD; Phase B, 2026-09-05: sourced from T
        % instead of VSS11/13/15/16/subsetGroup; F-test + pairwise math
        % moved into runWeightedTwoWayAnova. Converted to mixed-effects
        % 2026-09-05 per user direction (fitlme+(1|ID), pooled-variance
        % contrasts, Bonferroni correction). Diagnosis is 0/2 (CN/AD) in
        % this comparison, so it's remapped to 0/1 for the Factor2 input.)
        maskMCad=(T.Sex==0)&(T.Diagnosis==0); maskMAD=(T.Sex==0)&(T.Diagnosis==2);
        maskFCad=(T.Sex==1)&(T.Diagnosis==0); maskFAD=(T.Sex==1)&(T.Diagnosis==2);
        V=T.Value(maskMCad,:); Vv=T.Value(maskMAD,:); Vvv=T.Value(maskFCad,:); Vvvv=T.Value(maskFAD,:);
        R=T.Error(maskMCad,:); Rr=T.Error(maskMAD,:); Rrr=T.Error(maskFCad,:); Rrrr=T.Error(maskFAD,:);

        maskAllSxAD=maskMCad|maskMAD|maskFCad|maskFAD;
        Factor1AllSxAD=T.Sex(maskAllSxAD);
        Factor2AllSxAD=double(T.Diagnosis(maskAllSxAD)==2);
        IDAllSxAD=T.ID(maskAllSxAD);

        SexxADGridData=cell(1,size(V,2));
        for fac=1:size(V,2)
            V0 = V(:,fac);
            V1 = Vv(:,fac);
            V2 = Vvv(:,fac);
            V3 = Vvvv(:,fac);
            E0 = R(:,fac);
            E1 = Rr(:,fac);
            E2 = Rrr(:,fac);
            E3 = Rrrr(:,fac);

            W0=1./(E0.^2);
            W1=1./(E1.^2);
            W2=1./(E2.^2);
            W3=1./(E3.^2);

            FactorWeightAllSxAD=T.Value(maskAllSxAD,fac);
            % Skip metrics that don't apply to this scope at all (value=0
            % for every subject, a structural placeholder rather than a
            % real measurement) instead of handing fitlme an all-NaN-weight
            % table -- same guard the T-Test/AgeTrend blocks already use
            % (2026-09-05).
            if ~all(FactorWeightAllSxAD==0)
            WeightsAllSxAD=inverseVarianceWeights(T.Error(maskAllSxAD,fac));
            aov=runWeightedTwoWayAnova(FactorWeightAllSxAD,WeightsAllSxAD,Factor1AllSxAD,Factor2AllSxAD,IDAllSxAD,["M","F"],["CN","AD"]);

            if ~isnan(aov.pI)
                cfg=struct();
                cfg.groupLevels=["MCN","MAD","FCN","FAD"];
                cfg.demohex=demohex;
                cfg.democon=democon;
                cfg.demoorder=demoorder;
                cfg.democheat=democheat;
                cfg.clipFirstAtZero=(fac==1);
                cfg.figureWidthInches=4.7; % item 7, 2026-09-05: extra room for the box
                comparisons=struct( ...
                    'leftIdx',{1,3,1,2}, ...
                    'rightIdx',{2,4,3,4}, ...
                    'pValue',{aov.pA,aov.pB,aov.pC,aov.pD}, ...
                    'label',{"Male Significant Impairment","Female Significant Impairment","Control Sex","Significantly Impaired Sex"}, ...
                    'tStat',{aov.tStatA,aov.tStatB,aov.tStatC,aov.tStatD}, ...
                    'df',{aov.dfA,aov.dfB,aov.dfC,aov.dfD});
                cfg.overallEffects=struct('pValue',{aov.pCD,aov.pAB,aov.pI},'label',{"Significant Impairment","Sex","Interaction"}, ...
                    'tStat',{sqrt(aov.FStatCD),sqrt(aov.FStatAB),sqrt(aov.FStatI)},'df',{aov.DFR,aov.DFR,aov.DFR});
                if fac==1
                    callit=strcat(TitleName1,"/","Rate_AD_Sex_Box");
                    cfg.titleStr=strcat(coe,"_\rho");
                    cfg.logHeader="Rate AD and Sex Differences";
                else
                    callit=strcat(TitleName1,"/",longnames2{fac-1},"_AD_Sex_Box");
                    cfg.titleStr=strcat(coe,formulanames(fac-1));
                    cfg.logHeader=strcat(longnames2{fac-1}," AD and Sex Differences");
                end
                cfg.saveFile=resolveSavePath(Settings.outputDir,"Individual Plots",callit);

                SexxADGridData{fac}=struct('values',{{V0,V1,V2,V3}},'mu',[Mu11(fac),Mu15(fac),Mu13(fac),Mu16(fac)], ...
                    'sigma',[Sigma11(fac),Sigma15(fac),Sigma13(fac),Sigma16(fac)],'comparisons',comparisons,'overallEffects',cfg.overallEffects); % 2026-09-05, user request: carries the sig box through to the grid/aggregate plots

                cfg.drawPlot=Settings.tests.SexxAD.doViolinPlot;
                Effects=plotGroupViolin({V0,V1,V2,V3},[Mu11(fac),Mu15(fac),Mu13(fac),Mu16(fac)],[Sigma11(fac),Sigma15(fac),Sigma13(fac),Sigma16(fac)],comparisons,cfg);
                for e=1:height(Effects)
                    EffectsResults=appendEffectRow(EffectsResults,TitleName1,"SexxAD",Effects.Term(e),Effects.GroupName(e),MetricNames(fac),Effects.TStat(e),Effects.DF(e),Effects.PValue(e),Effects.Type(e));
                end
            end
            end
        end
        % Stashed for the Graphing section (below the main loop).
        Results.SexxAD.(scopeField)=SexxADGridData;
    end % closes "if notskip"
    fclose(fid);
end % closes "for i=indeces'"

% Sample-size diagnostic (Phase A addition, 2026-09-05): row count before
% and after NaN/Inf/outlier/zero-error removal, per scope, in Datasets.
scopeFields=fieldnames(Datasets);
SampleSizeSummary=table('Size',[numel(scopeFields),5], ...
    'VariableTypes',{'string','string','double','double','double'}, ...
    'VariableNames',{'Scope','SetName','NBefore','NAfter','PercentRemoved'});
for s=1:numel(scopeFields)
    d=Datasets.(scopeFields{s});
    SampleSizeSummary.Scope(s)=d.Label;
    SampleSizeSummary.SetName(s)=d.SetName;
    SampleSizeSummary.NBefore(s)=d.nBefore;
    SampleSizeSummary.NAfter(s)=d.nAfter;
    SampleSizeSummary.PercentRemoved(s)=100*(d.nBefore-d.nAfter)/d.nBefore;
end
writetable(SampleSizeSummary, strcat(Settings.outputDir,"SampleSizeSummary.csv"));

%% ===== GRAPHING (Phase C, 2026-09-05) =====
% Everything below reads from Results/Datasets, built up during the loop
% above -- no plot function itself changed, just when/from-what each gets
% called. Individual per-metric plots (violin, age-trend scatter) still
% get called inline during the loop above, since they're tightly bound to
% whichever single scope is being processed and don't benefit from
% deferring; only the grid, combined, and subnetwork-aggregate plots (which
% need every scope's data at once) are gathered here.

% Per-scope grid plots, one comparison at a time, over every scope Results
% has an entry for (whole-brain and subnetwork alike).
comparisonGridSpecs=struct( ...
    'name',{"Disorder","Sex","SexxDisorder","MCI","AD","Severity","SexxMCI","SexxAD"}, ...
    'groupLevels',{["CN","DS"],["M","F"],["MCN","MDS","FCN","FDS"],["CN","MCI"],["CN","AD"],["CN","MCI","AD"],["MCN","MMCI","FCN","FMCI"],["MCN","MAD","FCN","FAD"]}, ...
    'saveName',{"ViolinGridTest_Disorder","ViolinGridTest_Sex","ViolinGridTest_SexxDisorder","ViolinGridTest_MCI","ViolinGridTest_AD","ViolinGridTest_Severity","ViolinGridTest_SexxMCI","ViolinGridTest_SexxAD"}, ...
    'titlePrefix',{"Disorder: Control vs Impaired","Sex: Male vs Female","Sex x Disorder","MCI: Control vs Mild Impairment","AD: Control vs Significant Impairment","Severity: Control vs MCI vs AD","Sex x MCI","Sex x AD"});
for c=1:numel(comparisonGridSpecs)
    spec=comparisonGridSpecs(c);
    if ~Settings.tests.(spec.name).doGridPlot
        continue
    end
    scopeFields=fieldnames(Results.(spec.name));
    for s=1:numel(scopeFields)
        sf=scopeFields{s};
        label=Datasets.(sf).Label;
        plotFactorViolinGrid(Results.(spec.name).(sf),MetricNames,spec.groupLevels,demohex,democon,demoorder,democheat, ...
            resolveSavePath(Settings.outputDir,"Grid Plots and Heat Maps",strcat(label,"/",spec.saveName)),strcat(spec.titlePrefix," (",label,")"));
    end
end

% Combined Functional+Structural violin grid per comparison, built once
% both whole-brain scopes have been processed this run -- skipped (not
% errored) if either one wasn't, e.g. Settings.runFunctionalWholeBrain was off.
for c=1:numel(comparisonGridSpecs)
    spec=comparisonGridSpecs(c);
    % 2026-09-05, toggle audit: this loop (and the subnetwork-aggregate one
    % below) used to ignore doGridPlot entirely, so a comparison turned off
    % there (SexxMCI/SexxAD) still got its Combined and
    % All-Subnetworks grids generated -- only the per-scope grid actually
    % respected the toggle. Guarding all three the same way makes
    % doGridPlot control everything it visually looks like it should.
    if ~Settings.tests.(spec.name).doGridPlot
        continue
    end
    if isfield(Results.(spec.name),'Functional') && isfield(Results.(spec.name),'Structural')
        plotCombinedViolinGrid(Results.(spec.name).Functional,Results.(spec.name).Structural,MetricNames,spec.groupLevels,demohex,democon,demoorder,democheat, ...
            resolveSavePath(Settings.outputDir,"Aggregate Plots",strcat(spec.saveName,"_Combined")),strcat(spec.titlePrefix," (Functional vs Structural)"));
    end
end

% Subnetwork-aggregate violin grids: one comparison x one scope per figure,
% rows=subnetworks/columns=metrics (see plotSubnetAggregateViolinGrid).
% Functional and Structural are kept separate rather than combined side by
% side, since each scope has its own much smaller metric set at the
% subnetwork level. collectGridDataBySubnet re-keys Results by subnetwork
% name and naturally comes back empty if there's no subnetwork data this
% run (e.g. Settings.runFunctionalSubnets/runStructuralSubnets off).
for c=1:numel(comparisonGridSpecs)
    spec=comparisonGridSpecs(c);
    if ~Settings.tests.(spec.name).doGridPlot
        continue
    end
    bySubnetF=collectGridDataBySubnet(Results.(spec.name),Datasets,"Functional");
    if ~isempty(fieldnames(bySubnetF))
        plotSubnetAggregateViolinGrid(bySubnetF,SubnetLongNames,MetricNames,spec.groupLevels,demohex,democon,demoorder,democheat, ...
            resolveSavePath(Settings.outputDir,"Aggregate Plots",strcat(spec.saveName,"_Functional_AllSubnets")),strcat(spec.titlePrefix,", All Subnetworks (Functional)"));
    end
    bySubnetS=collectGridDataBySubnet(Results.(spec.name),Datasets,"Structural");
    if ~isempty(fieldnames(bySubnetS))
        plotSubnetAggregateViolinGrid(bySubnetS,SubnetLongNames,MetricNames,spec.groupLevels,demohex,democon,demoorder,democheat, ...
            resolveSavePath(Settings.outputDir,"Aggregate Plots",strcat(spec.saveName,"_Structural_AllSubnets")),strcat(spec.titlePrefix,", All Subnetworks (Structural)"));
    end
end

writetable(MuSigmaResults, strcat(Settings.outputDir,"MuSigmaResults.csv"));
writetable(EffectsResults, strcat(Settings.outputDir,"EffectsResults.csv"));

%% ===================== AGGREGATE FIGURES (1, 2, 3a, 3b) =====================
% Figures 1, 2, and 3a (below, commented out 2026-09-04) are superseded by
% the violin grids (Disorder/Sex/SexxDisorder/Severity) and the age-fit
% heatmap -- kept here, not deleted, in case the column-spec shape is
% useful reference later.
% Figure 1 columns: each 2-group comparison's main group-difference effect.
% Fig1Cols=struct( ...
%     'Comparison',{"Disorder","Sex","MCI","AD"}, ...
%     'Term',{"Impairment","Sex","Mild Impairment","Significant Impairment"}, ...
%     'GroupName',{"","","",""}, ...
%     'Label',{"Disorder","Sex","MCI","AD"});
%
% % Figure 2 columns: the 3 sex x impairment codings, 3 effect types each.
% Fig2Cols=struct( ...
%     'Comparison',{"SexxDisorder","SexxDisorder","SexxDisorder","SexxMCI","SexxMCI","SexxMCI","SexxAD","SexxAD","SexxAD"}, ...
%     'Term',{"Impairment","Sex","Interaction","Mild Impairment","Sex","Interaction","Significant Impairment","Sex","Interaction"}, ...
%     'GroupName',{"","","","","","","","",""}, ...
%     'Label',{"Disorder: Impairment","Disorder: Sex","Disorder: Interaction","MCI: Impairment","MCI: Sex","MCI: Interaction","AD: Impairment","AD: Sex","AD: Interaction"});
%
% % Figure 3a columns: per-subgroup age slopes, plus the two age-interaction terms.
% % (No "overall, ungrouped" column -- every age analysis in this pipeline is
% % already group-specific, so there's nothing to log for that.)
% Fig3aCols=struct( ...
%     'Comparison',{"Disorder","Disorder","MCI","AD","Sex","Sex","Disorder","Sex"}, ...
%     'Term',{"Age","Age","Age","Age","Age","Age","Interaction","Interaction"}, ...
%     'GroupName',{"CN","DS","MCI","AD","M","F","",""}, ...
%     'Label',{"CN","DS","MCI","AD","Male","Female","Age x Impairment","Age x Sex"});

% Age-fit heatmap columns: overall age effect standalone, then the
% Impairment x Age group (interaction + the two Disorder-comparison slopes)
% and the Sex x Age group (interaction + the two Sex-comparison slopes). No
% standalone sex main effect -- that's confounded with age-centering and is
% reported separately.
AgeHeatCols=struct( ...
    'Comparison',{"Disorder","Disorder","Disorder","Disorder","Sex","Sex","Sex"}, ...
    'Term',{"Age","Interaction","Age","Age","Interaction","Age","Age"}, ...
    'GroupName',{"","","CN","DS","","M","F"}, ...
    'Label',{"Age (Overall)","Interaction","Control (Age Slope)","Impaired (Age Slope)","Interaction","Male (Age Slope)","Female (Age Slope)"}, ...
    'GroupLabel',{"","Impairment x Age","Impairment x Age","Impairment x Age","Sex x Age","Sex x Age","Sex x Age"});

wholeBrainScopes=intersect(unique(EffectsResults.Subnetwork),["Functional","Structural"]);
scopesToRun=wholeBrainScopes;
if Settings.aggregateFigures.subnetFigures
    subnetScopes=setdiff(unique(EffectsResults.Subnetwork),["Functional","Structural"]);
    scopesToRun=[scopesToRun; subnetScopes];
end

if Settings.aggregateFigures.mainFigures || Settings.aggregateFigures.subnetFigures
    % The age-fit heatmap is produced once PER SCOPE (Functional and
    % Structural separately, not combined) -- so despite the historical
    % "Aggregate Figures" section header above (aggregated across
    % metrics/comparisons, not across scopes), it belongs in "Grid Plots
    % and Heat Maps" alongside the violin grids it plays the same role as,
    % not in "Aggregate Plots". Figures 1/2/3a (commented out 2026-09-04,
    % see Fig1Cols/Fig2Cols/Fig3aCols above) are no longer produced --
    % superseded by the violin grids and this heatmap.
    for s=1:numel(scopesToRun)
        scope=scopesToRun(s);
        plotAgeEffectHeatmap(EffectsResults,MuSigmaResults,scope,AgeHeatCols,MetricNames, ...
            resolveSavePath(Settings.outputDir,"Grid Plots and Heat Maps",strcat(scope,"/Figure_AgeFitHeatmap")),strcat("Age-Related Fits (",scope,")"));
    end

    % Combined Functional+Structural age-fit heatmap, stacked (whole-brain
    % only -- this is a cross-scope summary, not a per-subnet one, so it
    % goes in "Aggregate Plots" and follows mainFigures, not subnetFigures).
    if Settings.aggregateFigures.mainFigures && all(ismember(["Functional","Structural"],wholeBrainScopes))
        plotCombinedAgeEffectHeatmap(EffectsResults,MuSigmaResults,AgeHeatCols,MetricNames, ...
            resolveSavePath(Settings.outputDir,"Aggregate Plots","Figure_AgeFitHeatmap_Combined"),"Age-Related Fits (Functional vs Structural)");
    end

    % Subnetwork-aggregate age-fit heatmaps: one figure per scope (kept
    % separate, same reasoning as the subnetwork-aggregate violin grids),
    % rows=subnetworks stacked one metric-panel per figure. Follows
    % subnetFigures rather than mainFigures since it's specifically about
    % subnetworks, same as the per-subnetwork individual heatmaps above.
    if Settings.aggregateFigures.subnetFigures
        plotSubnetAggregateAgeHeatmap(EffectsResults,MuSigmaResults,"Functional",SubnetLongNames,AgeHeatCols,MetricNames, ...
            resolveSavePath(Settings.outputDir,"Aggregate Plots","Figure_AgeFitHeatmap_Functional_AllSubnets"),"Age-Related Fits, All Subnetworks (Functional)");
        plotSubnetAggregateAgeHeatmap(EffectsResults,MuSigmaResults,"Structural",SubnetLongNames,AgeHeatCols,MetricNames, ...
            resolveSavePath(Settings.outputDir,"Aggregate Plots","Figure_AgeFitHeatmap_Structural_AllSubnets"),"Age-Related Fits, All Subnetworks (Structural)");
    end

    % Figure 3b: auto-flag the strongest per-subgroup age relationships (from
    % the whole-brain scopes only) as candidates for individual scatter plots.
    % This doesn't draw anything new -- plotAgeTrendByGroup already saved a
    % PNG for every one of these; this just tells you which ones to feature.
    slopeRows=EffectsResults(ismember(EffectsResults.Subnetwork,wholeBrainScopes) ...
        & EffectsResults.Type=="slope" & EffectsResults.Group~="",:);
    slopeRows=sortrows(slopeRows,"PValue");
    nHeadline=min(Settings.aggregateFigures.numHeadlineScatters,height(slopeRows));
    fprintf("\nFigure 3b candidates (strongest per-group age relationships):\n");
    for k=1:nHeadline
        fprintf("  %s / %s / %s group: p=%5e, effect=%.2f  (look for the %s comparison's %s age-trend plot)\n", ...
            slopeRows.Subnetwork(k), slopeRows.Metric(k), slopeRows.Group(k), slopeRows.PValue(k), slopeRows.Effect(k), ...
            slopeRows.Comparison(k), slopeRows.Metric(k));
    end
end


%% ===================== LOCAL FUNCTIONS =====================

function Weights = inverseVarianceWeights(Error)
% inverseVarianceWeights  1./(Error.^2), except an exact-zero error maps
% to NaN instead of Inf (2026-09-05, replacing a too-blunt whole-row
% Phase-A drop that could wipe an entire scope's table when just one
% metric happened to be a structural placeholder -- value/error both
% zero for every subject -- confirmed with user). fitlme/fitlm already
% exclude any observation with a NaN in a predictor/response/weight, so
% this quietly excludes just the individual observations that lack a
% real error for THIS metric's own test, without touching T or any other
% metric's data.
    Weights=1./(Error.^2);
    Weights(Error==0)=NaN;
end

function [tStat, dfStat, p] = lmeContrastTest(lme, plusCoefs, minusCoefs)
% lmeContrastTest  Signed t-like test for a linear contrast (sum of
% plusCoefs' coefficients minus sum of minusCoefs') off an already-fitted
% LinearMixedModel -- the pooled-variance alternative to fitting a
% separate two-group model per comparison, used by every pairwise
% contrast in runWeightedOneWayAnova/runWeightedTwoWayAnova below
% (2026-09-05, per user direction to convert the ANOVA categories to
% mixed effects; confirmed as the more standard/textbook approach). Pass
% "" for a coefficient name to mean "the reference level" (no
% coefficient, contributes nothing to the contrast). coefTest returns an
% F-stat for a 1-row contrast (DF1=1), equivalent to a two-sided t-test
% with dfStat degrees of freedom; sign comes from the contrast's own
% estimate.
    names=lme.CoefficientNames;
    H=zeros(1,numel(names));
    for k=1:numel(plusCoefs)
        if strlength(plusCoefs(k))>0
            H(strcmp(names,plusCoefs(k)))=H(strcmp(names,plusCoefs(k)))+1;
        end
    end
    for k=1:numel(minusCoefs)
        if strlength(minusCoefs(k))>0
            H(strcmp(names,minusCoefs(k)))=H(strcmp(names,minusCoefs(k)))-1;
        end
    end
    estimate=H*lme.Coefficients.Estimate;
    [p,F,~,dfStat]=coefTest(lme,H);
    tStat=sign(estimate)*sqrt(F);
end

function result = runWeightedTwoWayAnova(FactorWeight, Weights, Factor1, Factor2, ID, factor1Levels, factor2Levels)
% runWeightedTwoWayAnova  Weighted 2x2 factorial mixed-effects test (main
% effects Factor1/Factor2 and their interaction, via fitlme+(1|ID)) plus
% the four pairwise contrasts needed for the grid's comparison brackets --
% replaces the old closed-form weighted-SS version (2026-09-05, per user
% direction to convert the ANOVA categories to mixed effects, matching
% the T-Test/Linear Fit categories already converted). Uses ONE
% pooled-variance model per comparison (not four separate two-group
% fits), with Bonferroni-corrected (x4) pairwise p-values -- confirmed
% with user as the more standard/textbook approach, needing less
% explanation of the statistical test than Tukey given the two-way
% designs' deliberately-partial 4-of-6 pairwise set. Same inverse-variance
% weighting as every other test in this file.
%
% Factor1/Factor2/ID are numeric/whole-group vectors over ALL rows in the
% comparison (not pre-split into groups); factor1Levels/factor2Levels name
% the 0/1 values of Factor1/Factor2, e.g. ["M","F"]/["CN","DS"] for
% SexxDisorder (index 1 = reference level). Cell layout: 0=(F1=ref,F2=ref),
% 1=(F1=ref,F2=other), 2=(F1=other,F2=ref), 3=(F1=other,F2=other) -- e.g.
% for SexxDisorder, 0=Male-Control, 1=Male-Diagnosed, 2=Female-Control,
% 3=Female-Diagnosed.
%
% Returns a struct with fields:
%   pAB, pCD, pI       - Factor1, Factor2, interaction omnibus p-values
%   FStatAB/CD/I, DFR  - their F-stats and shared residual DF (from anova(lme))
%   pA/pB/pC/pD        - four pairwise p-values (Bonferroni-corrected x4):
%                        A=0v1, B=2v3 (Factor2 effect within each Factor1
%                        level), C=0v2, D=1v3 (Factor1 effect within each
%                        Factor2 level)
%   tStatA/B/C/D, dfA/B/C/D - the matching contrast t-stats/DFs
    try
        F1=categorical(Factor1(:),[0,1],factor1Levels);
        F2=categorical(Factor2(:),[0,1],factor2Levels);
        tbl=table(FactorWeight(:),F1,F2,Weights(:),categorical(ID(:)), ...
            'VariableNames',{'FactorWeight','Factor1','Factor2','Weights','ID'});
        lme=fitlme(tbl,'FactorWeight ~ Factor1*Factor2 + (1|ID)','Weights',tbl.Weights);
        aovTbl=anova(lme);
        hitAB=find(strcmp(aovTbl.Term,'Factor1'),1);
        hitCD=find(strcmp(aovTbl.Term,'Factor2'),1);
        hitI=find(strcmp(aovTbl.Term,'Factor1:Factor2'),1);
        pAB=aovTbl.pValue(hitAB); FStatAB=aovTbl.FStat(hitAB); DFR=aovTbl.DF2(hitAB);
        pCD=aovTbl.pValue(hitCD); FStatCD=aovTbl.FStat(hitCD);
        pI=aovTbl.pValue(hitI); FStatI=aovTbl.FStat(hitI);

        coefF1=strcat("Factor1_",factor1Levels(2));
        coefF2=strcat("Factor2_",factor2Levels(2));
        coefI=strcat(coefF1,":",coefF2);
        [tA,dfA,pA]=lmeContrastTest(lme,"",coefF2);
        [tC,dfC,pC]=lmeContrastTest(lme,"",coefF1);
        [tB,dfB,pB]=lmeContrastTest(lme,"",[coefF2,coefI]);
        [tD,dfD,pD]=lmeContrastTest(lme,"",[coefF1,coefI]);

        result=struct('pAB',pAB,'pCD',pCD,'pI',pI,'FStatAB',FStatAB,'FStatCD',FStatCD,'FStatI',FStatI,'DFR',DFR, ...
            'pA',min(1,pA*4),'tStatA',tA,'dfA',dfA,'pB',min(1,pB*4),'tStatB',tB,'dfB',dfB, ...
            'pC',min(1,pC*4),'tStatC',tC,'dfC',dfC,'pD',min(1,pD*4),'tStatD',tD,'dfD',dfD);
    catch
        % Not enough usable data for this metric (e.g. every remaining
        % weight is NaN after inverseVarianceWeights excludes zero-error
        % rows) -- same "return NaN, let the caller's isnan(aov.pI) check
        % skip it" fallback runWeightedTTest already uses (2026-09-05).
        result=struct('pAB',NaN,'pCD',NaN,'pI',NaN,'FStatAB',NaN,'FStatCD',NaN,'FStatI',NaN,'DFR',NaN, ...
            'pA',NaN,'tStatA',NaN,'dfA',NaN,'pB',NaN,'tStatB',NaN,'dfB',NaN, ...
            'pC',NaN,'tStatC',NaN,'dfC',NaN,'pD',NaN,'tStatD',NaN,'dfD',NaN);
    end
end

function result = runWeightedOneWayAnova(FactorWeight, Weights, Group, ID, groupLevels)
% runWeightedOneWayAnova  Weighted one-way ANOVA mixed-effects test (via
% fitlme+(1|ID)) across any number of groups -- currently just Severity's
% Control/MCI/AD comparison. Replaces the old closed-form weighted-SS
% version (2026-09-05, per user direction to convert the ANOVA categories
% to mixed effects). Uses ONE pooled-variance model (not separate
% pairwise fits), with Bonferroni-corrected (x3) pairwise p-values in
% place of the old separate weightedPairwiseTTest calls. Same
% inverse-variance weighting as every other test in this file.
%
% Group/ID are numeric/whole-group vectors over ALL rows in the
% comparison (not pre-split into groups); Group takes values 0:(k-1),
% named by groupLevels (index 1 = reference level). The pairwise fields
% below assume exactly 3 groups, matching Severity's only current usage.
%
% Returns a struct with fields:
%   pOmnibus, FOmnibus, DF1Omnibus, DF2Omnibus - overall Group-effect F-test
%   pCM/pMA/pCA               - pairwise p-values (Bonferroni-corrected x3):
%                               group1 vs 2, group2 vs 3, group1 vs 3
%   tStatCM/MA/CA, dfCM/MA/CA - the matching contrast t-stats/DFs
    try
        GroupCat=categorical(Group(:),0:(numel(groupLevels)-1),groupLevels);
        tbl=table(FactorWeight(:),GroupCat,Weights(:),categorical(ID(:)), ...
            'VariableNames',{'FactorWeight','Group','Weights','ID'});
        lme=fitlme(tbl,'FactorWeight ~ Group + (1|ID)','Weights',tbl.Weights);
        aovTbl=anova(lme);
        hit=find(strcmp(aovTbl.Term,'Group'),1);
        pOmnibus=aovTbl.pValue(hit);
        FOmnibus=aovTbl.FStat(hit);
        DF1Omnibus=aovTbl.DF1(hit);
        DF2Omnibus=aovTbl.DF2(hit);

        coef2=strcat("Group_",groupLevels(2));
        coef3=strcat("Group_",groupLevels(3));
        [tCM,dfCM,pCM]=lmeContrastTest(lme,"",coef2);
        [tCA,dfCA,pCA]=lmeContrastTest(lme,"",coef3);
        [tMA,dfMA,pMA]=lmeContrastTest(lme,coef2,coef3);

        result=struct('pOmnibus',pOmnibus,'FOmnibus',FOmnibus,'DF1Omnibus',DF1Omnibus,'DF2Omnibus',DF2Omnibus, ...
            'pCM',min(1,pCM*3),'tStatCM',tCM,'dfCM',dfCM,'pMA',min(1,pMA*3),'tStatMA',tMA,'dfMA',dfMA, ...
            'pCA',min(1,pCA*3),'tStatCA',tCA,'dfCA',dfCA);
    catch
        % Not enough usable data for this metric -- same fallback as
        % runWeightedTwoWayAnova/runWeightedTTest (2026-09-05).
        result=struct('pOmnibus',NaN,'FOmnibus',NaN,'DF1Omnibus',NaN,'DF2Omnibus',NaN, ...
            'pCM',NaN,'tStatCM',NaN,'dfCM',NaN,'pMA',NaN,'tStatMA',NaN,'dfMA',NaN, ...
            'pCA',NaN,'tStatCA',NaN,'dfCA',NaN);
    end
end

function heightIn = violinRowHeightIn(nGroups)
% violinRowHeightIn  Per-row/panel physical height (inches) for the grid
% violin plots (plotFactorViolinGrid/plotCombinedViolinGrid/
% plotSubnetAggregateViolinGrid), as a function of group count (item 12,
% prior batch; exponent flattened further 2026-09-05 -- 3/4-group panels
% still looked a little squashed next to 2-group ones). The fixed 0.72in
% every row/panel used originally, regardless of nGroups, meant more
% groups packed into the same physical height compressed each violin's
% rendered width much more than fewer groups did (density ~ 1/nGroups,
% since the category-axis range grows with nGroups but the physical
% height didn't). This scales height so that compression is much less
% extreme (density ~ 1/nGroups^0.15 instead of ~1/nGroups), while
% reproducing the original 0.72in exactly at nGroups=2, so existing
% 2-group plots (Disorder, Sex, MCI, AD) don't change at all -- only the
% 3-4-group panels get noticeably more room.
    d0=0.3*2^0.15; % density (in per category-unit) at nGroups=2, matching
                   % the original fixed 0.72in/(2+0.4) exactly
    heightIn=d0*(nGroups+0.4)/nGroups^0.15;
end

function gridDataBySubnet = collectGridDataBySubnet(ComparisonResults, Datasets, setName)
% collectGridDataBySubnet  Filter one comparison's Results (keyed by full
% scopeField, e.g. "Functional_VN") down to just the subnetwork scopes for
% one Set ("Functional" or "Structural"), re-keyed by subnetwork name --
% the shape plotSubnetAggregateViolinGrid expects. Whole-brain scopes are
% excluded (Datasets.(scopeField).SubnetName=="" marks whole-brain).
    gridDataBySubnet=struct();
    scopeFields=fieldnames(ComparisonResults);
    for s=1:numel(scopeFields)
        sf=scopeFields{s};
        if ~isfield(Datasets,sf)
            continue
        end
        d=Datasets.(sf);
        if d.SetName==setName && strlength(d.SubnetName)>0
            gridDataBySubnet.(d.SubnetName)=ComparisonResults.(sf);
        end
    end
end

function [p, tStat, dfStat] = runWeightedTTest(V0, E0, ID0, V1, E1, ID1)
% runWeightedTTest  Weighted two-group test for one metric, via a
% mixed-effects model with a per-subject random intercept (fitlme,
% 'FactorWeight ~ Impairment + (1|ID)', weighted by inverse-variance).
% This is the same approach Disorder's T-Test already used; per user
% direction (2026-09-05), every T-Test comparison (Disorder/Sex/MCI/AD)
% now shares this one function instead of Sex/MCI/AD's previous
% closed-form weighted Welch test, for consistency.
%
%   V0, E0, ID0 - group 0's values/errors/subject-ID for this metric
%   V1, E1, ID1 - group 1's, same shape
%
% Returns NaN/NaN/NaN if the model can't be fit (not enough usable data).
    W0=inverseVarianceWeights(E0);
    W1=inverseVarianceWeights(E1);
    FactorWeight=[V0(:);V1(:)];
    Impairment=[zeros(numel(V0),1);ones(numel(V1),1)];
    Weights=[W0(:);W1(:)];
    ID_combined=[ID0(:);ID1(:)];
    tbl=table(FactorWeight,categorical(Impairment),Weights,categorical(ID_combined), ...
        'VariableNames',{'FactorWeight','Impairment','Weights','ID'});
    try
        lme=fitlme(tbl,'FactorWeight ~ Impairment + (1|ID)','Weights',Weights);
        p=lme.Coefficients.pValue(2);
        tStat=lme.Coefficients.tStat(2);
        dfStat=lme.Coefficients.DF(2);
    catch
        p=NaN;
        tStat=NaN;
        dfStat=NaN;
    end
end

function [V, R, Years, ID, Age, Sex, Diso] = subsetGroup(VAll, RAll, YearsAll, IDAll, AgeAll, SexAll, DisoAll, mask)
% subsetGroup  Slice every per-subject array down to one comparison
% group's rows by a logical mask -- the one operation every "XXX" group
% block in the group-subsetting section repeats seven times by hand
% (VSS/ESS/Years/ID/Age/Sex/Diso), factored out since copy-pasting all
% seven invites drift: two of the pre-refactor blocks (Years13/Years14)
% had exactly this, a stale mask left over from copy-pasting an earlier
% group and never caught since neither was ever actually read downstream.
    V=VAll(mask,:);
    R=RAll(mask,:);
    Years=YearsAll(mask);
    ID=IDAll(mask);
    Age=AgeAll(mask);
    Sex=SexAll(mask);
    Diso=DisoAll(mask);
end

function [Mu, Sigma] = weightedMuSigma(V, R)
% weightedMuSigma  Inverse-variance-weighted mean and standard error, per
% column, for a set of measurements V with known per-observation errors R
% (same size as V). Columns are independent; each row is one subject.
%
% This is the same measurement-error-aware pooling used throughout the
% analysis: a naive OLS variance estimate is computed and reduced by the
% average known measurement-error variance, clipped at zero, then used to
% inverse-variance-weight the mean.
    n=size(V,1);
    Mu=nan(1,size(V,2));
    Sigma=nan(1,size(V,2));
    for j=1:size(V,2)
        Vj=V(:,j);
        Rj=R(:,j);
        OLSVar=(1/(n-1))*sum((Vj-mean(Vj)).^2);
        MeasVar=Rj.^2;
        SigmaSqOLS=OLSVar-mean(MeasVar);
        if SigmaSqOLS<0
            SigmaSqOLS=0;
        end
        TotalVar=SigmaSqOLS+MeasVar;
        Weights=1./TotalVar;
        Mu(j)=sum(Vj.*Weights)/sum(Weights);
        Sigma(j)=1/sqrt(sum(Weights));
    end
end

function ResultsTable = appendMuSigmaRows(ResultsTable, Subnetwork, Comparison, GroupName, MetricNames, Mu, Sigma, N)
% appendMuSigmaRows  Append one row per metric to the running Mu/Sigma results table.
% Not every subnetwork/dataset determines the full set of named metrics, so
% Mu/Sigma can come back shorter than MetricNames -- align to whichever is shorter,
% and fall back to generic labels ("Metric15", ...) if Mu somehow has more entries
% than there are names for.
%
%   N - this group's row count (sample size diagnostic, added 2026-09-05):
%       same for every metric in this call, since a group is defined by
%       subject-level membership, not per-metric.
    n=numel(Mu);
    MetricNames=string(MetricNames(:));
    if numel(MetricNames)>=n
        MetricNames=MetricNames(1:n);
    else
        MetricNames=[MetricNames; ("Metric"+string((numel(MetricNames)+1):n))'];
    end
    NewRows=table(repmat(string(Subnetwork),n,1), repmat(string(Comparison),n,1), ...
        repmat(string(GroupName),n,1), MetricNames, Mu(:), Sigma(:), repmat(N,n,1), ...
        'VariableNames',{'Subnetwork','Comparison','Group','Metric','Mu','Sigma','N'});
    ResultsTable=[ResultsTable; NewRows];
end

function [MuSigmaResults, Mu, Sigma] = computeGroupMuSigma(V, R, TitleName1, Comparison, GroupName, MetricNames, MuSigmaResults, varSuffix, Efields, i, fid, logHeader)
% computeGroupMuSigma  Weighted mean/SE for one group -- the one calculation
% every "XXX MU/SIGMA" block in the Mu/Sigma section does, factored out
% since they were otherwise ~20-30 nearly-identical lines apiece.
%
%   V, R          - this group's raw values/errors (e.g. VSS0, ESS0)
%   varSuffix     - the numeric suffix used throughout the rest of the
%                   script for this group's Mu/Sigma (e.g. "0" for
%                   Mu0/Sigma0). Needed as a string, not just for naming:
%                   see below.
%   Efields, i, fid - passed through only so this function can replicate
%                   two things every original block did for you, by hand,
%                   every time:
%                   1) Besides returning Mu/Sigma as normal outputs (for
%                      this same iteration's use, e.g. Mu0(fac)), every
%                      block also stashed the value into the CALLER's
%                      workspace under a threshold-prefixed name (e.g.
%                      "II1Mu0"), via eval(strcat(Efields{i}(...),
%                      "Mu0=Mu0")). That's for a *different* purpose: some
%                      later code retrieves it back by that same prefixed
%                      name (still always for the SAME threshold i, not a
%                      different one -- every eval("Mu=",...) retrieval in
%                      this file uses Efields{i}, never Efields{l}). This
%                      function reproduces that exact storage via
%                      assignin('caller',...), the direct, non-eval way to
%                      write a variable into the calling workspace by a
%                      name built at runtime.
%                   2) The optional per-metric text log line (see
%                      logHeader below).
%   logHeader     - printed as a header line in the shared per-scope log,
%                   followed by one "Mu ± Sigma" line per metric. Pass ""
%                   to skip logging entirely (a few of the original blocks
%                   never had a log line either).
    [Mu,Sigma]=weightedMuSigma(V,R);
    MuSigmaResults=appendMuSigmaRows(MuSigmaResults,TitleName1,Comparison,GroupName,MetricNames,Mu,Sigma,size(V,1));

    prefix=Efields{i}(1:(strfind(Efields{i},"Thresh")-1));
    assignin('caller',strcat(prefix,"Mu",varSuffix),Mu);
    assignin('caller',strcat(prefix,"Sigma",varSuffix),Sigma);

    if strlength(logHeader)>0
        fprintf(fid,'%s\n',logHeader);
        formatSpec="%5e";
        for li=1:length(Mu)
            if li<length(Mu)
                fprintf(fid,strcat(formatSpec," ± "),Mu(li));
                fprintf(fid,strcat(formatSpec,","),Sigma(li));
            else
                fprintf(fid,strcat(formatSpec," ± "),Mu(li));
                fprintf(fid,strcat(formatSpec,"\n"),Sigma(li));
            end
        end
    end
end

function path = resolveSavePath(baseDir, category, scopeAndFile)
% resolveSavePath  Build (and ensure exists) the save path for one plot,
% centralizing the folder convention so every call site doesn't need its
% own mkdir/path-building logic:
%   baseDir\Individual Plots\<Set>\<Subnet>\<filename>
%   baseDir\Grid Plots and Heat Maps\<Set>\<Subnet>\<filename>
%   baseDir\Aggregate Plots\<filename>                (no Set/Subnet)
% <Subnet> is omitted entirely for a whole-brain scope.
%
%   category     - "Individual Plots", "Grid Plots and Heat Maps", or
%                  "Aggregate Plots"
%   scopeAndFile - for the two per-scope categories, "<TitleName1>/<filename>"
%                  (exactly what every comparison block already builds via
%                  callit=strcat(TitleName,"/",...), so call sites don't
%                  need to change how they build that string, only how
%                  they turn it into a path). TitleName1 is "Functional"/
%                  "Structural" alone for whole-brain, or "Functional
%                  <Subnet>" for a subnetwork. For "Aggregate Plots",
%                  scopeAndFile is just the filename (no scope, no "/").
    if category=="Aggregate Plots"
        dirPath=fullfile(baseDir,category);
        filename=scopeAndFile;
    else
        parts=strsplit(scopeAndFile,"/");
        titleName1=parts(1);
        filename=parts(2);
        setParts=strsplit(titleName1," ");
        setName=setParts(1);
        if numel(setParts)>1
            subnetName=strjoin(setParts(2:end)," ");
            dirPath=fullfile(baseDir,category,setName,subnetName);
        else
            dirPath=fullfile(baseDir,category,setName);
        end
    end
    if ~isfolder(dirPath)
        mkdir(dirPath);
    end
    path=fullfile(dirPath,filename);
end

function value = resolveToggle(Settings, path, defaultValue)
% resolveToggle  Look up a plot/test toggle by dotted path, most-specific
% wins: "the more specific one overrides the more general one" without
% every call site having to hand-write that precedence check itself.
%
%   path - dot-separated, BROADEST concept first, most SPECIFIC last, e.g.
%          "Individual.TTest.Impairment.Functional" (plot type -> test
%          family -> specific test -> scope). Checked from the full path
%          down to shorter and shorter prefixes -- "Individual.TTest.
%          Impairment.Functional", then "Individual.TTest.Impairment",
%          then "Individual.TTest", then "Individual" -- returning the
%          value of the FIRST one that's actually set. So setting only
%          Settings.toggles.Individual_TTest_Impairment_Functional = true
%          turns on just that one thing even if the broader
%          Settings.toggles.Individual = false says "no individual plots"
%          in general.
%   defaultValue - returned if nothing along the whole chain is set (i.e.
%          you've never touched the toggle system for this path at all).
%
% A toggle is "unset" simply by not being a field of Settings.toggles --
% there's no explicit tri-state (true/false/inherit) to manage, which
% means Settings.toggles only ever needs to list the exceptions, not the
% full ~80-toggle tree up front.
%
% Example: resolveToggle(Settings,"Individual.TTest.Impairment.Functional",true)
% with only Settings.toggles.Individual=false set returns false (nothing
% more specific overrides it); adding
% Settings.toggles.Individual_TTest_Impairment_Functional=true makes that
% same call return true, without touching Settings.toggles.Individual or
% any of the other Individual.* toggles it would otherwise still cover.
    segments=strsplit(path,".");
    for nSegments=numel(segments):-1:1
        fieldPath=strjoin(segments(1:nSegments),"_");
        if isfield(Settings.toggles,fieldPath)
            value=Settings.toggles.(fieldPath);
            return
        end
    end
    value=defaultValue;
end

function [TStat, DF, p] = weightedPairwiseTTest(Va, Wa, Vb, Wb)
% weightedPairwiseTTest  Inverse-variance-weighted, Welch-Satterthwaite
% two-sample t-test -- the same test already used inline throughout this
% file (e.g. the MCI and AD comparisons' own pairwise tests), factored out
% here since the Severity comparison needs it three times (Control-MCI,
% Control-AD, MCI-AD) for one figure.
%
% Note: the inline versions of this test elsewhere in the file use
% length(Vb)-1 in BOTH halves of the Welch-Satterthwaite denominator: this
% version uses each group's own length-1, which is the textbook formula.
    WMa=sum(Va.*Wa)/sum(Wa);
    WMb=sum(Vb.*Wb)/sum(Wb);
    WVa=sum(((Va-WMa).^2).*Wa)/sum(Wa);
    WVb=sum(((Vb-WMb).^2).*Wb)/sum(Wb);
    TStat=(WMa-WMb)/sqrt(WVa/length(Va)+WVb/length(Vb));
    DF=((WVa/length(Va)+WVb/length(Vb))^2) / ((WVa/length(Va))^2/(length(Va)-1) + (WVb/length(Vb))^2/(length(Vb)-1));
    p=2*(1-tcdf(abs(TStat),DF));
end

function model = fitAgeModel(tbl, formula, useMixedEffects)
% fitAgeModel  fitlme (with subject random intercept) or fitlm, both weighted
% by tbl.Weights, depending on useMixedEffects. Both return objects that
% expose .Coefficients.pValue and .CoefficientNames, which is all
% plotAgeTrendByGroup and its helpers rely on. (fitlm's LinearModel also
% exposes .ModelFitVsNullModel.Pvalue; fitlme's LinearMixedModel does not --
% see ageSlopePvalue for the model-vs-null p-value both types need.)
    if useMixedEffects
        model=fitlme(tbl,formula,'Weights',tbl.Weights);
    else
        model=fitlm(tbl,formula,'Weights',tbl.Weights);
    end
end

function p = ageSlopePvalue(model, useMixedEffects, tblk, isQuadratic)
% ageSlopePvalue  The model-vs-null p-value for one group's age-only refit
% (tblk ~ AgeCentered, or ~ AgeCentered + AgeCentered^2).
%
% Linear case: mathematically the same as the AgeCentered coefficient's own
% p-value (row 2), so that's used directly -- it works identically for both
% fitlm's LinearModel and fitlme's LinearMixedModel.
%
% Quadratic case ("does age matter at all, linear or quadratic"): a genuine
% joint test is needed. LinearModel exposes it directly via
% ModelFitVsNullModel; LinearMixedModel doesn't expose that property at all,
% so a null model (intercept + the same random-effects structure) is fit and
% compared via compare()'s likelihood-ratio test instead.
    if ~isQuadratic
        p=model.Coefficients.pValue(2);
        return
    end
    if useMixedEffects
        nullModel=fitlme(tblk,"FactorWeight ~ 1 + (1|ID)",'Weights',tblk.Weights);
        comparison=compare(nullModel,model);
        p=comparison.pValue(2);
    else
        p=model.ModelFitVsNullModel.Pvalue;
    end
end

function s = sigStars(p)
% sigStars  "*" / "**" / "***" for p<.05/.01/.001, "\dagger" for p<.1, else "".
    if isnan(p)
        s="";
    elseif p<0.001
        s="***";
    elseif p<0.01
        s="**";
    elseif p<0.05
        s="*";
    elseif p<0.1
        s=string(char(8224)); % dagger
    else
        s="";
    end
end

function s = boxSigLabel(p, baseFontSize)
% boxSigLabel  Same as sigStars, except: (1) every marker (stars AND a
% non-significant "n.s." fallback) is rendered in a softer sans-serif
% font (Calibri) than the box's own body text, reset back to the box's
% base font afterward via tex markup so it doesn't bleed into whatever
% text follows it in the same textbox (follow-up request, 2026-09-05);
% (2) a non-significant result (p>=0.1, where sigStars would return "")
% returns "n.s." instead, additionally shrunk and italicized within that
% same soft font. Used only in the summary/"box" text (plotAgeTrendByGroup's
% mainEffects box, plotGroupViolin's overallEffects box, the "Overall Fit"
% box) -- every other bracket/marker in the file keeps calling sigStars
% directly and is unaffected (daggers on marginal 0.05-0.1 results stay
% everywhere else). Boxes must set their base text 'FontName' to
% "Helvetica" (placeFitTextBox does) so this reset lands on a known font.
    stars=sigStars(p);
    softFont="Calibri";
    baseFont="Helvetica";
    if strlength(stars)==0
        smallSize=max(round(baseFontSize*0.7),7);
        % The space between \it and n.s. is just a tex command terminator,
        % not meant to be visible -- but the grid box (drawViolinGridPanel)
        % auto-wraps within a fixed-width textbox, and an ordinary space is
        % a valid wrap point to MATLAB's line-breaking, which could split
        % "n.s." away from the \it/\fontsize commands that style it
        % (2026-09-05, user report: wrapping was splitting a label from its
        % own marker at the space between them -- the same risk applies
        % here). A non-breaking space (char 160) still terminates the
        % command but is never treated as a wrap point.
        s=strcat("\fontname{",softFont,"}\fontsize{",string(smallSize),"}\it",string(char(160)),"n.s.\fontsize{",string(baseFontSize),"}\rm\fontname{",baseFont,"}");
    else
        s=strcat("\fontname{",softFont,"}",stars,"\fontname{",baseFont,"}");
    end
end

function textHandle = placeFitTextBox(x, y, str, fontSize, hAlign, vAlign)
% placeFitTextBox  Places a text label at (x,y) in the CURRENT axes' data
% coordinates with a known base FontName (Helvetica, so boxSigLabel's
% font-reset lands on a known font), no background or border (2026-09-05:
% scrapped after two attempts at an opaque background -- text()'s
% 'BackgroundColor' never rendered as a true opaque fill in this
% figure/export pipeline, and a manual patch behind the text introduced
% its own artifacts -- back to plain text, keeping the corner-placement/
% dynamic-clearance positioning work from those attempts). A forced
% drawnow keeps Extent reliable for callers that size axis limits around
% it. Returns the text handle so the caller can read its Extent.
    textHandle=text(x,y,str,'FontSize',fontSize,'FontName','Helvetica','Interpreter','tex', ...
        'HorizontalAlignment',hAlign,'VerticalAlignment',vAlign);
    drawnow;
end

function [ticks, niceStep] = niceTicksThroughZero(lo, hi, targetCount)
% niceTicksThroughZero  "Nice" evenly-spaced tick values (a round step --
% 1/2/5 x a power of ten -- like MATLAB's own auto-tick algorithm) that
% always include zero, anchored at zero instead of an arbitrary offset.
% Searches a small band of target counts (targetCount-1 .. targetCount+2)
% for a step that lands the actual tick count within [4,6], falling back
% to whichever candidate's count is closest to targetCount if none do --
% follow-up, 2026-09-05: a single fixed targetCount could swing wildly
% between e.g. 3 and 6 ticks depending on where lo/hi fell relative to a
% round step's boundaries, which read as "significant inconsistency"
% across plots even though each individually looked fine. Returns the
% chosen step too, so the caller can derive a matching decimal count
% (e.g. via decimalsForStep below) instead of a fixed format string.
    span=hi-lo;
    candSteps=zeros(1,4);
    candCounts=zeros(1,4);
    tcList=[targetCount-1,targetCount,targetCount+1,targetCount+2];
    for k=1:4
        tc=tcList(k);
        rawStep=span/max(tc-1,1);
        mag=10^floor(log10(rawStep));
        residual=rawStep/mag;
        if residual<1.5
            step=1*mag;
        elseif residual<3.5
            step=2*mag;
        elseif residual<7.5
            step=5*mag;
        else
            step=10*mag;
        end
        firstTick=ceil(lo/step)*step;
        candSteps(k)=step;
        candCounts(k)=numel(unique([0,firstTick:step:hi]));
    end
    inBand=find(candCounts>=4 & candCounts<=6,1);
    if ~isempty(inBand)
        pick=inBand;
    else
        [~,pick]=min(abs(candCounts-targetCount));
    end
    niceStep=candSteps(pick);
    firstTick=ceil(lo/niceStep)*niceStep;
    ticks=unique([0,firstTick:niceStep:hi]);
    ticks=ticks(ticks>=lo & ticks<=hi);
end

function d = decimalsForStep(step)
% decimalsForStep  Just enough decimal places to represent a "nice" tick
% step (1/2/5 x 10^k) exactly -- e.g. step=1 or 5 -> 0 decimals, step=0.5
% or 0.2 -> 1, step=0.05 -> 2 -- so a whole-number scale displays as
% whole numbers (was always forced to 2 decimals, so "10" showed as
% "10.00") while a fine-grained scale still gets the precision it needs.
% Combine with niceTicksThroughZero's returned step (2026-09-05).
    d=max(0,-floor(log10(step)+1e-9));
end

function label = groupDisplayLabel(code)
% groupDisplayLabel  Maps every group-level code used anywhere in this
% file (cfg.groupLevels entries) to the x-axis display string requested
% by the user (2026-09-05, item 2), for plotGroupViolin's individual
% violin plots. Falls back to the raw code itself for anything
% unrecognized, so a future new comparison doesn't silently lose its
% x-axis labels -- it just shows the bare code until a mapping is added
% here.
    switch code
        case "CN"
            label="Con.";
        case "DS"
            label="Impair.";
        case "MCI"
            label="Mild Impair.";
        case "AD"
            label="Sig. Impair.";
        case "M"
            label="Male";
        case "F"
            label="Female";
        case "MCN"
            label="Male Con.";
        case "MDS"
            label="Male Impair.";
        case "FCN"
            label="Female Con.";
        case "FDS"
            label="Female Impair.";
        case "MMCI"
            label="Male Mild Impair.";
        case "FMCI"
            label="Female Mild Impair.";
        case "MAD"
            label="Male Sig. Impair.";
        case "FAD"
            label="Female Sig. Impair.";
        otherwise
            label=code;
    end
end

function Effects = plotAgeTrendByGroup(tbl, cfg)
% plotAgeTrendByGroup  Scatter of every subject against Age, colored by group,
% with each group's fitted age trend line and a significance marker for that
% group's own slope. Handles any number of groups (2, 4, ...), mixed-effects
% (subject random intercept) or fixed-effects models, and linear-only or
% linear+quadratic age terms.
%
% tbl must contain: FactorWeight, Age, AgeCentered, Weights, GroupLabel
% (categorical, levels matching cfg.groupLevels), and ID if cfg.useMixedEffects.
%
% cfg fields:
%   formulaLinear, formulaQuadratic - model formula strings ('FactorWeight ~ ...'),
%                                     formulaQuadratic only used if includeQuadratic
%   includeQuadratic   - also fit the quadratic model and use it if it wins on AIC
%   useMixedEffects    - fitlme (subject random intercept) vs fitlm
%   groupLevels        - string array, plot order, matching tbl.GroupLabel and
%                         cfg.democheat entries
%   mainEffects        - struct array of (pattern,label,type): pattern is a
%                         regexp tested against the winning model's
%                         CoefficientNames, in the order the terms should
%                         appear in the plot's annotation and the text log.
%                         If the matched coefficient name contains "^2",
%                         "^2" is appended to the label automatically (so
%                         one "Age" entry covers both the linear and
%                         quadratic case). type is "diff" or "slope" (see
%                         tToStandardizedEffect) and controls how this term
%                         is converted to a standardized effect size.
%   demohex, democon, democheat - color lookup, as used throughout this file
%   titleStr           - plot title
%   markerSizeOffset   - scatter marker size floor ("plus" in the old code)
%   saveFile           - full output path, no extension
%   fid, logHeader     - text log file id and header line for this plot
%   drawPlot           - when false, skip every figure/save step and only
%                         fit the models + build+return the Effects table,
%                         so a test's result can be logged/reused regardless
%                         of whether this particular plot type is toggled on
%
% Returns a table with columns (Term, GroupName, TStat, DF, PValue, Type):
% one row per main effect (GroupName blank) and one row per group's own age
% slope (GroupName set). The caller converts these to standardized effects
% and logs them (via appendEffectRow) with whatever Subnetwork/Comparison/
% Metric context only the caller knows.

    Effects=table('Size',[0,6],'VariableTypes',{'string','string','double','double','double','string'}, ...
        'VariableNames',{'Term','GroupName','TStat','DF','PValue','Type'});

    try
        L=fitAgeModel(tbl,cfg.formulaLinear,cfg.useMixedEffects);
        if cfg.includeQuadratic
            Q=fitAgeModel(tbl,cfg.formulaQuadratic,cfg.useMixedEffects);
            models={L,Q};
            [~,which]=min([L.ModelCriterion.AIC,Q.ModelCriterion.AIC]);
        else
            models={L};
            which=1;
        end
    catch ME
        fprintf("Skipping %s (not enough usable data): %s\n",cfg.titleStr,ME.message);
        return
    end
    model=models{which};

    nGroups=numel(cfg.groupLevels);

    if cfg.drawPlot
        fig=figure('Visible','off');
        ax=gca; % captured once so the sig-box-only bonus save below can find every data child of THIS axes specifically
        set(ax,'box','off');
        set(fig,'units','inches','outerposition',[0 0 11 6.4],'windowstyle','normal') % item 4, 2026-09-05: a bit taller so the lifted title isn't clipped
        hold on

        if cfg.useMixedEffects
            tbl.Predicted=predict(model,tbl,'Conditional',false);
        else
            tbl.Predicted=predict(model,tbl);
        end

        lineHandles=gobjects(1,nGroups);
        for k=1:nGroups
            idx=tbl.GroupLabel==cfg.groupLevels(k);
            col=cfg.demohex(cfg.democheat==cfg.groupLevels(k));
            edgecol=cfg.democon(cfg.democheat==cfg.groupLevels(k),:);
            [sortedAge,sortIdx]=sort(tbl.Age(idx));
            pred=tbl.Predicted(idx);
            w=tbl.Weights(idx);
            scatter(tbl.Age(idx),tbl.FactorWeight(idx),30*normalize(w,'range')+cfg.markerSizeOffset,'filled', ...
                'MarkerFaceColor',col,'MarkerFaceAlpha',0.9,'MarkerEdgeColor',edgecol)
            lineHandles(k)=plot(sortedAge,pred(sortIdx),'-','LineWidth',2,'Color',col,'HandleVisibility','off');
        end
    end

    coefNames=model.CoefficientNames;
    pvals=model.Coefficients.pValue;
    ests=model.Coefficients.Estimate;
    ses=model.Coefficients.SE;
    mylabel="";
    % One companion text file per plot (same name, same folder, via
    % cfg.saveFile) instead of writing into the big shared per-scope log --
    % everything this specific plot's caption would need, in one place.
    if cfg.drawPlot
        txtFid=fopen(strcat(cfg.saveFile,".txt"),'wt');
        fprintf(txtFid,'%s\n',cfg.logHeader);
    end
    for m=1:numel(cfg.mainEffects)
        hit=find(~cellfun(@isempty,regexp(coefNames,cfg.mainEffects(m).pattern,'once')),1);
        if isempty(hit)
            continue
        end
        label=cfg.mainEffects(m).label;
        if contains(coefNames{hit},"^2")
            label=strcat(label,"^2");
        end
        p=pvals(hit);
        Effects(end+1,:)={label,"",ests(hit)/ses(hit),model.DFE,p,cfg.mainEffects(m).type}; %#ok<AGROW>
        if ~cfg.drawPlot
            continue
        end
        if strlength(mylabel)>0
            mylabel=strcat(mylabel,sprintf("\n"));
        end
        mylabel=strcat(mylabel,label," ",boxSigLabel(p,12));
        fprintf(txtFid,'%s',strcat(label,": p= "));
        fprintf(txtFid,strcat("%5e","\n"),p);
    end

    if cfg.drawPlot
        yline(0,'--','Color','k')
    end

    % Per-group age-effect significance: refit just that group's rows with the
    % same (linear or linear+quadratic) age terms as the winning overall model,
    % and use the model-vs-null p-value -- a joint test when quadratic, and
    % equivalent to the single coefficient's own p-value when linear.
    pGroup=nan(1,nGroups);
    for k=1:nGroups
        tblk=tbl(tbl.GroupLabel==cfg.groupLevels(k),:);
        try
            if cfg.useMixedEffects
                % fitlme errors out on a formula with no random-effect term
                % at all, so this needs the same (1|ID) the overall model
                % uses -- without it, every per-group refit here throws and
                % gets silently swallowed by the catch below, leaving that
                % group's slope missing everywhere it's used downstream.
                randomTerm=" + (1|ID)";
            else
                randomTerm="";
            end
            if which==1
                subModel=fitAgeModel(tblk,strcat("FactorWeight ~ AgeCentered",randomTerm),cfg.useMixedEffects);
            else
                subModel=fitAgeModel(tblk,strcat("FactorWeight ~ AgeCentered + AgeCentered^2",randomTerm),cfg.useMixedEffects);
            end
            pGroup(k)=ageSlopePvalue(subModel,cfg.useMixedEffects,tblk,which==2);
            % The "AgeCentered" coefficient is always the 2nd row (after the
            % intercept) regardless of linear or linear+quadratic, so this is
            % a consistent "slope near the mean age" for the forest plot even
            % when the plotted curve itself is quadratic.
            slopeT=subModel.Coefficients.Estimate(2)/subModel.Coefficients.SE(2);
            Effects(end+1,:)={"Age",cfg.groupLevels(k),slopeT,subModel.DFE,pGroup(k),"slope"}; %#ok<AGROW>
        catch ME
            pGroup(k)=NaN;
            fprintf("Skipping %s / %s group age slope (not enough usable data): %s\n",cfg.titleStr,cfg.groupLevels(k),ME.message);
        end
        if cfg.drawPlot
            fprintf(txtFid,'%s',strcat(cfg.groupLevels(k)," Age: p= "));
            fprintf(txtFid,strcat("%5e","\n"),pGroup(k));
        end
    end
    if cfg.drawPlot
        fclose(txtFid);
    end

    if ~cfg.drawPlot
        return
    end

    fontsize(14,'points')
    limx=xlim; limy=ylim;
    fracx=(limx(2)-limx(1))/20; fracy=(limy(2)-limy(1))/20;
    % Box pushed up and right, close to the true plot corner (item 8,
    % prior round: was 2*fracx/0.5*fracy in from the corner, now just a
    % sliver) -- placeFitTextBox handles the placement (no background,
    % per item 1 follow-up: opacity attempts scrapped).
    boxText=placeFitTextBox(limx(2)-0.2*fracx,limy(2)+3*fracy,mylabel,12,'right','top');
    ext=boxText.Extent; % read before hiding below, in case Visible='off' ever affects Extent computation
    boxText.Visible='off'; % 2026-09-05, user request: the box now lives only on the sig-box-only bonus image, turned back on right before that save
    % Extend the axes limits just enough to clear the box's own real
    % fit-to-text extent (whatever that turns out to be) plus a bit more
    % for a y-axis "x10^n" exponent label (item 17, prior batch), so the
    % box can sit close to the true plot edge without ever getting
    % clipped by it.
    xlim([limx(1),max(limx(2)+4*fracx,ext(1)+ext(3)+0.5*fracx)]);
    ylim([limy(1),max(limy(2)+4*fracy,ext(2)+ext(4)+0.5*fracy)]);
    % Item 5 (prior round): same "nice ticks through zero" treatment as the
    % individual violin plots' y-axis, computed off the real data range
    % (limy, before the box/exponent-label headroom was added above) so
    % ticks land on the actual data rather than the extra blank space.
    % Explicit tick label strings (item 2 follow-up, 2026-09-05) instead
    % of ytickformat+forced Exponent -- see the "Overall Fit" block for
    % why (that combination was unreliable both ways: showing a x10^4-scale
    % plot's decimal format applied to the scaled mantissa, e.g. "5.0000",
    % and rounding small-scale values to "0.00" elsewhere).
    [niceY,stepY]=niceTicksThroughZero(limy(1),limy(2),5);
    decY=decimalsForStep(stepY);
    yticks(niceY)
    yticklabels(arrayfun(@(v) sprintf(strcat('%.',num2str(decY),'f'),v), niceY, 'UniformOutput', false))
    xlabel("Age")
    ylabel("Factor Weight")
    t=title(cfg.titleStr);
    t.Units='normalized';
    t.Position(2)=t.Position(2)+0.02; % lift the title up slightly (item 5)
    set(gca,'box','off')

    % Per-group end-of-line significance marks (stars/dagger + a bounding
    % rectangle) removed 2026-09-05 per user direction (item 1) -- they
    % never rendered well and duplicated what the corner box above already
    % conveys via the Age/Interaction terms. Each group's own Age-slope
    % p-value (pGroup, computed above) is still written to the plot's
    % companion .txt file, just no longer drawn on the figure itself.

    saveas(fig,strcat(cfg.saveFile,".png"))

    % "Bonus" clean sig-box-only version (2026-09-05, user request,
    % same mechanism as plotGroupViolin): same figure, same axes/scaling/
    % ticks, but every scatter point/fitted trend line/zero-line hidden,
    % leaving only the mainEffects box against the correctly-scaled,
    % otherwise-empty axes -- for manually overlaying onto another
    % rendering of the same comparison later. Border framing just the
    % box text (not the axes' own outline, which with everything else
    % hidden just looks like one big rectangle connecting the axes to the
    % text inside it).
    dataChildren=ax.Children(ax.Children~=boxText);
    set(dataChildren,'Visible','off')
    set(boxText,'EdgeColor',[0.3,0.3,0.3],'Margin',4,'Visible','on') % turned back on here -- hidden on the main save above, per user request
    saveas(fig,strcat(cfg.saveFile,"_SigBoxOnly.png"))

    close(fig)
end

function h = clipViolinAtZero(h)
% clipViolinAtZero  If a violin's shape dips below zero, cut it off there.
% Only meaningful for the first metric plotted on a given axis (elsewhere in
% this file that's "fac==1", the "Rate" metric, which is why the original
% code only ever called this for that one case).
    if any(h.YData<=0)
        above=(h.XData>1).*(h.YData>=0)==1;
        below=(h.XData<1).*(h.YData>=0)==1;
        TheXData=[interp1(h.YData(h.XData>1),h.XData(h.XData>1),0);h.XData(above);h.XData(below);interp1(h.YData(h.XData<1),h.XData(h.XData<1),0)];
        TheYData=[0;h.YData(above);h.YData(below);0];
        h.XData=TheXData;
        h.YData=TheYData;
    end
end

function Effects = plotGroupViolin(groupValues, Mu, Sigma, comparisons, cfg)
% plotGroupViolin  One violin per group (subject-level distribution, with a
% weighted mean +/- SE marker), plus significance brackets for whichever
% group pairs are named in `comparisons`. Works for any number of groups.
%
%   groupValues - cell array {V1,V2,...}, one vector of subject-level values
%                 per group, in the same order as cfg.groupLevels
%   Mu, Sigma   - numeric arrays, one weighted mean/SE per group (e.g. from
%                 weightedMuSigma), same order as cfg.groupLevels
%   comparisons - struct array of (leftIdx,rightIdx,pValue,label,tStat,df):
%                 one significance bracket per entry. Non-significant
%                 entries (p>=0.1) draw nothing. Significant entries are
%                 packed into height levels via planComparisonLanes (same
%                 lane-packing the grid plots use): two brackets whose group
%                 spans don't overlap (e.g. "Male Impairment" 0v1 and
%                 "Female Impairment" 2v3) share the same level since they
%                 can never visually collide; anything that would collide
%                 gets its own level, stacked bottom-to-top. label is used
%                 in the text log and as the Term name in the returned
%                 Effects table; tStat/df are used only to compute that
%                 row's standardized effect size.
%
% cfg fields:
%   groupLevels                        - string array, plot order, matching
%                                         cfg.democheat entries
%   demohex, democon, demoorder, democheat - color lookup, as used throughout
%                                         this file
%   clipFirstAtZero                    - true for the "Rate" metric (fac==1),
%                                         matching the original code's
%                                         zero-floor clipping behavior
%   titleStr, saveFile, fid, logHeader
%   figureWidthInches                  - violin plots are narrow; this varies
%                                         by group count in the original code
%   overallEffects (optional)          - struct array of (pValue,label,tStat,df):
%                                         an extra text annotation for effects
%                                         that aren't a single pairwise
%                                         bracket (e.g. a 2x2 design's main
%                                         effects and interaction)
%   drawPlot                           - when false, skip every figure/save
%                                         step and only build+return the
%                                         Effects table, so a test's result
%                                         can be logged/reused regardless of
%                                         whether this particular plot type
%                                         is toggled on
% Whenever cfg.drawPlot is true, also saves a second "bonus" PNG (same
% file name + "_SigBoxOnly") with every violin/mean-line/errorbar/
% zero-line hidden, leaving only whatever significance indicator this
% comparison actually has, against the correctly-scaled, otherwise-empty
% axes: the overallEffects box for the ANOVA categories, or just the
% bracket+stars/dagger for a plain T-Test (which has no box) -- for
% manually overlaying onto another rendering of the same comparison later
% (2026-09-05, user request, extended from ANOVA-only to every individual
% violin plot).
%
% Returns a table with columns (Term, GroupName, TStat, DF, PValue, Type) --
% one row per entry in `comparisons` and, if given, per entry in
% cfg.overallEffects. GroupName is always blank (every effect here is a
% between-group difference, not a per-group value); Type is always "diff".
% The caller converts these to standardized effects and logs them (via
% appendEffectRow) with whatever Subnetwork/Comparison/Metric context only
% the caller knows.

    Effects=table('Size',[0,6],'VariableTypes',{'string','string','double','double','double','string'}, ...
        'VariableNames',{'Term','GroupName','TStat','DF','PValue','Type'});

    if cfg.drawPlot
        fig=figure('Visible','off');
        ax=gca; % captured once so the sig-box-only pass below can find every data child of THIS axes specifically
        sigHandles=gobjects(0); % every bracket/star/dagger/box handle -- kept visible in the sig-box-only bonus save below; everything else in ax.Children is "data" and gets hidden there
        boxTextHandle=gobjects(0,1); % stays empty for T-Test plots, which have no overallEffects box
        set(ax,'box','off');
        set(fig,'units','inches','outerposition',[0 0 cfg.figureWidthInches 12],'windowstyle','normal')
        hold on

        nGroups=numel(cfg.groupLevels);
        maxes=zeros(1,nGroups);
        mins=zeros(1,nGroups);
        for k=1:nGroups
            lvl=cfg.groupLevels(k);
            [h,vLine,~,~]=violin({groupValues{k}},'mc','','medc','', ...
                'facecolor',cfg.demoorder(cfg.democheat==lvl,:),'facealpha',0.75,'edgecolor',cfg.democon(cfg.democheat==lvl,:));
            vLine.Visible='off';
            if cfg.clipFirstAtZero
                h=clipViolinAtZero(h);
            end
            h.XData=h.XData+(k-1);
            plot([interp1(h.YData(h.XData<k),h.XData(h.XData<k),Mu(k)),interp1(h.YData(h.XData>k),h.XData(h.XData>k),Mu(k))], ...
                [Mu(k),Mu(k)],'Color',cfg.democon(cfg.democheat==lvl,:),'LineWidth',2)
            errorbar(k,Mu(k),Sigma(k),'Color',cfg.democon(cfg.democheat==lvl,:))
            maxes(k)=max(h.YData);
            mins(k)=min(h.YData);
        end

        % Extra 0.7 units of buffer at the right (item 7, 2026-09-05), only
        % when there's an overallEffects box that needs the room -- plain
        % pairwise T-Test plots have no box, so they keep the original
        % tight layout. figureWidthInches is widened to match at each
        % call site that has a box, so this doesn't just compress the
        % violins.
        if isfield(cfg,'overallEffects') && ~isempty(cfg.overallEffects)
            xlim([0.5,nGroups+1.2])
        else
            xlim([0.5,nGroups+0.5])
        end
        xticks(1:nGroups)
        xticklabels(arrayfun(@groupDisplayLabel,cfg.groupLevels)) % item 2, 2026-09-05
        xtickangle(30) % item 3, 2026-09-05: always diagonal
        lim=[min(mins),max(maxes)];
        frac=(lim(2)-lim(1))/15;

        % Lane-pack the brackets into height levels (item 8, 2026-09-05):
        % comparisons whose group spans don't overlap (e.g. "Male
        % Impairment" 0v1 and "Female Impairment" 2v3) share the same
        % level instead of stacking needlessly on top of each other --
        % the same interval-packing planComparisonLanes already uses for
        % the grid plots' bracket columns, reused here for vertical height
        % levels instead of horizontal lanes. Levels with no significant
        % comparison contribute no gap (only active lanes get a level).
        [laneOfComparison,nLanesActual]=planComparisonLanes(comparisons);
        activeMask=arrayfun(@(c) ~isnan(c.pValue) && c.pValue<0.1, comparisons);
        activeLaneNums=unique(laneOfComparison(activeMask));
        nLevels=numel(activeLaneNums);
        levelOfLane=zeros(1,max(nLanesActual,1));
        levelOfLane(activeLaneNums)=0:(nLevels-1);

        % Item 7: if the zero line would land inside the bracket zone
        % (e.g. mostly-negative data, where lim(2) sits just above 0),
        % push the whole bracket zone up an extra couple of frac so the
        % lowest bracket doesn't visually collide with yline(0).
        bracketBase=lim(2);
        if 0>lim(2)-frac && 0<lim(2)+2*frac
            bracketBase=bracketBase+2*frac;
        end
        ylim([lim(1),bracketBase+3*max(1,nLevels)*frac])

        % One companion text file per plot (same name, same folder, via
        % cfg.saveFile) instead of writing into the big shared per-scope log --
        % everything this specific plot's caption would need, in one place.
        txtFid=fopen(strcat(cfg.saveFile,".txt"),'wt');
        fprintf(txtFid,'%s\n',cfg.logHeader);
    end

    for c=1:numel(comparisons)
        p=comparisons(c).pValue;
        Effects(end+1,:)={comparisons(c).label,"",comparisons(c).tStat,comparisons(c).df,p,"diff"}; %#ok<AGROW>
        if ~cfg.drawPlot
            continue
        end
        x1=comparisons(c).leftIdx;
        x2=comparisons(c).rightIdx;
        fprintf(txtFid,'%s',strcat(comparisons(c).label,": p= "));
        fprintf(txtFid,strcat("%5e","\n"),p);
        if isnan(p) || p>=0.1
            continue
        end
        levelIdx=levelOfLane(laneOfComparison(c));
        y=bracketBase+frac+3*frac*levelIdx;
        tickLen=0.4*frac;
        % Drawn as an actual bracket (a horizontal spine plus two short
        % vertical end-ticks pointing down at the plot) rather than a
        % plain line, matching drawViolinGridPanel's bracket style
        % (item 8, 2026-09-05). Handles kept in sigHandles (2026-09-05
        % follow-up) so the sig-box-only bonus save below can leave the
        % bracket itself visible -- for a plain T-Test plot (no
        % overallEffects box), the bracket+stars IS the significance
        % indicator, so it's what stays instead of a box.
        hSpine=plot([x1,x2],[y,y],'Color','black','LineWidth',1.5);
        hTick1=plot([x1,x1],[y-tickLen,y],'Color','black','LineWidth',1.5);
        hTick2=plot([x2,x2],[y-tickLen,y],'Color','black','LineWidth',1.5);
        sigHandles=[sigHandles,hSpine,hTick1,hTick2]; %#ok<AGROW>
        xm=(x1+x2)/2;
        starY=y+0.5*frac;
        if p<0.05
            hStar1=plot(xm,starY,'k*');
            sigHandles=[sigHandles,hStar1]; %#ok<AGROW>
            if p<0.01
                hStar2=plot(xm-0.15,starY,'k*');
                sigHandles=[sigHandles,hStar2]; %#ok<AGROW>
                if p<0.001
                    hStar3=plot(xm+0.15,starY,'k*');
                    sigHandles=[sigHandles,hStar3]; %#ok<AGROW>
                end
            end
        else
            hDagger=text(xm,starY,char(8224),'FontSize',14,'Color','k','HorizontalAlignment','center');
            sigHandles=[sigHandles,hDagger]; %#ok<AGROW>
        end
    end

    if isfield(cfg,'overallEffects') && ~isempty(cfg.overallEffects)
        mylabel="";
        for m=1:numel(cfg.overallEffects)
            Effects(end+1,:)={cfg.overallEffects(m).label,"",cfg.overallEffects(m).tStat,cfg.overallEffects(m).df,cfg.overallEffects(m).pValue,"diff"}; %#ok<AGROW>
            if ~cfg.drawPlot
                continue
            end
            if m>1
                mylabel=strcat(mylabel,sprintf("\n"));
            end
            mylabel=strcat(mylabel,cfg.overallEffects(m).label," ",boxSigLabel(cfg.overallEffects(m).pValue,14));
            fprintf(txtFid,'%s',strcat(cfg.overallEffects(m).label,": p= "));
            fprintf(txtFid,strcat("%5e","\n"),cfg.overallEffects(m).pValue);
        end
        if cfg.drawPlot
            % Item 3/7 (prior round): positioned in DATA coordinates within
            % the xlim buffer reserved above (nGroups+1.2 upper limit) so
            % it can't cover the rightmost bracket, near the top of the
            % plot (rather than nGroups/2+0.7, which ran past the right
            % edge on narrow figures like Severity's 3-group/3in-wide or
            % SexxAD's 4-group/4in-wide plots). No background (item 1
            % follow-up, 2026-09-05: opacity attempts scrapped).
            yl=ylim;
            boxTextHandle=placeFitTextBox(nGroups+1.15,yl(2),mylabel,14,'right','top');
            boxTextHandle.Visible='off'; % 2026-09-05, user request: the box now lives only on the sig-box-only bonus image, turned back on right before that save below
            sigHandles=[sigHandles,boxTextHandle];
        end
    end

    if cfg.drawPlot
        fclose(txtFid);
        ylabel("Factor Weight")
        fontsize(14,'points')
        yl=ylim;
        % Item 6 (prior round), refined further: "nice" ticks that always
        % include 0, at a round step size close to a 5-tick target
        % (searched within a small band to avoid wildly different counts
        % across plots); decimal places derived from the step itself
        % (whole-number steps show as whole numbers, was always forced to
        % 2 decimals). Explicit tick label strings (item 2 follow-up,
        % 2026-09-05) instead of ytickformat+forced Exponent -- see the
        % "Overall Fit" block for why (that combination was unreliable
        % both ways: showing a x10^4-scale plot's decimal format applied
        % to the scaled mantissa, e.g. "5.0000", and rounding small-scale
        % values to "0.00" elsewhere).
        [niceY,stepY]=niceTicksThroughZero(yl(1),yl(2),5);
        decY=decimalsForStep(stepY);
        yticks(niceY)
        yticklabels(arrayfun(@(v) sprintf(strcat('%.',num2str(decY),'f'),v), niceY, 'UniformOutput', false))
        t=title(cfg.titleStr,'FontSize',19);
        t.Units='normalized';
        t.Position(2)=t.Position(2)+0.02; % item 5, 2026-09-05: lift title up slightly
        yline(0,'--','Color','k')
        set(gca,'box','off');

        saveas(fig,strcat(cfg.saveFile,".png"))

        % "Bonus" clean sig-box-only version (2026-09-05, extended to every
        % individual violin plot, T-Tests included): same figure, same
        % axes/scaling/ticks, but every violin/mean-line/errorbar/
        % zero-line hidden -- leaving only whatever significance
        % indicator this comparison actually has: the overallEffects box
        % for the ANOVA categories, or just the bracket+stars/dagger for a
        % plain T-Test (which has no box at all). sigHandles tracks
        % exactly those "keep" elements as they're drawn above; everything
        % else in ax.Children is "data" and gets hidden here. Border
        % framing just the significance TEXT (when there is one) rather
        % than set(ax,'box','on'), which drew the axes' own outline and,
        % with everything else hidden, just looked like one big rectangle
        % connecting the axes to the text inside it. text()'s own
        % EdgeColor/Margin are computed by MATLAB's renderer at draw time
        % (unlike our own earlier manual-patch attempt, which read a
        % possibly-stale Extent before the text had been laid out), so
        % this doesn't carry the same risk.
        dataChildren=ax.Children(~ismember(ax.Children,sigHandles));
        set(dataChildren,'Visible','off')
        if ~isempty(boxTextHandle)
            set(boxTextHandle,'EdgeColor',[0.3,0.3,0.3],'Margin',4,'Visible','on') % turned back on here -- hidden on the main save above, per user request
        end
        saveas(fig,strcat(cfg.saveFile,"_SigBoxOnly.png"))

        close(fig)
    end
end

function [effect, ciLow, ciHigh] = tToStandardizedEffect(t, df, type)
% tToStandardizedEffect  Convert a t-statistic + degrees of freedom into a
% standardized, scale-free effect size with an approximate 95% CI, so effects
% across metrics with very different raw scales (e.g. Rate ~100 vs a factor
% weight ~0.1) can sit on one forest-plot axis.
%
%   type "diff"  - a two-group difference (main effect, pairwise comparison,
%                  or an F-test with 1 numerator df converted via t=sqrt(F)):
%                  Cohen's d = 2t/sqrt(df), the standard t-to-d relationship.
%   type "slope" - a continuous predictor's coefficient (an age slope):
%                  standardized partial correlation r = t/sqrt(t^2+df).
%
% Returns NaN/NaN/NaN if t or df isn't usable (e.g. a skipped fit).
    if isnan(t) || isnan(df) || df<=0
        effect=NaN; ciLow=NaN; ciHigh=NaN;
        return
    end
    if type=="slope"
        effect=t/sqrt(t^2+df);
        se=sqrt((1-effect^2)/df);
    else
        effect=2*t/sqrt(df);
        se=sqrt(4/df+effect^2/(2*df));
    end
    ciLow=effect-1.96*se;
    ciHigh=effect+1.96*se;
end

function EffectsResults = appendEffectRow(EffectsResults, Subnetwork, Comparison, Term, GroupName, Metric, TStat, DF, PValue, Type)
% appendEffectRow  Convert one (t, df) pair to a standardized effect + CI and
% append it as one row to the running effects results table. The raw TStat
% and DF are kept alongside the standardized Effect -- Figures 1/2/3a plot
% the standardized version (comparable across metrics with very different
% raw scales), but some figures (e.g. the age-fit heatmap) want the raw
% t-statistic itself.
    [effect,ciLow,ciHigh]=tToStandardizedEffect(TStat,DF,Type);
    NewRow=table(string(Subnetwork), string(Comparison), string(Term), string(GroupName), string(Metric), ...
        TStat, DF, effect, ciLow, ciHigh, PValue, string(Type), ...
        'VariableNames',{'Subnetwork','Comparison','Term','Group','Metric','TStat','DF','Effect','CILow','CIHigh','PValue','Type'});
    EffectsResults=[EffectsResults; NewRow];
end

function cmap = divergingColormap()
% divergingColormap  A simple blue-white-red diverging colormap for the
% effect-size heatmap, written out directly rather than relying on
% redbluecmap (Bioinformatics Toolbox, not guaranteed to be installed).
    n=128;
    cmap=[linspace(0.15,1,n)',linspace(0.15,1,n)',ones(n,1); ...
          ones(n,1),linspace(1,0.15,n)',linspace(1,0.15,n)'];
end

function cmap = neutralDivergingColormap()
% neutralDivergingColormap  Purple-white-orange diverging colormap. Blue/red
% already means Control/Impaired and green/pink already means Sex elsewhere
% in this file, so the age-fit heatmap (which mixes both kinds of terms in
% one figure) uses this neutral third option instead.
    n=128;
    negColor=[0.90,0.57,0.05];
    posColor=[0.46,0.23,0.63];
    cmap=[linspace(negColor(1),1,n)',linspace(negColor(2),1,n)',linspace(negColor(3),1,n)'; ...
          linspace(1,posColor(1),n)',linspace(1,posColor(2),n)',linspace(1,posColor(3),n)'];
end

function [Effect, CILow, CIHigh, PValue, MetricOrder, TStat] = gatherEffectGrid(EffectsResults, subnetworkList, columnSpecs, MetricOrder)
% gatherEffectGrid  Build a Metric x Column grid of standardized effects
% (with CI, p-value and raw t-statistic grids the same size) from
% EffectsResults, for whichever (Comparison,Term,GroupName) triple each
% columnSpec names. Subnetwork rows are pooled across everything in
% subnetworkList (normally just one name). Rows where every column is NaN
% are dropped, so a scope that never determined a given metric (e.g.
% Assortativity for whole-brain) doesn't leave a blank row in the figure.
    nCols=numel(columnSpecs);
    nRows=numel(MetricOrder);
    Effect=nan(nRows,nCols);
    CILow=nan(nRows,nCols);
    CIHigh=nan(nRows,nCols);
    PValue=nan(nRows,nCols);
    TStat=nan(nRows,nCols);
    inScope=ismember(EffectsResults.Subnetwork,subnetworkList);
    for c=1:nCols
        inCol=inScope & EffectsResults.Comparison==columnSpecs(c).Comparison ...
            & EffectsResults.Term==columnSpecs(c).Term & EffectsResults.Group==columnSpecs(c).GroupName;
        colRows=EffectsResults(inCol,:);
        for m=1:nRows
            hit=find(colRows.Metric==MetricOrder(m),1);
            if ~isempty(hit) && ~isnan(colRows.Effect(hit))
                Effect(m,c)=colRows.Effect(hit);
                CILow(m,c)=colRows.CILow(hit);
                CIHigh(m,c)=colRows.CIHigh(hit);
                PValue(m,c)=colRows.PValue(hit);
                TStat(m,c)=colRows.TStat(hit);
            end
        end
    end
    keepRows=any(~isnan(Effect),2);
    Effect=Effect(keepRows,:);
    CILow=CILow(keepRows,:);
    CIHigh=CIHigh(keepRows,:);
    PValue=PValue(keepRows,:);
    TStat=TStat(keepRows,:);
    MetricOrder=MetricOrder(keepRows);
end

function plotForestSummary(EffectsResults, subnetworkList, columnSpecs, MetricOrder, saveFile, titleStr)
% plotForestSummary  One row per metric, one dot+CI per columnSpec, all on a
% standardized-effect axis so metrics with very different raw scales (Rate
% ~100 vs a factor weight ~0.1) sit on one comparable plot. A filled marker
% means p<0.05; an open marker means the CI is shown but the effect isn't
% significant at that threshold. Rows with no data in any column are
% dropped entirely (e.g. metrics that don't apply at this scope).
%
%   subnetworkList - string array; EffectsResults.Subnetwork values to pool
%                    (normally one name, e.g. "Functional" or "Functional VN")
%   columnSpecs    - struct array of (Comparison,Term,GroupName,Label): one
%                    column/series per entry, drawn left-to-right in legend
%                    order and as a small vertical offset within each row
%   MetricOrder    - string array giving the canonical row order (e.g.
%                    MetricNames): "Rate" then the named factors

    [Effect,CILow,CIHigh,PValue,MetricOrder]=gatherEffectGrid(EffectsResults,subnetworkList,columnSpecs,MetricOrder);
    nMetrics=numel(MetricOrder);
    nCols=numel(columnSpecs);
    if nMetrics==0
        fprintf("Skipping forest plot '%s': no data for any metric in this scope.\n",titleStr);
        return
    end

    fig=figure('Visible','off');
    set(fig,'units','inches','outerposition',[0 0 9 1.2+0.5*nMetrics],'windowstyle','normal')
    hold on

    colColors=lines(nCols);
    yOffsets=linspace(-0.3,0.3,nCols);
    if nCols==1
        yOffsets=0;
    end
    for m=1:nMetrics
        for c=1:nCols
            if isnan(Effect(m,c))
                continue
            end
            y=nMetrics-m+1+yOffsets(c);
            plot([CILow(m,c),CIHigh(m,c)],[y,y],'-','Color',colColors(c,:),'LineWidth',1.2,'HandleVisibility','off')
            if PValue(m,c)<0.05
                markerFace=colColors(c,:);
            else
                markerFace='none';
            end
            plot(Effect(m,c),y,'o','MarkerEdgeColor',colColors(c,:),'MarkerFaceColor',markerFace, ...
                'MarkerSize',6,'LineWidth',1.2,'HandleVisibility','off')
        end
    end

    xline(0,'--','Color',[0.5 0.5 0.5],'HandleVisibility','off')
    yticks(1:nMetrics)
    yticklabels(flip(MetricOrder))
    ylim([0.5,nMetrics+0.5])
    xlabel("Standardized Effect")
    title(titleStr,'Interpreter','none')
    box off

    for c=1:nCols
        plot(NaN,NaN,'o-','Color',colColors(c,:),'MarkerFaceColor',colColors(c,:),'DisplayName',columnSpecs(c).Label);
    end
    legend('Location','eastoutside')
    fontsize(11,'points')

    saveas(fig,strcat(saveFile,".png"))
    close(fig)
end

function plotEffectHeatmap(EffectsResults, subnetworkList, columnSpecs, MetricOrder, saveFile, titleStr)
% plotEffectHeatmap  Metric x (comparison,term) grid, color = standardized
% effect size. Rows with no data in any column are dropped. Cells show the
% numeric effect size; significance isn't overlaid on the cell (MATLAB's
% heatmap doesn't support a second text layer cleanly) -- check EffectsResults
% or the p-value grid directly for that.
    [Effect,~,~,~,MetricOrder]=gatherEffectGrid(EffectsResults,subnetworkList,columnSpecs,MetricOrder);
    nMetrics=numel(MetricOrder);
    nCols=numel(columnSpecs);
    if nMetrics==0
        fprintf("Skipping heatmap '%s': no data for any metric in this scope.\n",titleStr);
        return
    end

    fig=figure('Visible','off');
    set(fig,'units','inches','outerposition',[0 0 2+nCols*0.9 1.2+0.4*nMetrics],'windowstyle','normal')

    h=heatmap({columnSpecs.Label},cellstr(MetricOrder),Effect);
    h.Title=titleStr;
    h.Colormap=divergingColormap();
    lim=max(abs(Effect(:)),[],'omitnan');
    if isnan(lim) || lim==0
        lim=1;
    end
    h.ColorLimits=[-lim,lim];
    h.CellLabelFormat='%.2f';
    h.MissingDataColor=[0.92 0.92 0.92];

    saveas(fig,strcat(saveFile,".png"))
    close(fig)
end

function [DisplayVal, PValue, MetricOrder] = gatherAgeHeatmapGrid(EffectsResults, MuSigmaResults, subnetworkList, columnSpecs, MetricOrder)
% gatherAgeHeatmapGrid  Metric x age-related-term grid of DISPLAY values for
% the age-fit heatmap: t * sign(Mu) -- Mu being that row's own metric's
% overall (whole-sample, pooled) baseline mean from MuSigmaResults -- not
% the raw t-statistic. So positive/color-A always means "moving away from
% zero" and negative/color-B always means "moving toward zero" for that
% metric, consistently across every column in the row, rather than the raw
% coefficient sign (increasing vs decreasing), which alone doesn't say
% whether a negative slope is crossing toward/through zero or just getting
% more negative from an already-negative baseline. Magnitude and
% significance (sigStars(), from the un-flipped p-value) are unaffected --
% only the sign shown/colored changes. Rows this scope has no data for at
% all are dropped (see gatherEffectGrid), so MetricOrder returned here may
% be a subset of the one passed in.
    [~,~,~,PValue,MetricOrder,TStat]=gatherEffectGrid(EffectsResults,subnetworkList,columnSpecs,MetricOrder);
    nMetrics=numel(MetricOrder);

    % One baseline sign per metric (row), from the overall (whole-sample)
    % pooled mean already logged elsewhere as Comparison="Overall"/
    % Group="All", applied uniformly to every column in that row.
    rowMuSign=ones(nMetrics,1);
    muRows=MuSigmaResults(ismember(MuSigmaResults.Subnetwork,subnetworkList) ...
        & MuSigmaResults.Comparison=="Overall" & MuSigmaResults.Group=="All",:);
    for m=1:nMetrics
        hit=find(muRows.Metric==MetricOrder(m),1);
        if ~isempty(hit) && ~isnan(muRows.Mu(hit)) && muRows.Mu(hit)~=0
            rowMuSign(m)=sign(muRows.Mu(hit));
        end
    end
    DisplayVal=TStat.*rowMuSign;
end

function plotAgeEffectHeatmap(EffectsResults, MuSigmaResults, subnetworkList, columnSpecs, MetricOrder, saveFile, titleStr)
% plotAgeEffectHeatmap  Metric x age-related-term grid, one scope
% (Functional, Structural, or a subnetwork) per figure. See
% gatherAgeHeatmapGrid for what's plotted/colored in each cell.
%
% Columns are grouped under bracketed headers wherever consecutive
% columnSpecs share a non-blank GroupLabel (e.g. the interaction + two
% per-group slopes that make up "Impairment x Age"); a columnSpec with a
% blank GroupLabel stands alone with no bracket (e.g. "Age (Overall)").
% columnSpecs otherwise has the same fields as gatherEffectGrid expects
% (Comparison, Term, GroupName, Label), plus this GroupLabel.
    [DisplayVal,PValue,MetricOrder]=gatherAgeHeatmapGrid(EffectsResults,MuSigmaResults,subnetworkList,columnSpecs,MetricOrder);
    nMetrics=numel(MetricOrder);
    nCols=numel(columnSpecs);
    if nMetrics==0
        fprintf("Skipping age-fit heatmap '%s': no data for any metric in this scope.\n",titleStr);
        return
    end

    % One companion text file per heatmap image (2026-09-05, user request),
    % same name/folder as the .png via saveFile: every cell's display value
    % (t*sign(Mu), same quantity the color encodes) and p-value, organized
    % by metric row, so exact numbers are available without eyeballing
    % color intensity.
    txtFid=fopen(strcat(saveFile,".txt"),'wt');
    fprintf(txtFid,'%s\n',titleStr);
    for m=1:nMetrics
        fprintf(txtFid,'%s\n',MetricOrder(m));
        writeHeatmapCellLines(txtFid,"  ",columnSpecs,DisplayVal,PValue,m);
    end
    fclose(txtFid);

    cellFontSize=8;
    labelFontSize=9;

    colWidthIn=0.95;
    rowHeightIn=0.32;
    rightMarginIn=1.1; % colorbar + its label
    gapAboveAxesIn=0.12;
    bracketTickIn=0.05;
    bracketLabelHeightIn=0.22;
    titleGapIn=0.08;
    titleHeightIn=0.3;
    topMarginIn=gapAboveAxesIn+bracketTickIn+bracketLabelHeightIn+titleGapIn+titleHeightIn;

    % Measure the widest metric label (at the font size the axes will
    % actually use) so the left margin fits it exactly -- same approach as
    % plotFactorViolinGrid's row labels, and for the same reason: a guessed
    % margin either clips long names or wastes width that could go to cells.
    measureFig=figure('Visible','off');
    measureAx=axes(measureFig,'Units','inches');
    maxLabelWidthIn=0;
    for m=1:nMetrics
        t=text(measureAx,0,0,MetricOrder(m),'FontSize',labelFontSize,'Interpreter','none','Units','inches');
        maxLabelWidthIn=max(maxLabelWidthIn,t.Extent(3));
    end
    close(measureFig);
    leftMarginIn=maxLabelWidthIn+0.25;

    % Measure the tallest rotated column label the same way, so the bottom
    % margin fits it exactly instead of a fixed guess that clips long
    % column labels (item 10, 2026-09-05) -- rotation matches
    % ax.XTickLabelRotation=30 set on the real axes below.
    measureFig2=figure('Visible','off');
    measureAx2=axes(measureFig2,'Units','inches');
    maxColLabelHeightIn=0;
    for c=1:nCols
        t=text(measureAx2,0,0,columnSpecs(c).Label,'FontSize',labelFontSize,'Interpreter','none','Units','inches','Rotation',30);
        maxColLabelHeightIn=max(maxColLabelHeightIn,t.Extent(4));
    end
    close(measureFig2);
    bottomLabelPadIn=0.25; % gap below the rotated column labels; also how far short of the figure bottom the colorbar stops (item 4, 2026-09-05)
    bottomMarginIn=maxColLabelHeightIn+bottomLabelPadIn;

    figWidthIn=leftMarginIn+nCols*colWidthIn+rightMarginIn;
    figHeightIn=topMarginIn+nMetrics*rowHeightIn+bottomMarginIn;
    axesBottomIn=bottomMarginIn;

    fig=figure('Visible','off');
    set(fig,'units','inches','position',[0 0 figWidthIn figHeightIn],'windowstyle','normal')

    ax=axes(fig,'Units','inches','Position',[leftMarginIn,axesBottomIn,nCols*colWidthIn,nMetrics*rowHeightIn]);
    hold(ax,'on')

    lim=max(abs(DisplayVal(:)),[],'omitnan');
    if isnan(lim) || lim==0
        lim=1;
    end
    img=imagesc(ax,[1,nCols],[1,nMetrics],DisplayVal,[-lim,lim]);
    img.AlphaData=~isnan(DisplayVal);
    ax.Color=[0.92,0.92,0.92]; % shows through where AlphaData=0 (no data)
    colormap(ax,neutralDivergingColormap())
    ax.YDir='reverse'; % row 1 (MetricOrder(1), e.g. "Rate") at the top

    for gx=0.5:1:(nCols+0.5)
        plot(ax,[gx,gx],[0.5,nMetrics+0.5],'Color',[1,1,1],'LineWidth',1,'HandleVisibility','off')
    end
    for gy=0.5:1:(nMetrics+0.5)
        plot(ax,[0.5,nCols+0.5],[gy,gy],'Color',[1,1,1],'LineWidth',1,'HandleVisibility','off')
    end

    for m=1:nMetrics
        for c=1:nCols
            if isnan(DisplayVal(m,c))
                continue
            end
            stars=sigStars(PValue(m,c));
            label=sprintf('%.2f',DisplayVal(m,c));
            if strlength(stars)>0
                % Back on one line (follow-up, 2026-09-05, reverting the
                % two-stacked-line version) but with a couple of literal
                % spaces between the number and the stars instead of
                % nothing, so there's a little breathing room without a
                % full line break.
                label=strcat(label,"  ",stars);
            end
            if abs(DisplayVal(m,c))>0.6*lim
                txtColor=[1,1,1]; % readable over the colormap's darker ends
            else
                txtColor=[0,0,0];
            end
            text(ax,c,m,label,'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',cellFontSize,'Color',txtColor)
        end
    end

    ax.XLim=[0.5,nCols+0.5];
    ax.YLim=[0.5,nMetrics+0.5];
    ax.XTick=1:nCols;
    ax.XTickLabel={columnSpecs.Label};
    ax.XTickLabelRotation=30;
    ax.YTick=1:nMetrics;
    ax.YTickLabel=cellstr(MetricOrder);
    ax.TickLength=[0,0];
    ax.FontSize=labelFontSize;
    box(ax,'on')

    cb=colorbar(ax);
    cb.Label.String="Effect Direction (t-weighted)";
    % colorbar() auto-shrinks its target axes to make room for itself --
    % rightMarginIn already reserves that space, so put both back exactly
    % where they were meant to be. Without this, the axes ends up narrower
    % than colWidthIn*nCols actually accounts for, so the bracket headers
    % below (measured against the ORIGINAL, unshrunk column width) run past
    % where the columns actually ended up on screen.
    ax.Units='inches';
    ax.Position=[leftMarginIn,axesBottomIn,nCols*colWidthIn,nMetrics*rowHeightIn];
    cb.Units='inches';
    cb.Position=[leftMarginIn+nCols*colWidthIn+0.15,bottomLabelPadIn,0.25,axesBottomIn+nMetrics*rowHeightIn-bottomLabelPadIn]; % item 6 (prior round) + item 4 (2026-09-05): reaches down to the x-axis labels but stops bottomLabelPadIn short of the true figure bottom, same gap the labels themselves have

    % Bracketed group headers + the title, both drawn as fixed-position
    % annotations (figure-normalized coordinates) rather than an axes
    % title -- MATLAB auto-places an axes title right above ax.Position,
    % which would land at the same height as these brackets and overlap
    % them, so everything above the axes is laid out manually instead.
    colEdgesIn=leftMarginIn+(0:nCols)*colWidthIn; % left edge of column c is colEdgesIn(c)
    bracketLineYIn=axesBottomIn+nMetrics*rowHeightIn+gapAboveAxesIn+bracketTickIn;
    c=1;
    while c<=nCols
        gl=columnSpecs(c).GroupLabel;
        if strlength(gl)==0
            c=c+1;
            continue
        end
        c2=c;
        while c2<nCols && columnSpecs(c2+1).GroupLabel==gl
            c2=c2+1;
        end
        xStartFrac=colEdgesIn(c)/figWidthIn;
        xEndFrac=colEdgesIn(c2+1)/figWidthIn;
        yLineFrac=bracketLineYIn/figHeightIn;
        tickFrac=bracketTickIn/figHeightIn;
        annotation(fig,'line',[xStartFrac,xEndFrac],[yLineFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'line',[xStartFrac,xStartFrac],[yLineFrac-tickFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'line',[xEndFrac,xEndFrac],[yLineFrac-tickFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'textbox',[xStartFrac,yLineFrac,xEndFrac-xStartFrac,bracketLabelHeightIn/figHeightIn], ...
            'String',gl,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
            'FontSize',labelFontSize,'Interpreter','none')
        c=c2+1;
    end
    titleYIn=bracketLineYIn+bracketLabelHeightIn+titleGapIn;
    annotation(fig,'textbox',[0,titleYIn/figHeightIn,1,titleHeightIn/figHeightIn], ...
        'String',titleStr,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize',labelFontSize+3,'Interpreter','none')

    saveas(fig,strcat(saveFile,".png"))
    close(fig)
end

function ax = drawAgeHeatmapPanel(fig, DisplayVal, PValue, MetricOrder, columnSpecs, leftMarginIn, bottomIn, colWidthIn, rowHeightIn, cellFontSize, labelFontSize, lim, showXAxis)
% drawAgeHeatmapPanel  One scope's imagesc grid + row labels + per-cell
% text for plotCombinedAgeEffectHeatmap -- the same cell/gridline/label
% drawing plotAgeEffectHeatmap does for its single axes, factored out so
% it can be called once per stacked panel. showXAxis controls whether this
% panel gets its own column tick labels (only the bottom-most panel should
% -- see plotCombinedAgeEffectHeatmap for why the column axis is shared).
    nMetrics=numel(MetricOrder);
    nCols=numel(columnSpecs);

    ax=axes(fig,'Units','inches','Position',[leftMarginIn,bottomIn,nCols*colWidthIn,nMetrics*rowHeightIn]);
    hold(ax,'on')

    img=imagesc(ax,[1,nCols],[1,nMetrics],DisplayVal,[-lim,lim]);
    img.AlphaData=~isnan(DisplayVal);
    ax.Color=[0.92,0.92,0.92]; % shows through where AlphaData=0 (no data)
    ax.YDir='reverse'; % row 1 (MetricOrder(1), e.g. "Rate") at the top

    for gx=0.5:1:(nCols+0.5)
        plot(ax,[gx,gx],[0.5,nMetrics+0.5],'Color',[1,1,1],'LineWidth',1,'HandleVisibility','off')
    end
    for gy=0.5:1:(nMetrics+0.5)
        plot(ax,[0.5,nCols+0.5],[gy,gy],'Color',[1,1,1],'LineWidth',1,'HandleVisibility','off')
    end

    for m=1:nMetrics
        for c=1:nCols
            if isnan(DisplayVal(m,c))
                continue
            end
            stars=sigStars(PValue(m,c));
            label=sprintf('%.2f',DisplayVal(m,c));
            if strlength(stars)>0
                % Back on one line (follow-up, 2026-09-05, reverting the
                % two-stacked-line version) but with a couple of literal
                % spaces between the number and the stars instead of
                % nothing, so there's a little breathing room without a
                % full line break.
                label=strcat(label,"  ",stars);
            end
            if abs(DisplayVal(m,c))>0.6*lim
                txtColor=[1,1,1]; % readable over the colormap's darker ends
            else
                txtColor=[0,0,0];
            end
            text(ax,c,m,label,'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',cellFontSize,'Color',txtColor)
        end
    end

    ax.XLim=[0.5,nCols+0.5];
    ax.YLim=[0.5,nMetrics+0.5];
    ax.YTick=1:nMetrics;
    ax.YTickLabel=cellstr(MetricOrder);
    if showXAxis
        ax.XTick=1:nCols;
        ax.XTickLabel={columnSpecs.Label};
        ax.XTickLabelRotation=30;
    else
        ax.XTick=[];
    end
    ax.TickLength=[0,0];
    ax.FontSize=labelFontSize;
    box(ax,'on')
end

function plotCombinedAgeEffectHeatmap(EffectsResults, MuSigmaResults, columnSpecs, MetricOrder, saveFile, titleStr)
% plotCombinedAgeEffectHeatmap  Functional (top) and Structural (bottom)
% age-fit heatmaps stacked in one figure, since these are much wider than
% they are tall. Columns (the age-related terms) are identical for both
% scopes, so the column axis is only drawn once: the bracket group headers
% above the Functional panel, and the rotated column tick labels below the
% Structural panel -- repeating either would just be clutter since the
% columns never differ between the two panels. Metric ROWS can differ
% (e.g. Structural has no "Same Subnet"), so each panel keeps its own row
% labels, with a gap row (drawn as plain background, same as any other
% missing cell) wherever one side has no data for a metric the other side
% does. Both panels share one color scale and one colorbar.
    [DisplayValF,PValueF,MetricOrderF]=gatherAgeHeatmapGrid(EffectsResults,MuSigmaResults,"Functional",columnSpecs,MetricOrder);
    [DisplayValS,PValueS,MetricOrderS]=gatherAgeHeatmapGrid(EffectsResults,MuSigmaResults,"Structural",columnSpecs,MetricOrder);

    keepMetrics=ismember(MetricOrder,MetricOrderF)|ismember(MetricOrder,MetricOrderS);
    UnionMetrics=MetricOrder(keepMetrics);
    nMetrics=numel(UnionMetrics);
    nCols=numel(columnSpecs);
    if nMetrics==0
        fprintf("Skipping combined age-fit heatmap '%s': no data for any metric on either side.\n",titleStr);
        return
    end

    % Re-index each side onto the shared union row order -- NaN (=
    % background) wherever that side has no data for a metric the other
    % side does, so the two panels stay aligned metric-for-metric.
    GridF=nan(nMetrics,nCols);
    PF=nan(nMetrics,nCols);
    GridS=nan(nMetrics,nCols);
    PS=nan(nMetrics,nCols);
    for m=1:nMetrics
        hitF=find(MetricOrderF==UnionMetrics(m),1);
        if ~isempty(hitF)
            GridF(m,:)=DisplayValF(hitF,:);
            PF(m,:)=PValueF(hitF,:);
        end
        hitS=find(MetricOrderS==UnionMetrics(m),1);
        if ~isempty(hitS)
            GridS(m,:)=DisplayValS(hitS,:);
            PS(m,:)=PValueS(hitS,:);
        end
    end

    % One companion text file per heatmap image (2026-09-05, user request)
    % -- see plotAgeEffectHeatmap's own txtFid for the full explanation.
    txtFid=fopen(strcat(saveFile,".txt"),'wt');
    fprintf(txtFid,'%s\n',titleStr);
    for m=1:nMetrics
        fprintf(txtFid,'%s\n',UnionMetrics(m));
        if ~all(isnan(GridF(m,:)))
            fprintf(txtFid,'  Functional:\n');
            writeHeatmapCellLines(txtFid,"    ",columnSpecs,GridF,PF,m);
        end
        if ~all(isnan(GridS(m,:)))
            fprintf(txtFid,'  Structural:\n');
            writeHeatmapCellLines(txtFid,"    ",columnSpecs,GridS,PS,m);
        end
    end
    fclose(txtFid);

    cellFontSize=8;
    labelFontSize=9;

    colWidthIn=0.95;
    rowHeightIn=0.32;
    rightMarginIn=1.1; % colorbar + its label
    gapAboveAxesIn=0.12;
    bracketTickIn=0.05;
    bracketLabelHeightIn=0.22;
    titleGapIn=0.08;
    titleHeightIn=0.3;
    panelGapIn=0.28; % gap between the Functional and Structural panels
    panelLabelHeightIn=0.18; % "Functional"/"Structural" tag above each panel

    % Measure the widest metric label (at the font size the axes will
    % actually use) so the left margin fits it exactly -- same approach as
    % plotAgeEffectHeatmap's own row labels.
    measureFig=figure('Visible','off');
    measureAx=axes(measureFig,'Units','inches');
    maxLabelWidthIn=0;
    for m=1:nMetrics
        t=text(measureAx,0,0,UnionMetrics(m),'FontSize',labelFontSize,'Interpreter','none','Units','inches');
        maxLabelWidthIn=max(maxLabelWidthIn,t.Extent(3));
    end
    close(measureFig);
    leftMarginIn=maxLabelWidthIn+0.25;

    % Measure the tallest rotated column label the same way, so the bottom
    % margin (below the Structural panel) fits it exactly instead of a
    % fixed guess that clips long column labels (item 10, 2026-09-05).
    measureFig2=figure('Visible','off');
    measureAx2=axes(measureFig2,'Units','inches');
    maxColLabelHeightIn=0;
    for c=1:nCols
        t=text(measureAx2,0,0,columnSpecs(c).Label,'FontSize',labelFontSize,'Interpreter','none','Units','inches','Rotation',30);
        maxColLabelHeightIn=max(maxColLabelHeightIn,t.Extent(4));
    end
    close(measureFig2);
    bottomLabelPadIn=0.25; % gap below the rotated column labels; also how far short of the figure bottom the colorbar stops (item 4, 2026-09-05)
    bottomMarginIn=maxColLabelHeightIn+bottomLabelPadIn; % rotated column labels, once, below the Structural panel

    figWidthIn=leftMarginIn+nCols*colWidthIn+rightMarginIn;

    % Bottom-up layout: Structural axes, its label strip, the gap, the
    % Functional axes, its label strip, the bracket headers, then the title.
    sPanelBottomIn=bottomMarginIn;
    sLabelYIn=sPanelBottomIn+nMetrics*rowHeightIn;
    fPanelBottomIn=sLabelYIn+panelLabelHeightIn+panelGapIn;
    fLabelYIn=fPanelBottomIn+nMetrics*rowHeightIn;
    bracketLineYIn=fLabelYIn+panelLabelHeightIn+gapAboveAxesIn+bracketTickIn;
    titleYIn=bracketLineYIn+bracketLabelHeightIn+titleGapIn;
    figHeightIn=titleYIn+titleHeightIn;

    fig=figure('Visible','off');
    set(fig,'units','inches','position',[0 0 figWidthIn figHeightIn],'windowstyle','normal')

    lim=max(abs([GridF(:);GridS(:)]),[],'omitnan');
    if isnan(lim) || lim==0
        lim=1;
    end

    axS=drawAgeHeatmapPanel(fig,GridS,PS,UnionMetrics,columnSpecs,leftMarginIn,sPanelBottomIn,colWidthIn,rowHeightIn,cellFontSize,labelFontSize,lim,true);
    axF=drawAgeHeatmapPanel(fig,GridF,PF,UnionMetrics,columnSpecs,leftMarginIn,fPanelBottomIn,colWidthIn,rowHeightIn,cellFontSize,labelFontSize,lim,false);
    colormap(fig,neutralDivergingColormap())

    % colorbar() auto-shrinks its target axes to make room for itself --
    % rightMarginIn already reserves that space, so put both axes back
    % exactly where they were meant to be (same fix as plotAgeEffectHeatmap).
    cb=colorbar(axS);
    cb.Label.String="Effect Direction (t-weighted)";
    axS.Units='inches';
    axS.Position=[leftMarginIn,sPanelBottomIn,nCols*colWidthIn,nMetrics*rowHeightIn];
    axF.Units='inches';
    axF.Position=[leftMarginIn,fPanelBottomIn,nCols*colWidthIn,nMetrics*rowHeightIn];
    cb.Units='inches';
    cb.Position=[leftMarginIn+nCols*colWidthIn+0.15,bottomLabelPadIn,0.25,fLabelYIn-bottomLabelPadIn]; % item 6 (prior round) + item 4 (2026-09-05): reaches down to the x-axis labels but stops bottomLabelPadIn short of the true figure bottom, same gap the labels themselves have

    annotation(fig,'textbox',[leftMarginIn/figWidthIn,sLabelYIn/figHeightIn,nCols*colWidthIn/figWidthIn,panelLabelHeightIn/figHeightIn], ...
        'String',"Structural",'HorizontalAlignment','left','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize',labelFontSize+1,'FontWeight','bold','Interpreter','none')
    annotation(fig,'textbox',[leftMarginIn/figWidthIn,fLabelYIn/figHeightIn,nCols*colWidthIn/figWidthIn,panelLabelHeightIn/figHeightIn], ...
        'String',"Functional",'HorizontalAlignment','left','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize',labelFontSize+1,'FontWeight','bold','Interpreter','none')

    % Bracketed group headers + the title, once, above the Functional panel.
    colEdgesIn=leftMarginIn+(0:nCols)*colWidthIn; % left edge of column c is colEdgesIn(c)
    c=1;
    while c<=nCols
        gl=columnSpecs(c).GroupLabel;
        if strlength(gl)==0
            c=c+1;
            continue
        end
        c2=c;
        while c2<nCols && columnSpecs(c2+1).GroupLabel==gl
            c2=c2+1;
        end
        xStartFrac=colEdgesIn(c)/figWidthIn;
        xEndFrac=colEdgesIn(c2+1)/figWidthIn;
        yLineFrac=bracketLineYIn/figHeightIn;
        tickFrac=bracketTickIn/figHeightIn;
        annotation(fig,'line',[xStartFrac,xEndFrac],[yLineFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'line',[xStartFrac,xStartFrac],[yLineFrac-tickFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'line',[xEndFrac,xEndFrac],[yLineFrac-tickFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'textbox',[xStartFrac,yLineFrac,xEndFrac-xStartFrac,bracketLabelHeightIn/figHeightIn], ...
            'String',gl,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
            'FontSize',labelFontSize,'Interpreter','none')
        c=c2+1;
    end
    annotation(fig,'textbox',[0,titleYIn/figHeightIn,1,titleHeightIn/figHeightIn], ...
        'String',titleStr,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize',labelFontSize+3,'Interpreter','none')

    saveas(fig,strcat(saveFile,".png"))
    close(fig)
end

function plotSubnetAggregateAgeHeatmap(EffectsResults, MuSigmaResults, scopeName, subnetOrder, columnSpecs, MetricOrder, saveFile, titleStr)
% plotSubnetAggregateAgeHeatmap  Age-fit heatmaps for every subnetwork in
% ONE scope (Functional or Structural -- call once per scope, kept
% separate, same reasoning as plotSubnetAggregateViolinGrid), stacked
% vertically: one panel per SUBNETWORK (item 16, 2026-09-05 -- reversed
% from the original one-panel-per-metric layout per user direction),
% rows=metrics that apply at the subnetwork level (almost always a small
% subset of the whole-brain metric list, since whole-brain-only metrics
% like Assortativity or the Same-Subnet family were never computed per
% subnetwork), columns=the same age-related terms plotAgeEffectHeatmap
% uses. Reuses drawAgeHeatmapPanel as-is with metric names in its
% row-label slot (same as plotAgeEffectHeatmap itself), and the same
% bottom-up stacking layout plotCombinedAgeEffectHeatmap uses (column
% axis/brackets drawn once, one shared color scale/colorbar for every
% panel), generalized from that function's 2 fixed panels
% (Functional/Structural) to however many subnetworks apply at this scope.
    subnetworkList=strcat(scopeName," ",subnetOrder);

    % Gather each subnetwork's own Metric x Term grid.
    nSubnetsIn=numel(subnetOrder);
    perSubnetGrid=cell(1,nSubnetsIn);
    perSubnetPVal=cell(1,nSubnetsIn);
    perSubnetMetrics=cell(1,nSubnetsIn);
    hasSubnet=false(1,nSubnetsIn);
    for s=1:nSubnetsIn
        [DisplayVal,PValue,thisMetricOrder]=gatherAgeHeatmapGrid(EffectsResults,MuSigmaResults,subnetworkList(s),columnSpecs,MetricOrder);
        if ~isempty(thisMetricOrder)
            hasSubnet(s)=true;
            perSubnetGrid{s}=DisplayVal;
            perSubnetPVal{s}=PValue;
            perSubnetMetrics{s}=thisMetricOrder;
        end
    end
    subnetOrder=subnetOrder(hasSubnet);
    perSubnetGrid=perSubnetGrid(hasSubnet);
    perSubnetPVal=perSubnetPVal(hasSubnet);
    perSubnetMetrics=perSubnetMetrics(hasSubnet);
    nSubnets=numel(subnetOrder);
    if nSubnets==0
        fprintf("Skipping subnetwork-aggregate age-fit heatmap '%s': no subnetwork data.\n",titleStr);
        return
    end

    % Union of metrics across every subnetwork, canonical order -- now the
    % shared ROW order used inside every subnetwork panel (item 16,
    % 2026-09-05), not the panel-defining axis itself. Panels are now
    % subnetworks (nPanels=nSubnets, below).
    % FIXED 2026-09-05: was false(1,numel(MetricOrder)) -- a hardcoded row
    % vector, but MetricOrder (MetricNames, defined via ["Rate";longnames2])
    % is actually a column vector, so ismember(MetricOrder,...) returned a
    % column and the | below silently broadcast 1xN against Nx1 into an
    % NxN matrix, corrupting keepMetric and causing an out-of-bounds
    % logical index. false(size(MetricOrder)) always matches MetricOrder's
    % own shape.
    keepMetric=false(size(MetricOrder));
    for s=1:nSubnets
        keepMetric=keepMetric|ismember(MetricOrder,perSubnetMetrics{s});
    end
    UnionMetrics=MetricOrder(keepMetric);
    nRows=numel(UnionMetrics);
    nCols=numel(columnSpecs);
    nPanels=nSubnets;
    if nRows==0
        fprintf("Skipping subnetwork-aggregate age-fit heatmap '%s': no metric data for any subnetwork.\n",titleStr);
        return
    end

    % Re-index each SUBNETWORK's grid onto the shared union metric row
    % order -- one Metric x Term matrix per subnetwork PANEL now (item 16),
    % NaN (=background) wherever a subnetwork has no data for a given
    % metric. Same reindex-onto-a-union-row-order pattern
    % plotCombinedAgeEffectHeatmap uses for its GridF/GridS, applied here
    % per-subnetwork instead of per-scope.
    PanelGrids=cell(1,nPanels);
    PanelPVals=cell(1,nPanels);
    for s=1:nPanels
        G=nan(nRows,nCols);
        P=nan(nRows,nCols);
        for r=1:nRows
            hit=find(perSubnetMetrics{s}==UnionMetrics(r),1);
            if ~isempty(hit)
                G(r,:)=perSubnetGrid{s}(hit,:);
                P(r,:)=perSubnetPVal{s}(hit,:);
            end
        end
        PanelGrids{s}=G;
        PanelPVals{s}=P;
    end

    % One companion text file per heatmap image (2026-09-05, user request)
    % -- see plotAgeEffectHeatmap's own txtFid for the full explanation.
    % Panels here are subnetworks (item 16), rows within each panel are
    % metrics, so the nesting is subnetwork -> metric -> term.
    txtFid=fopen(strcat(saveFile,".txt"),'wt');
    fprintf(txtFid,'%s\n',titleStr);
    for s=1:nPanels
        fprintf(txtFid,'%s\n',subnetOrder(s));
        for r=1:nRows
            if all(isnan(PanelGrids{s}(r,:)))
                continue
            end
            fprintf(txtFid,'  %s\n',UnionMetrics(r));
            writeHeatmapCellLines(txtFid,"    ",columnSpecs,PanelGrids{s},PanelPVals{s},r);
        end
    end
    fclose(txtFid);

    cellFontSize=8;
    labelFontSize=9;

    colWidthIn=0.95;
    rowHeightIn=0.32;
    rightMarginIn=1.1; % colorbar + its label
    gapAboveAxesIn=0.12;
    bracketTickIn=0.05;
    bracketLabelHeightIn=0.22;
    titleGapIn=0.08;
    titleHeightIn=0.3;
    panelGapIn=0.28; % gap between subnetwork panels
    panelLabelHeightIn=0.18; % subnetwork-name tag above each panel

    % Row labels are now metric names (item 16), so measure UnionMetrics
    % instead of subnetOrder for the left margin.
    measureFig=figure('Visible','off');
    measureAx=axes(measureFig,'Units','inches');
    maxLabelWidthIn=0;
    for r=1:nRows
        t=text(measureAx,0,0,UnionMetrics(r),'FontSize',labelFontSize,'Interpreter','none','Units','inches');
        maxLabelWidthIn=max(maxLabelWidthIn,t.Extent(3));
    end
    close(measureFig);
    leftMarginIn=maxLabelWidthIn+0.25;

    % Measure the tallest rotated column label the same way, so the bottom
    % margin fits it exactly instead of a fixed guess that clips long
    % column labels (item 10, 2026-09-05).
    measureFig2=figure('Visible','off');
    measureAx2=axes(measureFig2,'Units','inches');
    maxColLabelHeightIn=0;
    for c=1:nCols
        t=text(measureAx2,0,0,columnSpecs(c).Label,'FontSize',labelFontSize,'Interpreter','none','Units','inches','Rotation',30);
        maxColLabelHeightIn=max(maxColLabelHeightIn,t.Extent(4));
    end
    close(measureFig2);
    bottomLabelPadIn=0.25; % gap below the rotated column labels; also how far short of the figure bottom the colorbar stops (item 4, 2026-09-05)
    bottomMarginIn=maxColLabelHeightIn+bottomLabelPadIn; % rotated column labels, once, below the bottom-most panel

    figWidthIn=leftMarginIn+nCols*colWidthIn+rightMarginIn;

    % Bottom-up layout: panel nPanels (bottom-most) up through panel 1
    % (top-most) -- generalizes plotCombinedAgeEffectHeatmap's fixed
    % Structural-then-Functional stack to however many panels there are.
    panelBottomIn=zeros(1,nPanels);
    y=bottomMarginIn;
    for p=nPanels:-1:1
        panelBottomIn(p)=y;
        y=y+nRows*rowHeightIn+panelLabelHeightIn;
        if p>1
            y=y+panelGapIn;
        end
    end
    bracketLineYIn=y+gapAboveAxesIn+bracketTickIn;
    titleYIn=bracketLineYIn+bracketLabelHeightIn+titleGapIn;
    figHeightIn=titleYIn+titleHeightIn;

    fig=figure('Visible','off');
    set(fig,'units','inches','position',[0 0 figWidthIn figHeightIn],'windowstyle','normal')

    lim=0;
    for p=1:nPanels
        lim=max(lim,max(abs(PanelGrids{p}(:)),[],'omitnan'));
    end
    if isnan(lim) || lim==0
        lim=1;
    end

    ax=gobjects(1,nPanels);
    for p=1:nPanels
        ax(p)=drawAgeHeatmapPanel(fig,PanelGrids{p},PanelPVals{p},UnionMetrics,columnSpecs,leftMarginIn,panelBottomIn(p),colWidthIn,rowHeightIn,cellFontSize,labelFontSize,lim,p==nPanels);
    end
    colormap(fig,neutralDivergingColormap())

    % colorbar() auto-shrinks its target axes to make room for itself --
    % rightMarginIn already reserves that space, so put every axes back
    % exactly where it was meant to be (same fix as plotAgeEffectHeatmap).
    cb=colorbar(ax(nPanels));
    cb.Label.String="Effect Direction (t-weighted)";
    for p=1:nPanels
        ax(p).Units='inches';
        ax(p).Position=[leftMarginIn,panelBottomIn(p),nCols*colWidthIn,nRows*rowHeightIn];
    end
    cb.Units='inches';
    cb.Position=[leftMarginIn+nCols*colWidthIn+0.15,bottomLabelPadIn,0.25,panelBottomIn(1)+nRows*rowHeightIn-bottomLabelPadIn]; % item 6 (prior round) + item 4 (2026-09-05): reaches down to the x-axis labels but stops bottomLabelPadIn short of the true figure bottom, same gap the labels themselves have

    for p=1:nPanels
        labelYIn=panelBottomIn(p)+nRows*rowHeightIn;
        annotation(fig,'textbox',[leftMarginIn/figWidthIn,labelYIn/figHeightIn,nCols*colWidthIn/figWidthIn,panelLabelHeightIn/figHeightIn], ...
            'String',subnetOrder(p),'HorizontalAlignment','left','VerticalAlignment','bottom','EdgeColor','none', ...
            'FontSize',labelFontSize+1,'FontWeight','bold','Interpreter','none')
    end

    % Bracketed group headers + the title, once, above the top-most panel.
    colEdgesIn=leftMarginIn+(0:nCols)*colWidthIn; % left edge of column c is colEdgesIn(c)
    c=1;
    while c<=nCols
        gl=columnSpecs(c).GroupLabel;
        if strlength(gl)==0
            c=c+1;
            continue
        end
        c2=c;
        while c2<nCols && columnSpecs(c2+1).GroupLabel==gl
            c2=c2+1;
        end
        xStartFrac=colEdgesIn(c)/figWidthIn;
        xEndFrac=colEdgesIn(c2+1)/figWidthIn;
        yLineFrac=bracketLineYIn/figHeightIn;
        tickFrac=bracketTickIn/figHeightIn;
        annotation(fig,'line',[xStartFrac,xEndFrac],[yLineFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'line',[xStartFrac,xStartFrac],[yLineFrac-tickFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'line',[xEndFrac,xEndFrac],[yLineFrac-tickFrac,yLineFrac],'Color','k','LineWidth',1)
        annotation(fig,'textbox',[xStartFrac,yLineFrac,xEndFrac-xStartFrac,bracketLabelHeightIn/figHeightIn], ...
            'String',gl,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
            'FontSize',labelFontSize,'Interpreter','none')
        c=c2+1;
    end
    annotation(fig,'textbox',[0,titleYIn/figHeightIn,1,titleHeightIn/figHeightIn], ...
        'String',titleStr,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize',labelFontSize+3,'Interpreter','none')

    saveas(fig,strcat(saveFile,".png"))
    close(fig)
end

function [laneOfComparison, nLanesActual, laneWidths, compactLaneWidthIn, wideLaneWidthIn, bracketColWidthIn] = planComparisonLanes(templateComparisons)
% planComparisonLanes  Pack a comparisons struct array (leftIdx,rightIdx,
% pValue) into lanes for plotFactorViolinGrid/plotCombinedViolinGrid: two
% brackets whose group spans never touch or overlap can share one lane "at
% the same level," since they can never visually collide (e.g. Male CN-DS
% and Female CN-DS -- disjoint groups entirely); anything that would
% collide gets its own lane, in order. This is just greedy interval
% packing. bracketColWidthIn is the worst-case total width needed (every
% lane active at once, including the compact-before-wide spacing bump).
    nComparisonsTemplate=numel(templateComparisons);
    laneOfComparison=zeros(1,nComparisonsTemplate);
    laneRangesLo={};
    laneRangesHi={};
    nLanesActual=0;
    for k=1:nComparisonsTemplate
        lo=min(templateComparisons(k).leftIdx,templateComparisons(k).rightIdx);
        hi=max(templateComparisons(k).leftIdx,templateComparisons(k).rightIdx);
        placedLane=0;
        for L=1:nLanesActual
            conflict=false;
            for j=1:numel(laneRangesLo{L})
                if ~(hi<laneRangesLo{L}(j) || laneRangesHi{L}(j)<lo)
                    conflict=true;
                    break
                end
            end
            if ~conflict
                placedLane=L;
                break
            end
        end
        if placedLane==0
            nLanesActual=nLanesActual+1;
            placedLane=nLanesActual;
            laneRangesLo{placedLane}=[]; %#ok<AGROW>
            laneRangesHi{placedLane}=[]; %#ok<AGROW>
        end
        laneOfComparison(k)=placedLane;
        laneRangesLo{placedLane}(end+1)=lo; %#ok<AGROW>
        laneRangesHi{placedLane}(end+1)=hi; %#ok<AGROW>
    end

    compactLaneWidthIn=0.22;
    wideLaneWidthIn=0.48;
    laneWidths=zeros(1,nLanesActual);
    for L=1:nLanesActual
        if all((laneRangesHi{L}-laneRangesLo{L})==1)
            laneWidths(L)=compactLaneWidthIn;
        else
            laneWidths(L)=wideLaneWidthIn;
        end
    end
    fullWidths=laneWidths;
    for L=1:nLanesActual-1
        if fullWidths(L)==compactLaneWidthIn && fullWidths(L+1)==wideLaneWidthIn
            fullWidths(L)=fullWidths(L)+0.05;
        end
    end
    bracketColWidthIn=sum(fullWidths);
end

function [activeLanes, activeWidths, totalWidthIn] = packActiveLanes(comparisons, laneOfComparison, nLanesActual, laneWidths, compactLaneWidthIn, wideLaneWidthIn)
% packActiveLanes  Given one row's comparisons and the grid-wide lane
% assignment (from planComparisonLanes), returns which lanes actually have
% a significant result THIS row (activeLanes), their packed widths
% (activeWidths, with the same compact-before-wide spacing bump
% planComparisonLanes itself applies to the worst case), and the total
% width needed. Shared (2026-09-05, item 9) by drawViolinGridPanel, which
% uses activeLanes/activeWidths to position this row's own brackets, and
% by each grid function's column-width sizing pass, which calls this once
% per row and takes the max totalWidthIn across the column -- replacing
% planComparisonLanes' own bracketColWidthIn (which assumed EVERY
% comparison is significant, the worst case that's structurally possible)
% with the worst case that actually occurs across the rows that exist,
% freeing up plotWidthIn on every row that needs fewer active lanes than
% that.
    laneHasStars=false(1,nLanesActual);
    for k=1:numel(comparisons)
        if strlength(sigStars(comparisons(k).pValue))>0
            laneHasStars(laneOfComparison(k))=true;
        end
    end
    activeLanes=find(laneHasStars);
    activeWidths=laneWidths(activeLanes);
    for L=1:numel(activeWidths)-1
        if activeWidths(L)==compactLaneWidthIn && activeWidths(L+1)==wideLaneWidthIn
            activeWidths(L)=activeWidths(L)+0.05;
        end
    end
    totalWidthIn=sum(activeWidths);
end

function boxHandle = drawViolinGridPanel(fig, ax, row, groupLevels, democon, demoorder, democheat, isRateRow, nGroups, rowWidthIn, bracketGapIn, bracketTickLenIn, starFontSize, laneOfComparison, nLanesActual, laneWidths, compactLaneWidthIn, wideLaneWidthIn, nTicks, tickFontSize, boxStartIn, boxWidthIn)
% drawViolinGridPanel  Draws one metric row's rotated violin set plus its
% significance brackets into an already-positioned axes -- the shared core
% of plotFactorViolinGrid's per-row rendering, factored out so
% plotCombinedViolinGrid can draw two of these (Functional and Structural)
% side by side in one row figure. Caller creates/positions ax (sets
% ax.Units='normalized', ax.Position, ax.FontSize) and handles anything
% outside the plot area itself (row label, x-axis label).
%
% Returns boxHandle, the overallEffects summary box's own annotation
% object (or empty if this row has none) -- it lives on fig rather than ax,
% so the caller can't discover it via ax.Children the way it can for the
% violin/errorbar/mean-line data. The caller needs this handle to toggle
% the box's visibility independently when building the grid's SigBoxOnly
% bonus composite (2026-09-05, user request): hidden on the original save,
% shown (with the violin/errorbar data hidden instead) on the bonus save --
% brackets stay visible in both, exactly like plotGroupViolin's own
% box/bracket treatment.
    boxHandle=gobjects(0);
    axes(ax); % violin() always draws into gca and takes no axes argument
    hold(ax,'on')
    xMinAll=Inf;
    xMaxAll=-Inf;
    for k=1:nGroups
        lvl=groupLevels(k);
        [h,vLine,~,~]=violin({row.values{k}},'mc','','medc','', ...
            'facecolor',demoorder(democheat==lvl,:),'facealpha',0.75,'edgecolor',democon(democheat==lvl,:));
        vLine.Visible='off';
        if isRateRow
            h=clipViolinAtZero(h);
        end
        h.XData=h.XData+(k-1);
        % Rotate this violin 90 degrees: swap the position axis and the
        % value axis, so the value reads left-to-right instead of
        % bottom-to-top.
        oldX=h.XData;
        oldY=h.YData;
        h.XData=oldY;
        h.YData=oldX;
        xMinAll=min([xMinAll,h.XData(:)']);
        xMaxAll=max([xMaxAll,h.XData(:)']);
        plot(ax,[row.mu(k),row.mu(k)], ...
            [interp1(h.XData(h.YData<k),h.YData(h.YData<k),row.mu(k)),interp1(h.XData(h.YData>k),h.YData(h.YData>k),row.mu(k))], ...
            'Color',democon(democheat==lvl,:),'LineWidth',1.5)
        errorbar(ax,row.mu(k),k,row.sigma(k),'horizontal','Color',democon(democheat==lvl,:))
    end
    if ~isfinite(xMinAll) || ~isfinite(xMaxAll) || xMinAll==xMaxAll
        xMinAll=-1;
        xMaxAll=1;
    end
    if isRateRow
        xSpan=xMaxAll-xMinAll;
        xlim(ax,[xMinAll-0.05*xSpan,xMaxAll+0.05*xSpan])
    else
        % Zero-center every non-Rate row on the same half-width scale
        % (rather than each row's own min/max), so the zero line lines up
        % at the same horizontal position across all rows -- a row skewed
        % mostly positive or negative then visibly leans to one side
        % instead of auto-filling the frame.
        halfWidth=max(abs(xMinAll),abs(xMaxAll))*1.05;
        xlim(ax,[-halfWidth,halfWidth])
    end
    % Same tick COUNT on every row (not the same numbers -- each row keeps
    % its own scale, computed independently, so this doesn't force
    % alignment across rows) so the axis doesn't look cluttered on rows
    % with a wide auto-chosen range. nTicks/tickFontSize are
    % caller-supplied (fewer, smaller for the narrow subnet-aggregate
    % panels) so numeric tick labels don't overlap; XTickLabelRotation is
    % forced to 0 (item 11, 2026-09-05) since MATLAB was otherwise
    % auto-rotating them diagonal when too many were crammed into a
    % narrow panel, and diagonal labels ran past the panel and got
    % clipped. niceTicksThroughZero (item 7, 2026-09-05) also gives the
    % Rate row a labeled zero tick whenever zero falls in its own range,
    % same as every other (already zero-centered) row.
    % These rows are stacked as separate small PNGs with a tight fixed
    % margin, and a middle row's auto power-of-ten exponent label had
    % nowhere reserved to render (only the bottom-most row's
    % xlabelExtraIn gives any extra room) -- it was silently
    % invisible/clipped (item 3, prior round). Explicit tick label strings
    % (item 2/3 follow-up, 2026-09-05, replacing ax.XAxis.Exponent=0 +
    % xtickformat) bypass MATLAB's exponent/ruler formatting entirely --
    % that combination turned out to be unreliable (rounding small-scale
    % values to "0.00" in some plots) -- so there's never an exponent
    % label to reserve room for, and our own sprintf('%.2g',...) is
    % exactly what gets shown, no second-guessing.
    tickValsX=niceTicksThroughZero(ax.XLim(1),ax.XLim(2),nTicks);
    xticks(ax,tickValsX)
    xticklabels(ax,arrayfun(@(v) sprintf('%.2g',v), tickValsX, 'UniformOutput', false))
    ax.XAxis.FontSize=tickFontSize;
    ax.XTickLabelRotation=0;
    ylim(ax,[0.3,nGroups+0.7])
    yticks(ax,[])
    xline(ax,0,'--','Color','k')
    box(ax,'off')

    % The bracket column normally starts at the fixed right edge of the
    % plot area, but when this row's data sits far from that edge (e.g. a
    % strongly negative-skewed metric, whose real values barely reach past
    % zero on the symmetric zero-centered axis), that leaves a wide empty
    % gap between the violins and their brackets, making it hard to tell
    % which belongs to which. Starting the column just past this row's own
    % actual rightmost data point instead keeps the brackets visually
    % anchored to their violins. This only READS the already-finalized
    % xlim/Position (both already set above) -- it never feeds back into
    % them, so it can't cause the zero-centering distortion a fully dynamic
    % xlim would.
    xl=ax.XLim;
    dataMaxFrac=ax.Position(1)+(xMaxAll-xl(1))/(xl(2)-xl(1))*ax.Position(3);
    rowBracketStartIn=dataMaxFrac*rowWidthIn+bracketGapIn;

    % Which lanes actually have something to draw this row -- a
    % non-significant test leaves its lane empty, and rather than reserving
    % its width anyway (a dead visual gap where a bracket would have been),
    % the remaining significant lanes are packed together starting right at
    % the plot edge, in lane order. packActiveLanes (item 9, 2026-09-05)
    % also drives the column-wide width sizing pass in each caller.
    [activeLanes,activeWidths]=packActiveLanes(row.comparisons,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn);
    activeOffsetsIn=cumsum([0,activeWidths(1:end-1)]);
    laneOffsetLookup=nan(1,nLanesActual);
    laneWidthLookup=nan(1,nLanesActual);
    laneOffsetLookup(activeLanes)=activeOffsetsIn;
    laneWidthLookup(activeLanes)=activeWidths;

    % Significance brackets, one per row.comparisons entry, drawn in the
    % lane laneOfComparison assigns it to (figure-normalized coordinates
    % via annotation, not data coordinates). Drawn as an actual bracket (a
    % vertical spine plus two short end-ticks pointing back at the plot)
    % rather than a plain line, spanning only its own two groups'
    % y-positions rather than assuming the full row height, so it
    % generalizes past the 2-group case.
    axBottomFrac=ax.Position(2); % already normalized -- caller sets ax.Units='normalized'
    axHeightFrac=ax.Position(4);
    tickLenFrac=bracketTickLenIn/rowWidthIn;
    for k=1:numel(row.comparisons)
        comp=row.comparisons(k);
        stars=sigStars(comp.pValue);
        if strlength(stars)==0
            continue
        end
        laneIdx=laneOfComparison(k);
        laneXStartIn=rowBracketStartIn+laneOffsetLookup(laneIdx);
        thisLaneWidthIn=laneWidthLookup(laneIdx);
        lineXFrac=laneXStartIn/rowWidthIn;
        textXFrac=(laneXStartIn+0.015)/rowWidthIn;
        textWidthFrac=(thisLaneWidthIn-0.015)/rowWidthIn;
        yLoFrac=axBottomFrac+(min(comp.leftIdx,comp.rightIdx)-0.3)/(nGroups+0.4)*axHeightFrac;
        yHiFrac=axBottomFrac+(max(comp.leftIdx,comp.rightIdx)-0.3)/(nGroups+0.4)*axHeightFrac;
        annotation(fig,'line',[lineXFrac,lineXFrac],[yLoFrac,yHiFrac],'Color','k','LineWidth',1.2);
        annotation(fig,'line',[lineXFrac-tickLenFrac,lineXFrac],[yLoFrac,yLoFrac],'Color','k','LineWidth',1.2);
        annotation(fig,'line',[lineXFrac-tickLenFrac,lineXFrac],[yHiFrac,yHiFrac],'Color','k','LineWidth',1.2);
        annotation(fig,'textbox',[textXFrac,(yLoFrac+yHiFrac)/2-0.5*axHeightFrac,textWidthFrac,axHeightFrac], ...
            'String',stars,'FontSize',starFontSize,'EdgeColor','none', ...
            'VerticalAlignment','middle','HorizontalAlignment','left');
    end

    % Grid/aggregate significance box: only the ANOVA rows carry
    % row.overallEffects (T-Test GridData has no such field), so this is a
    % no-op for those. Condensed single-line format per the user's explicit
    % choice, e.g. "Impair.* Sex n.s. Inter.**" -- reuses boxSigLabel so
    % stars/n.s./tex font markup match the individual-plot box exactly.
    % boxStartIn/boxWidthIn are caller-supplied FIXED positions (not
    % derived from this row's own dynamic rowBracketStartIn), so every
    % row's box lines up in the same column regardless of how far right
    % that row's own brackets happen to sit -- sized once by the caller
    % from the worst-case bracket reservation, so it's guaranteed clear of
    % every row's brackets. Drawn here with no border/background, same as
    % ever -- the caller now hides this handle on the original save and
    % shows only it (with everything else hidden) on the SigBoxOnly bonus
    % save, mirroring plotGroupViolin's own box treatment (2026-09-05, user
    % request).
    if isfield(row,'overallEffects') && ~isempty(row.overallEffects)
        boxStr="";
        for m=1:numel(row.overallEffects)
            if m>1
                boxStr=strcat(boxStr,"   ");
            end
            % Non-breaking space (char 160) between a label and its own
            % marker, not an ordinary space (2026-09-05, user report): this
            % box auto-wraps within its fixed width, and MATLAB's line
            % breaking treats any ordinary space as a valid wrap point --
            % which was splitting e.g. "Sex" onto one line and "n.s." onto
            % the next. A non-breaking space still renders as a normal gap
            % but is never chosen as a wrap point, so a label and its
            % marker now always wrap together as one unit; the "   "
            % separator between different effects above is left as ordinary
            % (breakable) spaces on purpose, since wrapping BETWEEN effects
            % is fine.
            boxStr=strcat(boxStr,row.overallEffects(m).label,string(char(160)),boxSigLabel(row.overallEffects(m).pValue,tickFontSize));
        end
        boxXFrac=boxStartIn/rowWidthIn;
        boxWidthFrac=boxWidthIn/rowWidthIn;
        boxHandle=annotation(fig,'textbox',[boxXFrac,axBottomFrac,boxWidthFrac,axHeightFrac], ...
            'String',boxStr,'FontSize',tickFontSize,'EdgeColor','none', ...
            'VerticalAlignment','middle','HorizontalAlignment','left','Interpreter','tex');
    end
end

function placeMissingPanelXLabel(fig, panelLeftIn, figWidthIn, plotWidthFrac, bottomMarginIn, thisRowHeightIn, metricName, baseFontSize)
% placeMissingPanelXLabel  plotCombinedViolinGrid's bottom row still needs
% its axis caption ("Rate (rho)"/"Factor Weight") even when a metric has
% no data on one side (on purpose) and so has no axes there for a normal
% xlabel to attach to (item 5, 2026-09-05). Uses annotation('textbox',...)
% in figure-normalized coordinates (not a plain text(), which would
% auto-create a stray default axes since none exists yet on this side) at
% the same span/height a real xlabel would have occupied.
    if metricName=="Rate"
        xLabelStr="Rate (\rho)";
    else
        xLabelStr="Factor Weight";
    end
    annotation(fig,'textbox',[panelLeftIn/figWidthIn,0,plotWidthFrac,bottomMarginIn/thisRowHeightIn], ...
        'String',xLabelStr,'HorizontalAlignment','center','VerticalAlignment','bottom','EdgeColor','none', ...
        'FontSize',baseFontSize,'Interpreter','tex')
end

function stackRowImages(rowFiles, saveFile)
% stackRowImages  Reads each row image, pads to a common width, stacks them
% top-to-bottom into one composite PNG at saveFile, then deletes the row
% files (caller still owns removing the now-empty temp directory). Shared
% by plotFactorViolinGrid and plotCombinedViolinGrid.
    rowImages=cell(numel(rowFiles),1);
    maxW=0;
    for r=1:numel(rowFiles)
        rowImages{r}=imread(rowFiles(r));
        maxW=max(maxW,size(rowImages{r},2));
    end
    for r=1:numel(rowFiles)
        im=rowImages{r};
        if size(im,2)<maxW
            im=[im,255*ones(size(im,1),maxW-size(im,2),size(im,3),class(im))]; %#ok<AGROW>
            rowImages{r}=im;
        end
    end
    imwrite(cat(1,rowImages{:}),strcat(saveFile,".png"));
    for r=1:numel(rowFiles)
        delete(rowFiles(r));
    end
end

function v = fieldOrEmpty(s, fieldName)
% fieldOrEmpty  s.(fieldName) if that field exists, else [] -- avoids an
% isfield check at every grid/aggregate call site that may or may not carry
% overallEffects (T-Test GridData never does; ANOVA GridData always does).
    if isfield(s,fieldName)
        v=s.(fieldName);
    else
        v=[];
    end
end

function writeComparisonEffectLines(txtFid, indent, comparisons, overallEffects)
% writeComparisonEffectLines  Writes one "<indent><label>: p= <value>" line
% per comparison and per overallEffects entry (either array can be empty)
% -- the shared line format for every violin grid/aggregate companion .txt
% file (2026-09-05, user request), matching plotGroupViolin's own per-plot
% text-file convention so the numbers behind a grid's brackets/box are
% available without eyeballing stars.
    for k=1:numel(comparisons)
        fprintf(txtFid,'%s',strcat(indent,comparisons(k).label,": p= "));
        fprintf(txtFid,strcat("%5e","\n"),comparisons(k).pValue);
    end
    for m=1:numel(overallEffects)
        fprintf(txtFid,'%s',strcat(indent,overallEffects(m).label,": p= "));
        fprintf(txtFid,strcat("%5e","\n"),overallEffects(m).pValue);
    end
end

function writeHeatmapCellLines(txtFid, indent, columnSpecs, DisplayVal, PValue, rowIdx)
% writeHeatmapCellLines  Writes one "<indent><column label>: value= <v>, p=
% <p>" line per column in this row, skipping columns with no data (NaN) --
% the shared line format for every age-fit-heatmap companion .txt file
% (2026-09-05, user request). value is DisplayVal (t*sign(Mu), the same
% quantity the heatmap cell's color encodes -- see gatherAgeHeatmapGrid),
% not the raw t-statistic, so it reads consistently with the plot.
    for c=1:numel(columnSpecs)
        v=DisplayVal(rowIdx,c);
        if isnan(v)
            continue
        end
        fprintf(txtFid,'%s',strcat(indent,columnSpecs(c).Label,": value= "));
        fprintf(txtFid,strcat("%5e",", p= "),v);
        fprintf(txtFid,strcat("%5e","\n"),PValue(rowIdx,c));
    end
end

function plotFactorViolinGrid(gridData, MetricOrder, groupLevels, demohex, democon, demoorder, democheat, saveFile, titleStr)
% plotFactorViolinGrid  TEST/PILOT. One horizontal row per metric, each row
% a small sideways version of the same violin plot plotGroupViolin draws for
% that metric on its own -- the same weighted mean +/- SE marker, the same
% zero line, the same significance test -- just compact enough to scan
% every metric in one figure, each on its own (unstandardized) scale.
% Horizontal rather than the usual vertical violins so the figure grows
% downward (adding metrics) rather than sideways (adding groups), which
% scales much better when there can be a dozen-plus metrics.
%
% Built as a stack of separate small figures (one per row, plus one for the
% title), each saved and reloaded as an image, rather than one figure with
% subplots/tiledlayout -- the third-party violin() function doesn't render
% reliably inside a tiledlayout, but this is the exact figure-per-plot
% pattern already used (and proven) everywhere else in this file.
%
%   gridData    - cell array, one entry per metric in MetricOrder, each a
%                 struct with fields:
%                   values - {V1,V2,...}: raw subject-level values per
%                            group, in groupLevels order
%                   mu, sigma - numeric arrays, one weighted mean/SE per
%                            group, same order as groupLevels
%                   comparisons - struct array of (leftIdx,rightIdx,pValue):
%                            one significance bracket per entry, in the same
%                            leftIdx/rightIdx/pValue shape plotGroupViolin
%                            uses. Entries are packed into lanes stepping
%                            right from the plot: two brackets whose group
%                            spans never touch or overlap (e.g. Male CN-DS
%                            and Female CN-DS) share one lane "at the same
%                            level," since they can never visually collide;
%                            anything that would collide gets its own lane.
%                            This packing is computed once from the first
%                            row, so a given lane is the same pairwise
%                            comparison(s) in every row -- order entries so
%                            related/adjacent-group comparisons come first
%                            if you want them packed together. Each
%                            bracket's line spans only its own leftIdx/
%                            rightIdx group positions, not the full row.
%                            Non-significant entries (sigStars(pValue)=="")
%                            are skipped, leaving that spot blank for that
%                            row. A 2-group comparison needs just one entry.
%                 An empty entry (metric not populated) is skipped.
%   groupLevels - string array matching demohex/democheat entries
%   saveFile    - full output path, no extension

    keep=~cellfun(@isempty,gridData);
    gridData=gridData(keep);
    MetricOrder=MetricOrder(keep);
    nMetrics=numel(MetricOrder);
    nGroups=numel(groupLevels);

    if nMetrics==0
        fprintf("Skipping violin grid '%s': no data for any metric.\n",titleStr);
        return
    end

    % One companion text file per grid image (2026-09-05, user request),
    % same name/folder as the .png via saveFile: every bracket/box p-value
    % that went into this grid, organized by metric row, so the numbers
    % behind the stars are available without eyeballing them off the plot.
    txtFid=fopen(strcat(saveFile,".txt"),'wt');
    fprintf(txtFid,'%s\n',titleStr);

    baseFontSize=11; % ylabel/xlabel/tick numbers; title and sig stars scale off this too
    starFontSize=max(baseFontSize-3,7); % small enough that multiple lanes fit side by side
    rowWidthIn=7.5;
    topMarginIn=0.08;
    plotHeightIn=violinRowHeightIn(nGroups); % item 12, 2026-09-05; must stay well
                        % above MATLAB's fixed axis-label/margin reservations, or
                        % the drawable plot area collapses to near-zero height
                        % and nothing real gets rendered
    tickMarginIn=0.28; % room for the x-tick numbers, every row
    xlabelExtraIn=0.30; % extra room below that, only on rows with an xlabel
    rightMarginIn=0.05;
    bracketGapIn=0.1;
    bracketTickLenIn=0.05; % length of each bracket's end-ticks, pointing back at the plot

    % Pack comparisons into lanes (computed once from the first row, since
    % every row shares the same comparisons shape) -- see planComparisonLanes.
    [laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn]=planComparisonLanes(gridData{1}.comparisons);

    % Reserve only as much bracket-column width as the row that actually
    % needs the most, not planComparisonLanes' own worst case (every
    % comparison significant at once) -- item 9, 2026-09-05. Frees up
    % plotWidthIn on every row with fewer active brackets than that.
    bracketColWidthIn=0;
    for m=1:numel(gridData)
        [~,~,w]=packActiveLanes(gridData{m}.comparisons,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn);
        bracketColWidthIn=max(bracketColWidthIn,w);
    end

    % Grid significance box, first pass (2026-09-05, user request): only
    % reserve extra width when at least one row actually carries
    % row.overallEffects (i.e. this grid is built from ANOVA GridData, not
    % T-Test GridData, which has no such field) -- so T-Test-only grids are
    % completely unaffected. boxWidthIn is a fixed first-pass guess sized
    % for the condensed single-line format ("Impair.* Sex n.s. Inter.**");
    % the user plans to review and correct sizing/placement once she sees
    % real output.
    hasBox=any(cellfun(@(r) isfield(r,'overallEffects') && ~isempty(r.overallEffects), gridData));
    boxWidthIn=0;
    if hasBox
        boxWidthIn=1.8;
    end
    % Follow-up, 2026-09-05: the box is meant to overlap the graph a lot
    % more than a fully separate column does -- only overhangFrac of the
    % box's own width should stick out past the row's original right edge,
    % with the rest of the box drawn on top of the existing plot/bracket
    % content (per user request).
    overhangFrac=0.075;

    % Measure the widest metric label (at the same font size the rows will
    % actually use) so the left margin is sized exactly to fit it -- no
    % more, no less -- freeing whatever's left over for the plot area
    % itself instead of sitting empty as excess margin.
    measureFig=figure('Visible','off');
    measureAx=axes(measureFig,'Units','inches');
    labelWidthsIn=zeros(1,nMetrics);
    for m=1:nMetrics
        t=text(measureAx,0,0,MetricOrder(m),'FontSize',baseFontSize,'Interpreter','none','Units','inches');
        labelWidthsIn(m)=t.Extent(3);
    end
    close(measureFig);
    % Cap the margin at the 75th-percentile label width rather than the
    % single widest one (2026-09-05): row labels already wrap to a second
    % line automatically (annotation textbox behavior) when they don't
    % fit, so letting just the longest handful wrap saves real width on
    % every other (shorter) row instead of sizing the whole column to fit
    % the single longest name.
    sortedWidthsIn=sort(labelWidthsIn);
    capIdx=max(1,ceil(0.75*numel(sortedWidthsIn)));
    maxLabelWidthIn=sortedWidthsIn(capIdx);
    leftMarginIn=maxLabelWidthIn+0.15;

    plotWidthIn=rowWidthIn-leftMarginIn-bracketGapIn-bracketColWidthIn-rightMarginIn;
    % The bracket column's fixed fallback start (leftMarginIn+plotWidthIn+
    % bracketGapIn) is no longer used directly -- each row now starts its
    % own bracket column dynamically, just past that row's own data (see
    % rowBracketStartIn below) -- but plotWidthIn above still reserves
    % enough total width for the worst case, since a row's dynamic start
    % can never land further right than this fixed position would.
    %
    % The sig box (if any) is anchored to the row's own original right edge
    % (rowWidthIn, still its pre-widened value here) so most of its width
    % overlaps the existing violin/bracket content instead of sitting in a
    % dedicated column -- only overhangFrac of it sticks out past that edge,
    % which is the only part rowWidthIn actually needs to grow by.
    boxStartIn=rowWidthIn-(1-overhangFrac)*boxWidthIn;
    if hasBox
        rowWidthIn=rowWidthIn+overhangFrac*boxWidthIn;
    end
    leftFrac=leftMarginIn/rowWidthIn;
    plotWidthFrac=plotWidthIn/rowWidthIn;

    tempDir=tempname;
    mkdir(tempDir);
    rowFiles=strings(nMetrics+1,1);
    % SigBoxOnly bonus composite (2026-09-05, user request): a second,
    % parallel row-file set -- same title row, but each metric row saved
    % with the violin/errorbar data hidden and the box shown instead of the
    % other way around -- stacked into its own composite alongside the
    % normal one. Only populated when hasBox, since a grid with no box has
    % nothing to make a bonus version out of.
    rowFilesBox=strings(nMetrics+1,1);

    % Row 0: a title-only figure, so the composite image has a caption.
    % Explicit full-figure axes (2026-09-05, consistency fix): a text()
    % call with no axes argument attaches to an auto-created DEFAULT axes,
    % which MATLAB insets from the figure edges -- coordinates given
    % without 'Units','normalized' additionally default to DATA units on
    % that axes' own auto XLim/YLim, not figure-normalized at all. A
    % dead-center point happens to still look roughly centered despite
    % this, which is why this title was never reported as visibly off,
    % but it's the same underlying issue found on plotCombinedViolinGrid/
    % plotSubnetAggregateViolinGrid's left-aligned headers.
    fig=figure('Visible','off');
    set(fig,'units','inches','position',[0 0 rowWidthIn 0.5],'windowstyle','normal')
    ax0=axes(fig,'Position',[0,0,1,1],'Visible','off');
    text(ax0,0.5,0.5,titleStr,'FontSize',baseFontSize+4,'HorizontalAlignment','center','Interpreter','none','Units','normalized')
    axis off
    rowFiles(1)=fullfile(tempDir,"row0.png");
    print(fig,rowFiles(1),'-dpng','-r200')
    close(fig)
    if hasBox
        % Identical title row, just under its own filename so the two
        % composites don't fight over (and delete) the same temp file.
        rowFilesBox(1)=fullfile(tempDir,"row0_box.png");
        copyfile(rowFiles(1),rowFilesBox(1));
    end

    for m=1:nMetrics
        row=gridData{m};
        fprintf(txtFid,'%s\n',MetricOrder(m));
        writeComparisonEffectLines(txtFid,"  ",row.comparisons,fieldOrEmpty(row,'overallEffects'));
        hasXLabel=(MetricOrder(m)=="Rate")||(m==nMetrics);
        bottomMarginIn=tickMarginIn;
        if hasXLabel
            % Rows that carry an xlabel need extra room below the tick
            % numbers so the label itself doesn't get clipped -- rather
            % than giving every row that room (wasting it on the rest),
            % only these rows' figures grow taller for it.
            bottomMarginIn=bottomMarginIn+xlabelExtraIn;
        end
        thisRowHeightIn=topMarginIn+plotHeightIn+bottomMarginIn;
        fig=figure('Visible','off');
        set(fig,'units','inches','position',[0 0 rowWidthIn thisRowHeightIn],'windowstyle','normal')
        ax=gca;
        % Fix the axes' own footprint in absolute inches (converted to the
        % normalized fractions Position needs) rather than leaving it to
        % MATLAB's auto layout, which reserves enough pixels for
        % labels/margins to swallow a short figure like this entirely.
        ax.Units='normalized';
        ax.Position=[leftFrac,bottomMarginIn/thisRowHeightIn,plotWidthFrac,plotHeightIn/thisRowHeightIn];
        ax.FontSize=baseFontSize; % tick numbers -- set explicitly rather than
                                   % relying on fontsize() inheritance timing

        boxHandle=drawViolinGridPanel(fig,ax,row,groupLevels,democon,demoorder,democheat,MetricOrder(m)=="Rate",nGroups, ...
            rowWidthIn,bracketGapIn,bracketTickLenIn,starFontSize,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn,5,baseFontSize,boxStartIn,boxWidthIn);
        % The box lives only on the SigBoxOnly bonus row below -- hidden
        % here on the original save (2026-09-05, user request, mirroring
        % plotGroupViolin's own box treatment).
        set(boxHandle,'Visible','off')

        % A native ylabel's anchor sits right at the axes' own edge (since
        % yticks are hidden, there's no y-tick-label width to push it
        % further left) -- 'HorizontalAlignment','center' centers the text
        % ON that anchor, so half of it spilled rightward over the axis
        % instead of sitting in the leftFrac margin column (2026-09-05,
        % follow-up: this differs from plotCombinedViolinGrid/
        % plotSubnetAggregateViolinGrid, which already draw their row
        % labels as a fixed-width annotation textbox spanning the WHOLE
        % margin column, so 'center' there centers within that column
        % instead). Switched to the same annotation textbox approach here
        % for the same centered-in-the-margin result.
        annotation(fig,'textbox',[0,bottomMarginIn/thisRowHeightIn,leftFrac,plotHeightIn/thisRowHeightIn], ...
            'String',MetricOrder(m),'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'EdgeColor','none','FontSize',baseFontSize,'Interpreter','none')
        if MetricOrder(m)=="Rate"
            xlabel(ax,"Rate (\rho)",'FontSize',baseFontSize)
        elseif m==nMetrics
            % Shared axis label for every non-Rate row, shown once on the
            % bottom-most row only -- the standard convention for a stacked
            % small-multiples grid where every panel shares the same units.
            xlabel(ax,"Factor Weight",'FontSize',baseFontSize)
        end

        rowFiles(m+1)=fullfile(tempDir,strcat("row",num2str(m),".png"));
        print(fig,rowFiles(m+1),'-dpng','-r200')
        if hasBox
            % SigBoxOnly bonus row: violin/errorbar/mean-line data hidden,
            % box shown instead -- brackets are untouched (annotation
            % objects on fig, not ax.Children) so they stay visible in both
            % saves, same as plotGroupViolin's own box/bracket treatment.
            dataChildren=ax.Children;
            set(dataChildren,'Visible','off')
            set(boxHandle,'Visible','on')
            rowFilesBox(m+1)=fullfile(tempDir,strcat("row",num2str(m),"_box.png"));
            print(fig,rowFilesBox(m+1),'-dpng','-r200')
        end
        close(fig)
    end

    fclose(txtFid);
    stackRowImages(rowFiles,saveFile);
    if hasBox
        stackRowImages(rowFilesBox,strcat(saveFile,"_SigBoxOnly"));
    end
    rmdir(tempDir);
end

function plotCombinedViolinGrid(gridDataF, gridDataS, MetricOrder, groupLevels, demohex, democon, demoorder, democheat, saveFile, titleStr)
% plotCombinedViolinGrid  Functional and Structural violin grids side by
% side in one figure: one shared column of metric-name labels on the left,
% a Functional/Structural header pair above the two panel columns, and one
% row per metric shared by both sides. A metric present on only one side
% (e.g. Structural has no "Assortativity", Functional has no "Same
% Subnet") still gets a row, with that side's panel left blank -- so the
% two columns stay aligned metric-for-metric instead of drifting apart.
%
% gridDataF, gridDataS - cell arrays, the same shape gridData takes in
%                         plotFactorViolinGrid, each indexed by the SAME
%                         MetricOrder position (an empty entry = that side
%                         has no data for that metric).
    keepF=~cellfun(@isempty,gridDataF);
    keepS=~cellfun(@isempty,gridDataS);
    keepEither=keepF|keepS;
    if ~any(keepEither)
        fprintf("Skipping combined violin grid '%s': no data for any metric on either side.\n",titleStr);
        return
    end
    gridDataF=gridDataF(keepEither);
    gridDataS=gridDataS(keepEither);
    hasF=keepF(keepEither);
    hasS=keepS(keepEither);
    MetricOrder=MetricOrder(keepEither);
    nMetrics=numel(MetricOrder);
    nGroups=numel(groupLevels);

    % One companion text file per grid image (2026-09-05, user request) --
    % see plotFactorViolinGrid's own txtFid for the full explanation.
    txtFid=fopen(strcat(saveFile,".txt"),'wt');
    fprintf(txtFid,'%s\n',titleStr);

    if hasF(1)
        templateComparisons=gridDataF{1}.comparisons;
    else
        templateComparisons=gridDataS{1}.comparisons;
    end
    [laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn]=planComparisonLanes(templateComparisons);

    % Reserve only as much bracket-column width as the row (checking BOTH
    % sides) that actually needs the most, not planComparisonLanes' own
    % worst case -- item 9, 2026-09-05.
    bracketColWidthIn=0;
    for m=1:nMetrics
        if hasF(m)
            [~,~,w]=packActiveLanes(gridDataF{m}.comparisons,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn);
            bracketColWidthIn=max(bracketColWidthIn,w);
        end
        if hasS(m)
            [~,~,w]=packActiveLanes(gridDataS{m}.comparisons,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn);
            bracketColWidthIn=max(bracketColWidthIn,w);
        end
    end

    baseFontSize=11;
    starFontSize=max(baseFontSize-3,7);
    panelTotalWidthIn=5; % each side's plot+bracket budget -- narrower than a standalone grid's
                          % full 7.5in (2026-09-05: with two panels side by side, that felt
                          % overly stretched out; more compact reads better)

    % Grid significance box, first pass (2026-09-05, user request): only
    % reserve the extra column when at least one row on either side
    % actually carries overallEffects (ANOVA GridData; T-Test GridData has
    % no such field), so a T-Test-only combined grid is unaffected.
    hasBox=false;
    for m=1:nMetrics
        if hasF(m) && isfield(gridDataF{m},'overallEffects') && ~isempty(gridDataF{m}.overallEffects)
            hasBox=true;
        end
        if hasS(m) && isfield(gridDataS{m},'overallEffects') && ~isempty(gridDataS{m}.overallEffects)
            hasBox=true;
        end
    end
    boxWidthIn=0;
    if hasBox
        boxWidthIn=1.8;
    end
    % Follow-up, 2026-09-05: the box is meant to overlap the graph a lot
    % more than a fully separate column does -- only overhangFrac of the
    % box's own width should stick out past the panel's original right
    % edge, with the rest drawn on top of the existing plot/bracket content.
    overhangFrac=0.075;
    panelFullWidthIn=panelTotalWidthIn+(hasBox*overhangFrac*boxWidthIn); % each side's total
                          % footprint, including its own box's small overhang when present --
                          % this is what actually spaces the two panels apart below, NOT
                          % panelTotalWidthIn directly
    topMarginIn=0.08;
    plotHeightIn=violinRowHeightIn(nGroups); % item 12, 2026-09-05
    tickMarginIn=0.28;
    xlabelExtraIn=0.30;
    rightMarginIn=0.05;
    bracketGapIn=0.1;
    bracketTickLenIn=0.05;
    panelGapIn=0.35; % gap between the Functional and Structural panels

    % Measure the widest metric label (checking BOTH sides, since a metric
    % only present on one side still needs to fit) so the shared left
    % margin fits it exactly.
    measureFig=figure('Visible','off');
    measureAx=axes(measureFig,'Units','inches');
    labelWidthsIn=zeros(1,nMetrics);
    for m=1:nMetrics
        t=text(measureAx,0,0,MetricOrder(m),'FontSize',baseFontSize,'Interpreter','none','Units','inches');
        labelWidthsIn(m)=t.Extent(3);
    end
    close(measureFig);
    % Cap the margin at the 75th-percentile label width rather than the
    % single widest one (2026-09-05): row labels already wrap to a second
    % line automatically (annotation textbox behavior) when they don't
    % fit, so letting just the longest handful wrap saves real width on
    % every other (shorter) row instead of sizing the whole column to fit
    % the single longest name.
    sortedWidthsIn=sort(labelWidthsIn);
    capIdx=max(1,ceil(0.75*numel(sortedWidthsIn)));
    maxLabelWidthIn=sortedWidthsIn(capIdx);
    leftMarginIn=maxLabelWidthIn+0.15;

    plotWidthIn=panelTotalWidthIn-bracketGapIn-bracketColWidthIn-rightMarginIn;
    figWidthIn=leftMarginIn+2*panelFullWidthIn+panelGapIn;
    fPanelLeftIn=leftMarginIn;
    sPanelLeftIn=leftMarginIn+panelFullWidthIn+panelGapIn;
    plotWidthFrac=plotWidthIn/figWidthIn;
    % Each side's box is anchored to that panel's own original right edge
    % (panelTotalWidthIn) so most of its width overlaps the existing
    % violin/bracket content instead of sitting in a dedicated column --
    % only overhangFrac of it sticks out past that edge, same reasoning as
    % plotFactorViolinGrid's boxStartIn.
    fBoxStartIn=fPanelLeftIn+panelTotalWidthIn-(1-overhangFrac)*boxWidthIn;
    sBoxStartIn=sPanelLeftIn+panelTotalWidthIn-(1-overhangFrac)*boxWidthIn;

    tempDir=tempname;
    mkdir(tempDir);
    rowFiles=strings(nMetrics+1,1);
    % SigBoxOnly bonus composite (2026-09-05, user request): see
    % plotFactorViolinGrid's own rowFilesBox for the full explanation --
    % same idea, just with both panels' data hidden per row instead of one.
    rowFilesBox=strings(nMetrics+1,1);

    % Row 0: title, plus "Functional"/"Structural" left-aligned at the
    % start of their respective panels (item 12, prior round, matching the
    % heatmap panel-label convention -- centering on just plotWidthIn/2
    % sat visibly off-center of the FULL panel once the bracket lane/
    % right-margin width was included); bolded, and the row shrunk/pulled
    % down so there's less dead space above the plots below. Row height
    % and title/header y raised again (item 8, 2026-09-05: the prior
    % round's shrink pulled the title down far enough to overlap the
    % header). The title's x centering was ALSO still wrong (2026-09-05,
    % item 6): figWidthIn (leftMarginIn+2*panelTotalWidthIn+panelGapIn)
    % includes each panel's FULL allotment, but the actually-visible
    % violin content is only plotWidthIn wide within each panel
    % (panelTotalWidthIn also reserves bracketGapIn+bracketColWidthIn+
    % rightMarginIn per side for the bracket lane, which isn't part of the
    % visible plot) -- so centering on [leftMarginIn,figWidthIn] put the
    % midpoint well to the right of the true content center (matching the
    % report that Structural, the rightmost panel, looked closer to
    % aligned than Functional). Centering on the true span instead
    % (leftMarginIn through the Structural panel's own real right edge)
    % fixes both sides at once.
    fig=figure('Visible','off');
    % Row0 grown (0.62->0.85in) and the header y raised (0.2->0.45) --
    % 2026-09-05 follow-up: at y=0.2 in a 0.62in-tall row, "Functional"/
    % "Structural"'s own text height put its bottom edge right at (or
    % past) row0's own bottom border, so it was getting clipped right
    % where row0 meets the first data row below -- reading as the header
    % "hiding behind" that row's plot. Title also nudged up a hair so the
    % now-higher headers still clear it.
    set(fig,'units','inches','position',[0 0 figWidthIn 0.85],'windowstyle','normal')
    % A plain text() call with no axes argument attaches to an
    % auto-created DEFAULT axes, which MATLAB insets from the figure edges
    % (roughly the classic [0.13 0.11 0.775 0.815] box) -- so
    % 'Units','normalized' was never actually relative to the whole
    % figure, it was relative to that smaller inset box. A left-aligned
    % label is very sensitive to that (the whole inset offset shows up
    % directly), while the centered title happened to look close anyway
    % (near-center values are much less affected by a roughly-symmetric
    % inset) -- which is why only "Functional"/"Structural" looked
    % noticeably pushed over (2026-09-05 follow-up). Fixed by creating an
    % explicit axes spanning the true full figure and using it for every
    % text() call below.
    ax0=axes(fig,'Position',[0,0,1,1],'Visible','off');
    trueContentRightIn=sPanelLeftIn+plotWidthIn;
    titleXFrac=(leftMarginIn+trueContentRightIn)/(2*figWidthIn);
    text(ax0,titleXFrac,0.97,titleStr,'FontSize',baseFontSize+4,'HorizontalAlignment','center','VerticalAlignment','top','Interpreter','none','Units','normalized')
    text(ax0,fPanelLeftIn/figWidthIn,0.45,"Functional",'FontSize',baseFontSize+2,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top','Interpreter','none','Units','normalized')
    text(ax0,sPanelLeftIn/figWidthIn,0.45,"Structural",'FontSize',baseFontSize+2,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top','Interpreter','none','Units','normalized')
    axis off
    rowFiles(1)=fullfile(tempDir,"row0.png");
    print(fig,rowFiles(1),'-dpng','-r200')
    close(fig)
    if hasBox
        rowFilesBox(1)=fullfile(tempDir,"row0_box.png");
        copyfile(rowFiles(1),rowFilesBox(1));
    end

    for m=1:nMetrics
        fprintf(txtFid,'%s\n',MetricOrder(m));
        if hasF(m)
            fprintf(txtFid,'  Functional:\n');
            writeComparisonEffectLines(txtFid,"    ",gridDataF{m}.comparisons,fieldOrEmpty(gridDataF{m},'overallEffects'));
        end
        if hasS(m)
            fprintf(txtFid,'  Structural:\n');
            writeComparisonEffectLines(txtFid,"    ",gridDataS{m}.comparisons,fieldOrEmpty(gridDataS{m},'overallEffects'));
        end
        hasXLabel=(MetricOrder(m)=="Rate")||(m==nMetrics);
        bottomMarginIn=tickMarginIn;
        if hasXLabel
            bottomMarginIn=bottomMarginIn+xlabelExtraIn;
        end
        thisRowHeightIn=topMarginIn+plotHeightIn+bottomMarginIn;
        fig=figure('Visible','off');
        set(fig,'units','inches','position',[0 0 figWidthIn thisRowHeightIn],'windowstyle','normal')
        boxHandles=gobjects(0);

        % Shared metric-name label, independent of either panel's own axes
        % so it always shows even when one side has no data to plot.
        annotation(fig,'textbox',[0,bottomMarginIn/thisRowHeightIn,leftMarginIn/figWidthIn-0.005,plotHeightIn/thisRowHeightIn], ...
            'String',MetricOrder(m),'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'EdgeColor','none','FontSize',baseFontSize,'Interpreter','none')

        if hasF(m)
            axF=axes(fig,'Units','normalized', ...
                'Position',[fPanelLeftIn/figWidthIn,bottomMarginIn/thisRowHeightIn,plotWidthFrac,plotHeightIn/thisRowHeightIn]);
            axF.FontSize=baseFontSize;
            boxHandleF=drawViolinGridPanel(fig,axF,gridDataF{m},groupLevels,democon,demoorder,democheat,MetricOrder(m)=="Rate",nGroups, ...
                figWidthIn,bracketGapIn,bracketTickLenIn,starFontSize,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn,5,baseFontSize,fBoxStartIn,boxWidthIn);
            boxHandles=[boxHandles,boxHandleF]; %#ok<AGROW>
            if MetricOrder(m)=="Rate"
                xlabel(axF,"Rate (\rho)",'FontSize',baseFontSize)
            elseif m==nMetrics
                xlabel(axF,"Factor Weight",'FontSize',baseFontSize)
            end
        elseif hasXLabel
            % Item 5, 2026-09-05: this metric has no Functional data (on
            % purpose), so there's no axes here for a normal xlabel to
            % attach to -- but the bottom row still needs its axis
            % caption, so place one directly where it would have sat had
            % the panel existed.
            placeMissingPanelXLabel(fig,fPanelLeftIn,figWidthIn,plotWidthFrac,bottomMarginIn,thisRowHeightIn,MetricOrder(m),baseFontSize)
        end

        if hasS(m)
            axS=axes(fig,'Units','normalized', ...
                'Position',[sPanelLeftIn/figWidthIn,bottomMarginIn/thisRowHeightIn,plotWidthFrac,plotHeightIn/thisRowHeightIn]);
            axS.FontSize=baseFontSize;
            boxHandleS=drawViolinGridPanel(fig,axS,gridDataS{m},groupLevels,democon,demoorder,democheat,MetricOrder(m)=="Rate",nGroups, ...
                figWidthIn,bracketGapIn,bracketTickLenIn,starFontSize,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn,5,baseFontSize,sBoxStartIn,boxWidthIn);
            boxHandles=[boxHandles,boxHandleS]; %#ok<AGROW>
            if MetricOrder(m)=="Rate"
                xlabel(axS,"Rate (\rho)",'FontSize',baseFontSize)
            elseif m==nMetrics
                xlabel(axS,"Factor Weight",'FontSize',baseFontSize)
            end
        elseif hasXLabel
            % Same as the Functional case above, for a metric missing on
            % the Structural side instead (item 5, 2026-09-05).
            placeMissingPanelXLabel(fig,sPanelLeftIn,figWidthIn,plotWidthFrac,bottomMarginIn,thisRowHeightIn,MetricOrder(m),baseFontSize)
        end

        % The box lives only on the SigBoxOnly bonus row below -- hidden
        % here on the original save (2026-09-05, user request).
        set(boxHandles,'Visible','off')
        rowFiles(m+1)=fullfile(tempDir,strcat("row",num2str(m),".png"));
        print(fig,rowFiles(m+1),'-dpng','-r200')
        if hasBox
            % SigBoxOnly bonus row: both panels' violin/errorbar/mean-line
            % data hidden, boxes shown instead -- brackets untouched (fig
            % annotations, not ax.Children) so they stay visible in both.
            dataChildren=gobjects(0);
            if hasF(m)
                dataChildren=[dataChildren;axF.Children]; %#ok<AGROW>
            end
            if hasS(m)
                dataChildren=[dataChildren;axS.Children]; %#ok<AGROW>
            end
            set(dataChildren,'Visible','off')
            set(boxHandles,'Visible','on')
            rowFilesBox(m+1)=fullfile(tempDir,strcat("row",num2str(m),"_box.png"));
            print(fig,rowFilesBox(m+1),'-dpng','-r200')
        end
        close(fig)
    end

    fclose(txtFid);
    stackRowImages(rowFiles,saveFile);
    if hasBox
        stackRowImages(rowFilesBox,strcat(saveFile,"_SigBoxOnly"));
    end
    rmdir(tempDir);
end

function plotSubnetAggregateViolinGrid(gridDataBySubnet, subnetOrder, MetricOrder, groupLevels, demohex, democon, demoorder, democheat, saveFile, titleStr)
% plotSubnetAggregateViolinGrid  One comparison's violin grid aggregated
% across every subnetwork, for ONE scope (Functional or Structural -- call
% this once per scope, kept separate, since each scope has a different and
% much smaller metric set at the subnetwork level than at whole-brain).
% Rows are subnetworks, columns are metrics -- the opposite of
% plotFactorViolinGrid's own row axis, because at the subnetwork level
% there are more subnetworks than metrics (whole-brain-only metrics like
% Assortativity or the Same-Subnet family aren't computed per subnetwork),
% so subnetworks get the axis that's free to grow. Each cell is the same
% drawViolinGridPanel violin+bracket panel plotFactorViolinGrid and
% plotCombinedViolinGrid already use, just tiled into a 2-D grid instead
% of one column (plotFactorViolinGrid) or two (plotCombinedViolinGrid).
%
%   gridDataBySubnet - struct, one field per subnetwork (field name =
%                       subnetOrder entry), each holding a gridData cell
%                       array in plotFactorViolinGrid's shape, indexed
%                       1:N against MetricOrder(1:N) for whatever N that
%                       subnetwork's data actually has (every subnetwork
%                       in the same scope shares the same N and the same
%                       MetricOrder(1:N) meaning -- only whole-brain adds
%                       the extra whole-brain-only metrics beyond N).
%   subnetOrder      - string array, the canonical row order; only entries
%                       actually present as a field are kept
    hasSubnet=isfield(gridDataBySubnet,cellstr(subnetOrder));
    subnetOrder=subnetOrder(hasSubnet);
    nSubnets=numel(subnetOrder);
    if nSubnets==0
        fprintf("Skipping subnetwork-aggregate violin grid '%s': no subnetwork data.\n",titleStr);
        return
    end

    % Every subnetwork in this scope shares the same metric-column count;
    % take the max defensively in case one came up short.
    nMetricsThisScope=0;
    for s=1:nSubnets
        nMetricsThisScope=max(nMetricsThisScope,numel(gridDataBySubnet.(subnetOrder(s))));
    end
    keepMetric=false(1,nMetricsThisScope);
    for s=1:nSubnets
        gd=gridDataBySubnet.(subnetOrder(s));
        thisKeep=~cellfun(@isempty,gd);
        keepMetric(1:numel(thisKeep))=keepMetric(1:numel(thisKeep))|thisKeep;
    end
    keptIdx=find(keepMetric); % maps kept-column position -> original MetricOrder index
    MetricOrder=MetricOrder(keptIdx);
    nMetrics=numel(MetricOrder);
    nGroups=numel(groupLevels);
    if nMetrics==0
        fprintf("Skipping subnetwork-aggregate violin grid '%s': no metric data for any subnetwork.\n",titleStr);
        return
    end

    % One companion text file per grid image (2026-09-05, user request) --
    % see plotFactorViolinGrid's own txtFid for the full explanation.
    txtFid=fopen(strcat(saveFile,".txt"),'wt');
    fprintf(txtFid,'%s\n',titleStr);

    % Lane packing (brackets) is the same shape for every cell -- same
    % comparisons struct regardless of subnetwork or metric -- so it's
    % computed once from the first populated cell found.
    templateComparisons=[];
    for s=1:nSubnets
        gd=gridDataBySubnet.(subnetOrder(s));
        for c=1:nMetrics
            if keptIdx(c)<=numel(gd) && ~isempty(gd{keptIdx(c)})
                templateComparisons=gd{keptIdx(c)}.comparisons;
                break
            end
        end
        if ~isempty(templateComparisons)
            break
        end
    end
    [laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn]=planComparisonLanes(templateComparisons);

    % Reserve only as much bracket-column width as the cell (checking
    % every subnetwork x metric combination) that actually needs the
    % most, not planComparisonLanes' own worst case -- item 9, 2026-09-05.
    bracketColWidthIn=0;
    hasBox=false;
    for s=1:nSubnets
        gd=gridDataBySubnet.(subnetOrder(s));
        for c=1:nMetrics
            if keptIdx(c)<=numel(gd) && ~isempty(gd{keptIdx(c)})
                [~,~,w]=packActiveLanes(gd{keptIdx(c)}.comparisons,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn);
                bracketColWidthIn=max(bracketColWidthIn,w);
                if isfield(gd{keptIdx(c)},'overallEffects') && ~isempty(gd{keptIdx(c)}.overallEffects)
                    hasBox=true;
                end
            end
        end
    end

    baseFontSize=11;
    starFontSize=max(baseFontSize-3,7);
    plotWidthIn=2.2; % per-column violin plot area; generous rather than print-page-constrained,
                      % since this figure is meant to be viewed/zoomed digitally, not printed to a fixed page
    colGapIn=0.35; % gap between metric columns
    topMarginIn=0.08;
    rowHeightIn=violinRowHeightIn(nGroups); % item 12, 2026-09-05
    tickMarginIn=0.28;
    xlabelExtraIn=0.30; % extra room on the bottom row only
    rightMarginIn=0.05;
    bracketGapIn=0.1;
    bracketTickLenIn=0.05;

    colWidthIn=bracketGapIn+bracketColWidthIn+plotWidthIn+rightMarginIn;
    % Grid significance box, first pass (2026-09-05, user request): only
    % reserve extra width when at least one cell in this scope actually
    % carries overallEffects (ANOVA GridData; T-Test GridData has no such
    % field), so a T-Test-only aggregate is unaffected.
    boxWidthIn=0;
    if hasBox
        boxWidthIn=1.8;
    end
    % Follow-up, 2026-09-05: the box is meant to overlap the graph a lot
    % more than a fully separate column does -- only overhangFrac of the
    % box's own width should stick out past the column's original right
    % edge, with the rest drawn on top of the existing plot/bracket content.
    overhangFrac=0.075;
    colFullWidthIn=colWidthIn+(hasBox*overhangFrac*boxWidthIn);
    colTotalWidthIn=colFullWidthIn+colGapIn;

    % Subnetwork row-name labels get their own (possibly smaller) font size
    % when there are many of them, and a larger safety pad beyond the
    % measured text width -- item 15, 2026-09-05. The row labels are drawn
    % via annotation('textbox',...) below, which carries its own internal
    % margin/padding not accounted for by a plain text() extent
    % measurement, so the plain-text-based margin alone could still clip
    % long names.
    if nSubnets>8
        subnetLabelFontSize=baseFontSize-1;
    else
        subnetLabelFontSize=baseFontSize;
    end

    % Measure the widest subnetwork label so the shared left margin fits
    % it exactly -- same approach every other grid in this file uses.
    measureFig=figure('Visible','off');
    measureAx=axes(measureFig,'Units','inches');
    maxLabelWidthIn=0;
    for s=1:nSubnets
        t=text(measureAx,0,0,subnetOrder(s),'FontSize',subnetLabelFontSize,'Interpreter','none','Units','inches');
        maxLabelWidthIn=max(maxLabelWidthIn,t.Extent(3));
    end
    close(measureFig);
    leftMarginIn=maxLabelWidthIn+0.3; % wider safety pad (was +0.15) to cover the
                                       % annotation textbox's own margin

    figWidthIn=leftMarginIn+nMetrics*colTotalWidthIn-colGapIn;
    colLeftIn=leftMarginIn+(0:nMetrics-1)*colTotalWidthIn;
    plotWidthFrac=plotWidthIn/figWidthIn;
    % Each column's box is anchored to that column's own original right
    % edge (colWidthIn) so most of its width overlaps the existing
    % violin/bracket content instead of sitting in a dedicated column --
    % only overhangFrac of it sticks out past that edge, same reasoning as
    % plotFactorViolinGrid's boxStartIn.
    colBoxStartIn=colLeftIn+colWidthIn-(1-overhangFrac)*boxWidthIn;

    tempDir=tempname;
    mkdir(tempDir);
    rowFiles=strings(nSubnets+1,1);
    % SigBoxOnly bonus composite (2026-09-05, user request): see
    % plotFactorViolinGrid's own rowFilesBox for the full explanation --
    % same idea, just with every column's data hidden per row instead of one.
    rowFilesBox=strings(nSubnets+1,1);

    % Row 0: title, plus one metric-name header per column. Headers are
    % left-aligned at each column's own start (item 12, prior round,
    % matching the heatmap panel-label convention) rather than centered
    % on just the plot area (plotWidthIn/2), which sat visibly off-center
    % of the FULL column once the bracket lane/right-margin width was
    % included; bolded, and the row shrunk/pulled down so there's less
    % dead space between the headers and the plots below. Row height and
    % title/header y raised again (item 8, prior round: the earlier
    % round's shrink pulled the title down far enough to overlap the
    % headers). The title's x centering was ALSO still wrong (2026-09-05,
    % item 6): figWidthIn (leftMarginIn+nMetrics*colTotalWidthIn-colGapIn)
    % includes every column's FULL allotment, but the actually-visible
    % violin content per column is only plotWidthIn wide (colTotalWidthIn
    % also reserves bracketGapIn+bracketColWidthIn+rightMarginIn+colGapIn,
    % which isn't part of the visible plot) -- so centering on
    % [leftMarginIn,figWidthIn] put the midpoint well to the right of the
    % true content center (matching the report that a middle-to-last
    % column looked closer to aligned than the true center). Centering on
    % the true span instead (leftMarginIn through the last column's own
    % real right edge) fixes it.
    fig=figure('Visible','off');
    % Row0 grown (0.62->0.85in) and the header y raised (0.2->0.45) --
    % 2026-09-05 follow-up: at y=0.2 in a 0.62in-tall row, each column
    % header's own text height put its bottom edge right at (or past)
    % row0's own bottom border, so it was getting clipped right where
    % row0 meets the first subnetwork row below -- reading as the header
    % "hiding behind" that row's plot. Title also nudged up a hair so the
    % now-higher headers still clear it.
    set(fig,'units','inches','position',[0 0 figWidthIn 0.85],'windowstyle','normal')
    % A plain text() call with no axes argument attaches to an
    % auto-created DEFAULT axes, which MATLAB insets from the figure edges
    % (roughly the classic [0.13 0.11 0.775 0.815] box) -- so
    % 'Units','normalized' was never actually relative to the whole
    % figure, it was relative to that smaller inset box. A left-aligned
    % label is very sensitive to that (the whole inset offset shows up
    % directly), while the centered title happened to look close anyway
    % (near-center values are much less affected by a roughly-symmetric
    % inset) -- which is why only the per-column headers looked noticeably
    % pushed over, not the title (2026-09-05 follow-up). Fixed by creating
    % an explicit axes spanning the true full figure and using it for
    % every text() call below.
    ax0=axes(fig,'Position',[0,0,1,1],'Visible','off');
    trueContentRightIn=colLeftIn(nMetrics)+plotWidthIn;
    titleXFrac=(leftMarginIn+trueContentRightIn)/(2*figWidthIn);
    text(ax0,titleXFrac,0.97,titleStr,'FontSize',baseFontSize+4,'HorizontalAlignment','center','VerticalAlignment','top','Interpreter','none','Units','normalized')
    for c=1:nMetrics
        text(ax0,colLeftIn(c)/figWidthIn,0.45,MetricOrder(c),'FontSize',baseFontSize+1,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top','Interpreter','none','Units','normalized')
    end
    axis off
    rowFiles(1)=fullfile(tempDir,"row0.png");
    print(fig,rowFiles(1),'-dpng','-r200')
    close(fig)
    if hasBox
        rowFilesBox(1)=fullfile(tempDir,"row0_box.png");
        copyfile(rowFiles(1),rowFilesBox(1));
    end

    for s=1:nSubnets
        gd=gridDataBySubnet.(subnetOrder(s));
        fprintf(txtFid,'%s\n',subnetOrder(s));
        hasXLabel=(s==nSubnets);
        bottomMarginIn=tickMarginIn;
        if hasXLabel
            bottomMarginIn=bottomMarginIn+xlabelExtraIn;
        end
        thisRowHeightIn=topMarginIn+rowHeightIn+bottomMarginIn;
        fig=figure('Visible','off');
        set(fig,'units','inches','position',[0 0 figWidthIn thisRowHeightIn],'windowstyle','normal')
        boxHandles=gobjects(0);
        rowAxes=gobjects(0);

        % Shared subnetwork-name label, independent of any single column's
        % own axes so it always shows even if every column here happened
        % to be empty for this subnetwork.
        annotation(fig,'textbox',[0,bottomMarginIn/thisRowHeightIn,leftMarginIn/figWidthIn-0.005,rowHeightIn/thisRowHeightIn], ...
            'String',subnetOrder(s),'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'EdgeColor','none','FontSize',subnetLabelFontSize,'Interpreter','none')

        for c=1:nMetrics
            origIdx=keptIdx(c);
            if origIdx>numel(gd) || isempty(gd{origIdx})
                continue
            end
            fprintf(txtFid,'  %s\n',MetricOrder(c));
            writeComparisonEffectLines(txtFid,"    ",gd{origIdx}.comparisons,fieldOrEmpty(gd{origIdx},'overallEffects'));
            ax=axes(fig,'Units','normalized', ...
                'Position',[colLeftIn(c)/figWidthIn,bottomMarginIn/thisRowHeightIn,plotWidthFrac,rowHeightIn/thisRowHeightIn]);
            ax.FontSize=baseFontSize;
            % Fewer, smaller x-axis ticks here than the other grid callers
            % (item 11, 2026-09-05): this panel is much narrower
            % (plotWidthIn) than the other grids' rowWidthIn/
            % panelTotalWidthIn, so 5 ticks at the normal font size
            % overlapped and got auto-rotated diagonal by MATLAB, running
            % past the panel and clipping.
            boxHandle=drawViolinGridPanel(fig,ax,gd{origIdx},groupLevels,democon,demoorder,democheat,MetricOrder(c)=="Rate",nGroups, ...
                figWidthIn,bracketGapIn,bracketTickLenIn,starFontSize,laneOfComparison,nLanesActual,laneWidths,compactLaneWidthIn,wideLaneWidthIn,3,max(baseFontSize-3,7),colBoxStartIn(c),boxWidthIn);
            boxHandles=[boxHandles,boxHandle]; %#ok<AGROW>
            rowAxes=[rowAxes,ax]; %#ok<AGROW>
            if hasXLabel
                if MetricOrder(c)=="Rate"
                    xlabel(ax,"Rate (\rho)",'FontSize',baseFontSize)
                else
                    xlabel(ax,"Factor Weight",'FontSize',baseFontSize)
                end
            end
        end

        % The box lives only on the SigBoxOnly bonus row below -- hidden
        % here on the original save (2026-09-05, user request).
        set(boxHandles,'Visible','off')
        rowFiles(s+1)=fullfile(tempDir,strcat("row",num2str(s),".png"));
        print(fig,rowFiles(s+1),'-dpng','-r200')
        if hasBox
            % SigBoxOnly bonus row: every column's violin/errorbar/
            % mean-line data hidden, boxes shown instead -- brackets
            % untouched (fig annotations, not ax.Children) so they stay
            % visible in both.
            dataChildren=gobjects(0);
            for a=1:numel(rowAxes)
                dataChildren=[dataChildren;rowAxes(a).Children]; %#ok<AGROW>
            end
            set(dataChildren,'Visible','off')
            set(boxHandles,'Visible','on')
            rowFilesBox(s+1)=fullfile(tempDir,strcat("row",num2str(s),"_box.png"));
            print(fig,rowFilesBox(s+1),'-dpng','-r200')
        end
        close(fig)
    end

    fclose(txtFid);
    stackRowImages(rowFiles,saveFile);
    if hasBox
        stackRowImages(rowFilesBox,strcat(saveFile,"_SigBoxOnly"));
    end
    rmdir(tempDir);
end

