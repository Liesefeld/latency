function [res,cfgNew] = latency(cfg,avgs,sign)
% [res, cfgNew] = latency(cfg,avgs,sign)
% cfg and avgs are mandatory, cfg can be empty if sign is entered as third
%   argument
% avgs is a subjects x channels x time array, contains such an array in the 
%   field 'individual' or 'data' or as specified in cfg.datafield, or
%   contains a field 'avg' which is a cell array of structures (one
%   structure for each subject)
% sign can be indicated as cfg.sign or as the third input argument
%  (argument overrules struct field); it is either 'pos'/1 or 'neg'/-1
% all timining information is in samples or in the units of cfg.times
% optional fields for cfg are (all indices as in avgs.individual):
%   extract:     which measures to extract; see the description of the
%                output [res] below
%   aggregate:   individual, GA (= grand average), of one of two jackknife
%                estimates, namely jackMiller (Miller, Patterson,& Ulrich, 
%                1998, Psycholphysiology), or jackSmulders(Smulders, 2010, 
%                Psychophysiology)
%   fig:         if true, individual ERPs and estimates are plotted
%   subs:        subject included; indices or, if cfg.subNum is given,
%                subject numbers(default = all)
%   subNum:      list of subject numbers (same order as in avgs)
%   chans:       channel included into the average; indices or, if
%                cfg.chanName is given, channel names(default = all)
%   chanName:    list of channel names (same order as in avgs)
%   peakWin:     search window for peak detection [start end]; indices or,
%                if cfg.times is given, time points(default = full range);
%                can also be a subjects x 2 matrix for individual search
%                ranges
%   meanTime:    time range across which activity is averaged for mean
%                (default = peakWin); can also be a subjects x 2 matrix
%                for individual mean windows
%   times:       all individual time points OR start, end and sampling rate
%                of the submitted data;
%                if set, temporal inputs (peakWidth, peakWin, meanTime)
%                are interpreted in time units instead of sampling points and
%                latency estimates (peakLat, onset, offset, width, areaLat
%                counterLat) are returned in time units
%   peakWidth:   time around peak averaged for peak amplitude; time range = 
%                peak time +- peakWidth; default = 5 sampling points or ms))
%   cWinStart:   time point from which on the counter peak is searched
%                ('peak'  or 'peakWin' (default))
%   cWinWidth:   width of the counter-peak search interval starting from
%                cfg.cWinStart (negative numbers for counter peaks
%                preceding the actual peak; indices or times);
%   cWin:        alternative to cfg.cWinStart and cfg.cWinWidth: search
%                for the counter peak [start end]
%   percArea:    percentage of total area for percent-area latency
%                (default = 0.5)
%   percAmp:     percentage of the amplitude (peak-to-peak if counterpeak
%                is used);
%   areaWin:     determines temporal boundaries for the area calculation
%                'peakWin' (default), 'ampLat', 'fullRange', or [start end]
%   areaBase:    determines one boundary in amplitude space (the other is
%                determined by the component; 'zero' (default) or 'percAmp'
%   ampLatWin:   where on- and offsets are searched: 'fullRange' (default),
%                'peakWin' or [start end]
%   cBound:      Boolean that indicats whether the counter peak is one
%                border for onset/offset
%                
% res contains a (latency) measure for each subject (is a struct if 
%   several measures are requested); latencies are indices or, if cfg.time
%   is given, time points
%   mean:        mean over cfg.meanTime;
%   peakLat:     latency of the detected peak, which can also be at one border
%   onset:       time point where the ERP has reached cfg.percAmp
%   offset:      time point where the ERP has fallen back to cfg.percAmp
%   width:       difference between on- and offset
%   areaLat:     percent-area latency
%   peakAmp:     mean amplitude in the window peakLat+-peakWidth;
%   peak2peak:   difference in peak amplitude between searched peak and
%                counter peak
%   baseline:    amplitude cfg.percAmp% of the peak amplitude or in between 
%                the two peaks if counter peak is used
%   area:        sum of values above/below the cfg.areaBase
%   foundLocal:  Boolean that indicates whether a local peak was found;
%                if false, peakLat is one boundary of the search interval
%   foundOn:     Boolean that indicates whether start of the component
%                (ERP crossed the baseline) was found
%   foundOff:    Boolean that indicates whether end of the component was found
%   foundOff:    Boolean that indicates whether area latency was extracted
%                successfully
%   counterLat:  same as peakLat for the counter peak
%   counterAmp:  same as peakAmp for the counter peak
%   foundLocalC: same as foundLocal for the counter peak
% cfgNew is the same as the input cfg plus default parameters for not
% indicated parameters and parameter-sensitive explanations for each field
% of res in cfgNew.ann
%
% Typical problems:
% - instead of finding a local peak - the algorithm returns a global maximum/
%   minimum and res.foundLocal is false; increasing cfg.peakWin should help
% - the component does not return to baseline before or after the peak -
%   the respective border defined by cfg.ampLatWin is selected as on- or
%   offset and res.foundOn or res.foundOff is false, respectively
%
% Originally created by Heinrich René Liesefeld in September 2014
% Last modified by Heinrich René Liesefeld in April 2018
%
% Please send bug reports to Heinrich.Liesefeld@psy.lmu.de

if isempty(cfg)
    cfg=struct;
end
if exist('sign','var')
    %sign of the component can also be supplied as third optional input argument
    cfg.sign=sign;
end
if isstruct(avgs)
    if isfield(cfg,'datafield')
        avgs=avgs.(cfg.datafield);
    elseif isfield(avgs,'individual')
        avgs=avgs.individual;
    elseif isfield(avgs,'data')
        avgs=avgs.data;
    else
        error('avgs must either be a subject x electrode x time array or contain a field "individual" (Fieldtrip) or "data" (EEGlab) or a user-specified field (cfg.datafield), which is such an array.')
    end
elseif iscell(avgs)
    avgsOld=avgs;
    nSubs=length(avgs);
    avgs=nan([nSubs,size(avgs{1}.avg)]);
    for subi=1:nSubs
        avgs(subi,:,:) = avgsOld{subi}.avg;
    end
end
cfg=checkConfig(cfg,avgs);

if strcmp(cfg.aggregate,'GA')
    nSubs=1;
elseif strcmp(cfg.aggregate,'jackMiller')
    nSubs=length(cfg.subs)+1;
elseif ismember(cfg.aggregate,{'individual','jackSmulders'})
    nSubs=length(cfg.subs);
end


%preallocate all variables
[comp.mean]=deal(nan(nSubs,1));
if cfg.get.peak
    [comp.peakAmp,comp.peakLat]=deal(nan(nSubs,1));
    comp.foundLocal=false(nSubs,1);
end
[comp.area,comp.baseline]=deal(nan(nSubs,1));
if cfg.get.ampLat
    [comp.onset,comp.offset,comp.width]=deal(nan(nSubs,1));
    [comp.foundOn,comp.foundOff]=deal(false(nSubs,1));
end
if cfg.get.counterPeak
    [comp.counterAmp,comp.counterLat,comp.peak2peak]=deal(nan(nSubs,1));
    comp.isLocalC=false(nSubs,1);
end
if cfg.get.areaLat
    [comp.areaLat,comp.foundArea]=deal(nan(nSubs,1));
end


if size(cfg.peakWin,1)>1 && size(cfg.peakWin,2)>1
    allPeakWin=cfg.peakWin;    
    %transpose array if it is the wrong way round
    if size(allPeakWin,1)<size(allPeakWin,2)
        allPeakWin=allPeakWin';
    end
end
if size(cfg.meanTime,1)>1 && size(cfg.meanTime,2)>1
    allMeanTimes=cfg.meanTime;    
    %transpose array if it is the wrong way round
    if size(allMeanTimes,1)<size(allMeanTimes,2)
        allMeanTimes=allMeanTimes';
    end
end


for subi=1:nSubs %cycle through subjects
    if strcmp(cfg.aggregate,'GA')
        subNum=1;
    elseif strcmp(cfg.aggregate,'jackMiller') && subi==nSubs
        subNum=99;
    else
        subNum=cfg.subNum(subi);
    end
    if exist('allPeakWin','var')
        cfg.peakWin=allPeakWin(subi,:);
    end
    switch cfg.aggregate
        case 'GA'
            ERP = squeeze(mean(avgs(cfg.subs,cfg.chans,:),1));
        case {'jackMiller','jackSmulders'}
            if strcmp(cfg.aggregate,'jackMiller') && subi==nSubs
                %last value is estimate from full grand average (the value
                %actually reported)
                ERP = squeeze(mean(avgs(cfg.subs,cfg.chans,:),1));
            else
                %other values are leave-one-out (used to estimate standard
                %error with the Miller approach)
                ERP = squeeze(mean(avgs(cfg.subs(cfg.subs~=subNum),cfg.chans,:),1));
            end
        case 'individual'
            ERP = squeeze(avgs(cfg.subs(subi),cfg.chans,:));
    end
    if numel(cfg.chans) > 1 %pool over several channels
        ERP=mean(ERP,1)';
    end
    
    %mean area
    if exist('allMeanTimes','var')
        comp.mean(subi) = mean(ERP(allMeanTimes(subi,1):allMeanTimes(subi,2)));
    else
        comp.mean(subi) = mean(ERP(cfg.meanTime(1):cfg.meanTime(2)));
    end
    peak=peakDetection(cfg,ERP);
    comp.peakLat(subi)=cfg.times(peak.lat);
    comp.foundLocal(subi)=peak.foundLocal;
    if ~peak.foundLocal
        warning('No local peak found for subject %g',subNum)
    end
    comp.peakAmp(subi)=peak.amp;
    
    if cfg.get.counterPeak
        if isfield(cfg,'cWin')
            cPeak=peakDetection(cfg,ERP,cfg.cWin);
        else
            cPeak=peakDetection(cfg,ERP,peak.lat);
        end
        comp.counterLat(subi)=cfg.times(cPeak.lat);
        comp.counterAmp(subi)=cPeak.amp;
        comp.peak2peak(subi)=peak.amp-cPeak.amp;
        comp.isLocalC(subi)=cPeak.foundLocal;
        if ~cPeak.foundLocal
            warning('No local counter-peak found for subject %g',subNum)
        end
    end
    
    
    
    %percent amplitude
    if cfg.get.counterPeak
        %calculate threshold relative to peak-to-peak amplitude when
        %counterPeak is used
        if cfg.sign==1
            comp.baseline(subi)=peak.amp-(peak.amp-cPeak.amp)*(1-cfg.percAmp);
        elseif cfg.sign==-1
            comp.baseline(subi)=peak.amp+(cPeak.amp-peak.amp)*(1-cfg.percAmp);
        end
    else
        %calculate threshold relative to component amplitude when
        %counterPeak is not used
        comp.baseline(subi)=peak.amp-peak.amp*(1-cfg.percAmp);
    end

    if cfg.get.ampLat
        if isfield(comp,'counterLat')
            ampLat=amplitudeLatency(cfg,ERP,comp.baseline(subi),peak.lat,cPeak.lat);
        else
            ampLat=amplitudeLatency(cfg,ERP,comp.baseline(subi),peak.lat);
        end
        if cfg.warnings
            if ~ampLat.foundOn
                warning('Could not find onset for subject %g. Increase cfg.percAmp or cfg.ampLatWin to solve this problem!',subNum)
            end
            if ~ampLat.foundOff
                warning('Could not find offset for subject %g. Increase cfg.percAmp or cfg.ampLatWin to solve this problem!',subNum)
            end
        end

        comp.foundOn(subi)=ampLat.foundOn;
        comp.foundOff(subi)=ampLat.foundOff;
        comp.onset(subi)=cfg.times(ampLat.start);
        comp.offset(subi)=cfg.times(ampLat.end);
        comp.width(subi)=(ampLat.end-ampLat.start)/cfg.sampRate;

    end
    if cfg.get.area
        %get sum of all values above/below the baseline and areaLat
        switch cfg.areaBase
            case 'percAmp'
                baseline=comp.baseline(subi);
            case 'zero'
                baseline=0;
        end
        if isnumeric(cfg.areaWin)
            time=cfg.areaWin;
        else
            switch cfg.areaWin
                case 'ampLat'
                    time=[ampLat.start,ampLat.end];
                case 'peakWin'
                    time=cfg.peakWin;
                case 'fullRange'
                    time=[1,length(ERP)];
                otherwise
                    error('cfg.areaWin must be either "ampLat", "peakWin", or "fullRange", it is %s',cfg.areaWin)
            end
        end
        comp.area(subi)=totalArea(cfg,ERP,baseline,time);
        if cfg.get.areaLat
            areaLat=areaLatency(cfg,ERP,comp.area(subi),baseline,time);
        end
    end
    if cfg.get.areaLat
        comp.areaLat(subi)=cfg.times(areaLat.lat);
        comp.foundArea(subi)=areaLat.foundLat;
        if ~areaLat.foundLat
            warning('Did not find fractional area latency for subject %g; set latency to average of the time range',subNum);
        end
    end
    
    
    
    if cfg.fig
        if subi==1
            figure('Name','Individual ERPs and extracted estimates');
        end
        ncols=ceil(nSubs/10);
        subplot(ceil(nSubs/ncols),ncols,subi);
        if isfield(cfg,'times')
            xRange=[cfg.times(1),cfg.times(end)];
        else
            xRange=[0 length(ERP)];
        end
        yRange=[min(ERP) max(ERP)];
        plot(cfg.times,ERP,'k');
        hold on
        plot(comp.peakLat(subi),comp.peakAmp(subi),'o');
        if cfg.get.counterPeak
            plot(comp.counterLat(subi),comp.counterAmp(subi),'ro');
        end
        plot(xRange,[comp.baseline(subi) comp.baseline(subi)]);
        if cfg.get.ampLat
            plot([comp.onset(subi) comp.onset(subi)],yRange);
            plot([comp.offset(subi) comp.offset(subi)],yRange);
        end
        if cfg.get.areaLat
            plot([comp.areaLat(subi) comp.areaLat(subi)],yRange,'g');
        end
        plot(cfg.times([cfg.peakWin(1) cfg.peakWin(1)]),yRange,'k--');
        plot(cfg.times([cfg.peakWin(2) cfg.peakWin(2)]),yRange,'k--');
        axis('tight')
        if ismember(cfg.aggregate,{'jackMiller','jackSmulders'})
            if strcmp(cfg.aggregate,'jackMiller') && subi==nSubs
                title('Average');
            else
                title(sprintf('left out subject #%g',subNum));
            end
        else
            title(sprintf('subject #%g',subNum));
        end
    end
    
end %end of loop over subjects

if exist('allPeakWin','var')
    cfg.peakWin=allPeakWin;
end

if mean(comp.foundLocal) < 0.5
    %check the proportion of local peaks found
    fprintf('Attention: In less than half of the cases a local peak was found.\n');
    fprintf('Consider increasing the search interval.\n');
end
if isfield(comp,'isLocalC') && mean(comp.isLocalC) < 0.5
    %check the proportion of local counter peaks found
    fprintf('Attention: In less than half of the cases a local counter peak was found.\n');
    fprintf('Consider increasing the search interval.\n');
end
if isfield(comp,'foundOn') && mean(comp.foundOn) <0.5
    fprintf('Attention: In less than half of the cases the onset of the component was found.\n');
    fprintf('Consider increasing the on-/offset search interval (cfg.ampLatWin), changing cfg.percAmp, and/or using the counter-peak method.\n');
end
if isfield(comp,'foundOff') && mean(comp.foundOff) <0.5
    fprintf('Attention: In less than half of the cases the offset of the component was found.\n');
    fprintf('Consider increasing the on-/offset search interval (cfg.ampLatWin), changing cfg.percAmp, and/or using the counter-peak method.\n');
end
if isfield(cfg,'chanName')
    cfg.chans=cfg.chanName(cfg.chans);
elseif any(ischar(cfg.chans)) || any(mod(cfg.chans,1))
    error('If cfg.chanName is not set, cfg.chans must contain indices, not channel names.')
end
if iscell(cfg.extract)
    for outi=1:length(cfg.extract)
        if strcmp(cfg.aggregate,'jackSmulders')
            %oi = n*mean(J) - (n - 1)*ji, see Smulders (2010)
            for subi=1:nSubs
                res.(cfg.extract{outi})(subi)=nSubs*mean(comp.(cfg.extract{outi}))-(nSubs-1)*comp.(cfg.extract{outi})(subi);
            end
        else
            res.(cfg.extract{outi})=comp.(cfg.extract{outi});
        end
    end
elseif strcmp(cfg.aggregate,'jackSmulders')
    %oi = n*mean(J) - (n - 1)*ji, see Smulders (2010)
    for subi=1:nSubs
        res(subi)=nSubs*mean(comp.(cfg.extract))-(nSubs-1)*comp.(cfg.extract)(subi);
    end
    
else
    res=comp.(cfg.extract);
end
cfgNew = cfg;
end

function cfg=checkConfig(cfg,avgs)
%cfg.sign must be indicated or function returns an error
%avgs is necessary when cfg.subs or cfg.chans or cfg.peakWin is
%missing (to determine the maximum extend in each of these dimensions)

if ~isfield(cfg,'sign')
    error('Indicating the sign of the component (cfg.sign = ''pos'' or 1 or ''neg'' or -1) is mandatory.')
end
if ischar(cfg.sign)
    switch cfg.sign
        case 'pos'
            cfg.sign=1;
        case 'neg'
            cfg.sign=-1;
    end
end

if ~isfield(cfg,'warnings'),cfg.warnings=1;end
if ~isfield(cfg,'subs') || strcmp(cfg.subs,'all')
    cfg.subs=1:size(avgs,1);
end
if isfield(cfg,'subs')
    if all(ismember(cfg.subs,[1,0]))
        cfg.subs = logical(cfg.subs);
    end
    if islogical(cfg.subs)
        %cfg.subs can also be set as a filter
        cfg.subs=find(cfg.subs);
    end
elseif isfield(cfg,'subNum')
    allSubs=1:length(cfg.subNum);
    cfg.subs=allSubs(ismember(cfg.subNum,cfg.subs));
end
if ~isfield(cfg,'subNum')
    cfg.subNum=cfg.subs;
end

if ~isfield(cfg,'chans')
    if cfg.warnings
        warning('No channels (cfg.chans) indicated; will average across all of them.')
    end
    cfg.chans=1:size(avgs,2);
elseif isfield(cfg,'chanName') && (iscell(cfg.chans) || ischar(cfg.chans))
    allChans=1:length(cfg.chanName);
    cfg.chans=allChans(ismember(cfg.chanName,cfg.chans));    
end
if ~isfield(cfg,'peakWin')
    if cfg.warnings
        warning('No time range (cfg.peakWin) indicated; will use full time range.');
    end
    cfg.peakWin=[1,size(avgs,3)];
end
if isfield(cfg,'times')
    %transform all input to sampling points
    if length(cfg.times)==3
        cfg.sampRate=cfg.times(3);
        timesTemp=cfg.times(1):1/cfg.sampRate:cfg.times(2);
    elseif length(cfg.times)==2
        if ~isfield(cfg,'sampRate')
            error('Sampling rate must be specified as third value in cfg.times or in cfg.sampRate')
        end
        timesTemp=cfg.times(1):1/cfg.sampRate:cfg.times(2);
    elseif length(cfg.times)==1
        error('cfg.times must include at least two values (start and end times).')
    elseif ~isfield(cfg,'sampRate') % and if all time points are indicated
        cfg.sampRate=1/mean(diff(cfg.times));
        if mod(cfg.sampRate,1)
            cfg.sampRate=round(cfg.sampRate);
            warning('No sampling rate indicated; cfg.sampRate set to %g',cfg.sampRate);
        end
    end
    if length(cfg.times)<4
        if timesTemp(1) ~= cfg.times(1) || abs(timesTemp(end)-cfg.times(2)) > cfg.sampRate
            error('Indicated sampling rate does not match indicated times, please correct cfg.sampRate or cfg.times')
        end        
        cfg.times=timesTemp;
    else
        if ~all(diff(cfg.times)-mean(diff(cfg.times))<1e6) %deviations in sub-microsecond range do not matter
            error('Apparently non-constant sampling rate - check cfg.times!')
        end
    end
    if length(cfg.times) ~= size(avgs,3)
        error('Something went wrong, please check cfg.times and (if in use)) cfg.sampRate (all correct, all in same unit [ms or s]?)')
    end
    
    
    %transform to sampling points
    cfg.peakWin=nearestN(cfg.times,cfg.peakWin);
    if isfield(cfg,'cWin') && isnumeric(cfg.cWin)
        cfg.cWin=nearestN(cfg.times,cfg.cWin);
    end
    if isfield(cfg,'cWinWidth')
        cfg.cWinWidth=round(cfg.cWinWidth*cfg.sampRate);
    end
    if isfield(cfg,'areaWin') && isnumeric(cfg.areaWin)
        cfg.areaWin=nearestN(cfg.times,cfg.areaWin);
    end
    if isfield(cfg,'ampLatWin') && isnumeric(cfg.ampLatWin)
        cfg.ampLatWin=nearestN(cfg.times,cfg.ampLatWin);
    end
else
    %check whether all timing information is indicated in sampling points
    msg='if cfg.times is not set, all timing information is interpreted as sampling points and must consequently be in integers. Please check: ';
    throwError=false;
    if isfield(cfg,'peakWidth') && mod(cfg.peakWidth,1)
        msg=[msg,' peakWdith'];
        throwError=true;
    end
    if isfield(cfg,'peakWin') && any(mod(cfg.peakWin,1))
        msg=[msg,' peakWin'];
        throwError=true;
    end
    if isfield(cfg,'cWinWidth') && mod(cfg.cWinWidth,1)
        msg=[msg,' cWinWidth'];
        throwError=true;
    end
    if isfield(cfg,'meanTime') && any(mod(cfg.meanTime,1))
        msg=[msg,' meanTime'];
        throwError=true;
    end
    if throwError
        error(msg)
    end
end
if isfield(cfg,'peakWidth')
    if isfield(cfg,'times')
        cfg.peakWidth=round(cfg.peakWidth*cfg.sampRate);
    end
else
    if isfield(cfg,'times') && cfg.sampRate >= 50 %likely in samples/s (Hz)
        peakWidthS=0.005;
        cfg.peakWidth=round(peakWidthS*cfg.sampRate);
        peakWidthS=cfg.peakWidth/cfg.sampRate;
        warning('No peak width (cfg.peakWidth) indicated; peakWidth set to %g s.',peakWidthS)              
    elseif isfield(cfg,'times') %if sampRate < 50, it is likely in samples/ms
        peakWidthMs=5;
        cfg.peakWidth=round(peakWidthMs*cfg.sampRate);
        peakWidthMs=cfg.peakWidth/cfg.sampRate;
        warning('No peak width (cfg.peakWidth) indicated; peakWidth set to %g ms.',peakWidthMs) 
    else %in sampling points
        cfg.peakWidth=5;
        if cfg.warnings
            warning('No peak width (cfg.peakWidth) indicated; peakWidth set to %g asmpling points.',cfg.peakWidth)
        end
    end
end
if ~isfield(cfg,'times')
    cfg.times=1:size(avgs,3);
    cfg.sampRate=1;
end
if ~isfield(cfg,'extract')
    cfg.extract={'peakLat','onset','offset','areaLat','mean','peakAmp','area','width','peak2peak','baseline'};
    if cfg.warnings
        warning('No output measures chosen; will extract all possible output measures:')
        fprintf('%s, ',cfg.extract{:});
        fprintf('\n');
    end
end
if strcmp(cfg.extract,'all')
    cfg.extract={'peakLat','onset','offset','areaLat','mean','peakAmp','area','width','peak2peak','baseline'};
end

if ~isfield(cfg,'meanTime'),cfg.meanTime=cfg.peakWin;end
cfg.get.peak=true; %needs nothing in addition

if ismember('areaLat',cfg.extract)
    if ~isfield(cfg,'areaWin'),cfg.areaWin='peakWin';end
    if ~isfield(cfg,'percArea'),cfg.percArea=.5;end
    if ~isfield(cfg,'areaBase'),cfg.areaBase='zero';end
    cfg.get.areaLat=true;
else
    cfg.get.areaLat=false;
end

if any(ismember({'onset','offset','width'},cfg.extract)) || (cfg.get.areaLat && strcmp(cfg.areaWin,'ampLat'))
    cfg.get.ampLat=true; %needs only peak
    if ~isfield(cfg,'ampLatWin'),cfg.ampLatWin='fullRange';end
else
    cfg.get.ampLat=false;
end

if any(ismember({'areaLat','area'},cfg.extract))
    cfg.get.area=true;
else
    cfg.get.area=false;
end
%searches for the point where the set percentage  of the total area is
%reached and therefore critically depends on area
if ~isfield(cfg,'fig'),cfg.fig=false;end
if ~isfield(cfg,'aggregate'),cfg.aggregate='individual';  end

if isfield(cfg,'cWinWidth')
    if ~isfield(cfg,'cWinStart')
        cfg.cWinStart='peakWin';
    end
    %if set, percAmp determines the peak-to-peak threshold instead of the
    %x-axis-to-peak threshold
    if ~isfield(cfg,'cBound'),cfg.cBound=true;end
    %search window for releative latency is bounded on one side by the counter peak
    if cfg.cBound,cfg.get.counterPeak=1;end
elseif isfield(cfg,'cWin')
    cfg.get.counterPeak=1;
else
    cfg.get.counterPeak=0;
end

if isfield(cfg,'percAmp')
    cfg.get.newBaseline=1;
    %the area is bound in amplitude by percAmp of the peak amplitude or of
    %the peak-to-peak amplitude (when counterPeak is calculated)
else
    %otherwise area uses the typical baseline
    cfg.get.newBaseline=0;
    if any(ismember({'onset','offset','width'},cfg.extract))
        cfg.percAmp=0.5; %because 0 does not make sense (will be caught by noise)
    else
        cfg.percAmp=0; %makes sense for areaLat
    end
end


if cfg.sign==1
    direction='posi';
elseif cfg.sign==-1
    direction='nega';
end
%annotations that explain in some detail what each measure reflects
cfg.ann.mean=sprintf('Mean activity in the interval %g-%g',...
    cfg.meanTime(1),cfg.meanTime(2));
cfg.ann.peak=sprintf('Most %stive peak in the search interval %g-%g',...
    direction,cfg.peakWin(1),cfg.peakWin(2));
if cfg.get.counterPeak
    if isfield(cfg,'cWin')
        cfg.ann.counterPeak=sprintf('Strongest peak of opposite polarity in the time window %g - %g',cfg.cWin);
    else
        if cfg.cWinWidth<0 && strcmp(cfg.cWinStart,'peak')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polarity in between %g before the actual peak and the actual peak',-cfg.cWinWidth);
        elseif cfg.cWinWidth>0 && strcmp(cfg.cWinStart,'peak')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polartity in between the actual peak and %g after the actual peak',cfg.cWinWidth);
        elseif cfg.cWinWidth<0 && strcmp(cfg.cWinStart,'peakWin')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polarity in between %g before the start of the search window and the actual peak',-cfg.cWinWidth);
        elseif cfg.cWinWidth>0 && strcmp(cfg.cWinStart,'peakWin')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polartity in between the actual peak and %g after the end of the search window',cfg.cWinWidth);
        end
    end
end
if cfg.get.newBaseline
    if cfg.get.counterPeak
        cfg.ann.baseline=sprintf('Activity %g%% of the peak-to-counter-peak distance away from the peak\n(lower values are closer to the peak)',...
            cfg.percAmp*100);
    else
        cfg.ann.baseline=sprintf('Activity %g%% of the peak amplitude away from the peak amplitude\n(lower values are closer to the peak)',...
            cfg.percAmp*100);
    end
    cfg.ann.onset=sprintf('Point before the peak, where activity crosses the baseline\n(see cfg.ann.baseline)');
    cfg.ann.offset=sprintf('Point after the peak, where activity crosses the baseline\n(see cfg.ann.baseline)');
    
    cfg.ann.area=sprintf('Cumulative sum over all %stive samples in between on- and offset (see cfg.ann.onset),\nthe ERP and the %g%% threshold (see cfg.ann.baseline)',...
        direction,cfg.percAmp*100);
else
    cfg.ann.area=sprintf('Cumulative sum over all %stive samples in between on- and offset (see cfg.ann.onset),\nthe ERP and the pre-stimulus baseline activity',...
        direction);
end
if cfg.get.areaLat
    cfg.ann.areaLat=sprintf('Point in time where %g%% of the total area (see cfg.ann.area) is reached',cfg.percArea*100);
end
end

function samplesOut=nearestN(allTimes,timesIn) %search for nearest neighbor
% latency and time in same units
% allTimes   - array of equidistant time values (ideally)
% timesIn    - sample indices 
[samplesOut,delta] = deal(zeros(size(timesIn)));
maxDelta      = max(diff(allTimes)); %because it is not always perfectly equidistant
for i=1:length(timesIn)
    [delta(i),samplesOut(i)] = min(abs(allTimes-timesIn(i)));
    if delta(i) > maxDelta
        error('Time point %g not found in time range [%g, %g] or something is wrong about cfg.times',...
            timesIn(i),allTimes(1),allTimes(end));
    end
end
end

function peak=peakDetection(cfg,ERP,peakLat)
%will search for counterpeak if peakLat is set

if nargin<3
    peakWin=cfg.peakWin;
    direction=cfg.sign;
else
    %counter peak has the opposite polarity
    direction = cfg.sign*-1;
    if length(peakLat)==1
        switch cfg.cWinStart
            case 'peak'
                if cfg.cWinWidth < 0 %preceding counter peak
                    peakWin=[peakLat+cfg.cWinWidth,peakLat];
                elseif cfg.cWinWidth >= 0 %following counter peak
                    peakWin=[peakLat,peakLat+cfg.cWinWidth];
                end
            case 'peakWin'
                if cfg.cWinWidth < 0 %preceding counter peak
                    peakWin=[cfg.peakWin(1)+cfg.cWinWidth,peakLat];
                elseif cfg.cWinWidth >= 0 %following counter peak
                    peakWin=[peakLat,cfg.peakWin(2)+cfg.cWinWidth];
                end
        end
    elseif length(peakLat)==2
        peakWin=peakLat;
    end
end

%make sure not to get a peak that lies directly at one of the ends of the
%ERP so that amplitude can always be calculated
peakWin(1)=max(peakWin(1),cfg.peakWidth+1);
peakWin(2)=min(peakWin(2),length(ERP)-cfg.peakWidth);
    
ERPcut = ERP(peakWin(1):peakWin(2));

if direction==-1
    ERPcut=ERPcut*-1;
    %if a negative peak is searched, search for local positive peak on
    %the reversed  polarity ERP
end
if length(ERPcut)>2
    [dummy,lats] = findpeaks(ERPcut,'SORTSTR', 'descend');
    if numel(lats)<1
        %if no local peak is found (this is what findpeaks does), simply take the maximum
        [dummy,lat]=max(ERPcut);
        %track whether a local peak was found or not
        peak.foundLocal= false;
    else %take the strongest local peak
        lat = lats(1);
        peak.foundLocal = true;
    end
else
    if cfg.warnings
        warning('Search range for peak is too narrow (less than 3 samples)!')
    end
    [dummy,lat]=max(ERPcut);
    peak.foundLocal= false;
end

%convert from relative to search range TO relative to episode onset
peak.lat = peakWin(1)+lat-1;

peak.amp = mean(ERP(peak.lat-cfg.peakWidth:peak.lat+cfg.peakWidth));
%use original data instead of cutERP, because the peakWin might overlap
%with start/end of cutERP
end

function ampLat=amplitudeLatency(cfg,ERP,baseline,peakLat,counterLat)
%ampLat.start: where ERP reaches cfg.percAmp before the peak
%ampLat.end:   where ERP reaches cfg.percAmp after the peak


%determine the minimal and maximal values for ampLat
%in case of peakWin, one side is bounded by the counter peak (if there is one)
if ischar(cfg.ampLatWin)
    switch cfg.ampLatWin
        case 'peakWin'
            time=cfg.peakWin;
        case 'fullRange'
            time=[1,length(ERP)];
        otherwise
            error('cfg.ampLatWin must either be "peakWin", "fullRange" or start and end time indices; it is %s',cfg.ampLatWin)            
    end
    if exist('counterLat','var') && cfg.cBound
        %counter peak is one border of the interval
        if counterLat<peakLat
            time(1)=counterLat;
        elseif counterLat>=peakLat
            time(2)=counterLat;
        end
        
    end
elseif isnumeric(cfg.ampLatWin) && numel(cfg.ampLatWin)==2;
    time=cfg.ampLatWin;
else
    error('cfg.ampLatWin must either be "peakWin", "fullRange" or start and end time indices')
end

%make sure to avoid samples that lie directly at one of the ends of the
%ERP so that we certainly can calculate an amplitude
time(1)=max(time(1),cfg.peakWidth+1);
time(2)=min(time(2),length(ERP)-cfg.peakWidth);


%first search for early crossing of baseline by going backwards from
%peak
ampLat=struct;
for sample=peakLat:-1:time(1)+cfg.peakWidth
    %go through the search window and calculate mean area (peakWidth
    %broad) and each time check whether baseline has been reached
    currentAmp = mean(ERP(sample-cfg.peakWidth:sample+cfg.peakWidth));
    if (cfg.sign==1 && currentAmp<=baseline) || (cfg.sign==-1 && currentAmp>=baseline)
        ampLat.start=sample;
        break;
    end
    
end

%now search for the later crossing of the threshold
for sample=peakLat:time(2)-cfg.peakWidth
    currentAmp = mean(ERP(sample-cfg.peakWidth:sample+cfg.peakWidth));
    if (cfg.sign==1 && currentAmp<=baseline) || (cfg.sign==-1 && currentAmp>=baseline)
        ampLat.end=sample;
        break;
    end
end

%When ERP did not cross the baseline in the search interval,
%percent-amplitude latency is set to the borders of the interval
if isfield(ampLat,'start')
    ampLat.foundOn=true;
else
    ampLat.foundOn=false;
    ampLat.start=time(1);
end
if isfield(ampLat,'end')
    ampLat.foundOff=true;
else
    ampLat.foundOff=false;
    ampLat.end=time(2);
end

end %end of sub-function ampLat

function area=totalArea(cfg,ERP,baseline,time)
cutERP = (ERP(time(1):time(2))-baseline)*cfg.sign;
area = sum(cutERP(cutERP>0))*cfg.sign; % area is negative if the component is negative
end %end of sub-function totalArea

function areaLat=areaLatency(cfg,ERP,area,baseline,time)
%time in samples (taken care of by checkConfig)
cutERP=(ERP(time(1):time(2))-baseline)*cfg.sign; % ERP & area need to be flipped if they are negative
latidx = find(cumsum(cutERP) >= area * cfg.sign * cfg.percArea);
if ~isempty(latidx)
   areaLat.lat = latidx(1) + time(1) -1;
   areaLat.foundLat=true;
else
    areaLat.lat = NaN;
    areaLat.foundLat=false;
end
end %end of sub-function areaLat



