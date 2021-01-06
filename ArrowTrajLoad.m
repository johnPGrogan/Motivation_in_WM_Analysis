function traj = ArrowTrajLoad(dataFiles)
% load up trajectories from Arrow Task
% look at some trials
% verify angles are correct (inverted?)
% look at initial traj angle instead of final angle


%% load traj and stim info
nPP = length(dataFiles);

for iPP = 1:nPP
    d = load( dataFiles{iPP} );
    result = d.result;
    
    %% load up the things that are fixed across trials
    
    traj(iPP).centre = result.params.screenSize / 2; % get screen centre
    
    if isfield(result.params,'maxMouseRadius')
        traj(iPP).maxMouseRad = result.params.maxMouseRadius; % max radius for mouse resps
        traj(iPP).mouseBound = [sin(0:.1:2*pi); cos(0:.1:2*pi)]' .* traj(iPP).maxMouseRad; % make into circle
    end

    traj(iPP).colours = result.params.targetColours ./255; % colours in order
    
    traj(iPP).targetShape = result.params.targetShape; % shape of arrow
    traj(iPP).targetShape(end+1,:) = traj(iPP).targetShape(1,:);%complete the shape
    
    
    %get name of probe start time
    if isfield(result.data(1), 'probeStart')
        probeStart = 'probeStart';
    elseif isfield(result.data(1),'startProbe')
        probeStart = 'startProbe';
    else
        error('probeStart or startProbe field not found for probe onset time')
    end
    %% now load up things that vary on each trial
    
    nTrials = length(result.data);
    
    for iTrial = 1:nTrials
        traj(iPP).trials(iTrial).nTargets = result.data(iTrial).nTargets;%num targets
        traj(iPP).trials(iTrial).targInd = result.data(iTrial).targetIndex;
        
        %get stim locations
        stimLoc = result.data(iTrial).locations; % coords of stim
        traj(iPP).trials(iTrial).stimLoc = stimLoc - traj(iPP).centre; % centre them
        %also get these as complex coords
        traj(iPP).trials(iTrial).stimLocCompl = real(traj(iPP).trials(iTrial).stimLoc(:,1)) + real(traj(iPP).trials(iTrial).stimLoc(:,2))*1j;
        
        %get stim angles and colours
        traj(iPP).trials(iTrial).stimAngles = mod(result.data(iTrial).angles,2*pi);%angles of stim
        traj(iPP).trials(iTrial).stimColours = result.data(iTrial).colourOrder;%colour inds for trial
        
        %get trajectory
        traj(iPP).trials(iTrial).traj = result.data(iTrial).trajectory(:,2); % trajectory of responses
        traj(iPP).trials(iTrial).trajTime = real(result.data(iTrial).trajectory(:,1)); % time of each mouse sample
        
        %startProbe
        traj(iPP).trials(iTrial).trajTime = traj(iPP).trials(iTrial).trajTime - result.data(iTrial).(probeStart); % time relative to probe onset
        traj(iPP).trials(iTrial).trajTimeRel = traj(iPP).trials(iTrial).trajTime - traj(iPP).trials(iTrial).trajTime(1); % time relative to movement onset
        
        %get angle from centre of screen at each time point
        traj(iPP).trials(iTrial).angles = angle(traj(iPP).trials(iTrial).traj); % get angles of resps
        
        %get angle of first movement
        traj(iPP).trials(iTrial).initAngle = mod(traj(iPP).trials(iTrial).angles(1),2*pi);
        %get final angle
        traj(iPP).trials(iTrial).finalAngle = mod(traj(iPP).trials(iTrial).angles(end),2*pi);
        
        %place angles onto radius circle - need to invert angle here and
        %shift by 90deg
        if isfield(traj(iPP),'maxMouseRad')
            traj(iPP).trials(iTrial).angleCircle = [sin(pi/2-traj(iPP).trials(iTrial).angles),cos(pi/2-traj(iPP).trials(iTrial).angles)] .* traj(iPP).maxMouseRad; % place angles around a circle
        end
        
        traj(iPP).trials(iTrial).finalPrec = (mod(traj(iPP).trials(iTrial).finalAngle - traj(iPP).trials(iTrial).stimAngles(traj(iPP).trials(iTrial).targInd) + pi, pi*2) - pi); % precision
        traj(iPP).trials(iTrial).initPrec = (mod(traj(iPP).trials(iTrial).initAngle - traj(iPP).trials(iTrial).stimAngles(traj(iPP).trials(iTrial).targInd) + pi, pi*2) - pi); % precision
        
        if isfield(result.data(iTrial),'showProbeInitially')
            traj(iPP).trials(iTrial).arrowProbe = result.data(iTrial).showProbeInitially;%whether arrow probe shown
        else
            traj(iPP).trials(iTrial).arrowProbe = [];%whether arrow probe shown
        end
        traj(iPP).trials(iTrial).probeAngle = result.data(iTrial).targetAngle + result.data(iTrial).responseOffset;
    end
    
    
end