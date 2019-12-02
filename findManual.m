function [blinks, goBack] = findManual(data, manualCheck,nrSegments, escapeKey, ploti)

% e.g. findManual(rebinnedData, manualCheck, trialDur/EEG.srate)

if manualCheck == 1
blinks = [];
countf = 0;
countb = 0;
segmenti = 0;

while (segmenti < nrSegments)
segmenti = segmenti +1;
% for segmenti = 1:nrSegments

subplot(2,1,ploti)
if size(data,1) < 500 
    plot(data(:,:,segmenti)') 
else
    plot(data(:,segmenti))    
end

keyIsDown = 0;  
    while keyIsDown == 0
        title('***')
        pause
        [ keyIsDown, timeSecs, keyCode ] = KbCheck;
       if keyIsDown
        KeyName = keyCode;
          keyIsDown = keyIsDown +1;
            if keyCode(escapeKey)
                break;
            end
       end                                                                                                                                               % jääb nupuvajutust ootama
    end
    
% blinks{subi, eventIndx, thisCueEvent, segmenti} = 1;      
         
    if KeyName(40) == 1
        blinks(segmenti,1) = 0; % down arrow
        title('leaving out')
        pause(0.2)
    elseif KeyName(38) == 1 % up arrow
        blinks(segmenti,1) = 1;
        title('leaving in')
        pause(0.2)
    else
        title('wrong key!')
    end
    
    keyIsDown = 0;  
    while keyIsDown == 0
        pause
        [ keyIsDown, timeSecs, keyCode ] = KbCheck;
       if keyIsDown
        KeyName = keyCode;
          keyIsDown = keyIsDown +1;
            if keyCode(escapeKey)
                break;
            end
       end                                                                                                                                               % jääb nupuvajutust ootama
    end
    
        if KeyName(37) == 1 % left arrow
            title('going back')
            goBack = -2 * countb;
            countb = countb + 1;
        elseif KeyName(32) == 1 % space
            title('going back')
            goBack = -4 * countb;
            countb = countb + 1;
        elseif KeyName(39) == 1 % right arrow
%             title('going forward')
%             goBack = 2 * countf;
%             countf = countf + 1;
              blinks = ones(nrSegments,1);
              segmenti = nrSegments;
              goBack=0;
        else
            goBack = 0; 
        end
       
        
end
else
    
% blinks = ones(nrSegments,1);

for ji = 1:nrSegments
blinks(ji,1) = randi(0:1);
end
goBack = 0;

end
hold off
end