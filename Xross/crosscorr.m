function [pivlayers, pivmeanvalid] = crosscorr(inkymo,slit, pixelsize, frameinterval)

%INPUTS: (matrix), (scalar), (scalar), (vector)
%%%% filters %%%%
precorrImg_threshold = 1; % std units above mean
pairwisecorr_threshold = 0 ; % std units above mean
fullcorr_threshold = 1 ; % std units above mean
conewidthcriterium = .5 ; % of max

bckfcorrflag = 0; % for now, use always 0
%%%%%%    %%%%%%%

kym3 = kym3dim(inkymo,slit);

s = size(kym3) ;
nframes = s(3) ;
framesize = size (kym3(:,:,1)) ;
pixelvelocity = 60*pixelsize/frameinterval ;  % calibration: velocity (in micron/min !!)
%  which corresponds to moving 1 pixel in 1 frame

prompt = { sprintf(['Kymo dimensions: ', num2str(s(1)), ' (V) x ', num2str(s(2)), ' (H) x ', num2str(s(3)),' (T) [pixels]', ...
    '\n                                ', ...
    num2str(pixelsize*s(1)), ' (V) x ', num2str(pixelsize*s(2)), ' (H) x ', num2str(frameinterval*s(3)),' (T) [um x um x sec] \n\n', ...
    ' 1 pixel shift per frame would be ', num2str(pixelvelocity), ' um/min\n\nMax. velocity is (um/min):']),'Starting frame : ','Stack size' };
def = {'25','1',num2str(nframes)};
dlgTitle = '';
dialog = inputdlg(prompt,dlgTitle,2,def);
maxvelocity = str2num(dialog{1})
startingframe=str2num(dialog{2})
blocksize=str2num(dialog{3})

semiwindow = ceil(maxvelocity/pixelvelocity) % WINDOW TO SEARCH NEIGHBORHOOD in pixels!

% Define the kernel as the next image cropped to allow a valid (weighted)
% region of 'semiwindow' neighbourhood
beginy = round( framesize(1)/2-semiwindow );
endy =   round( framesize(1)/2+semiwindow );
beginx = round( framesize(2)/2-semiwindow );
endx =   round( framesize(2)/2+semiwindow );

kym3crop = kym3( beginy:endy , beginx:endx , startingframe:(startingframe+blocksize-1) ) ;
nframes = blocksize;
s = size(kym3) ;


%%% uncorrelated background correlation
box1=kym3crop(:,:,round(nframes*3/4));
box2=kym3(:,:,round(nframes*1/4));
bckgcorr  = normxcorr2(box1,box2) ;
%%% END OF uncorrelated background correlation


% Correlation between immediate neighbouring frames
for i = 1:nframes-1
    
    box1 = kym3crop(:,:,i+1);
    box2 = kym3(:,:,i);
    boxhelp = kym3crop(:,:,i);
    
    %% PRE: Constraints at the single-layer level
    level = mean(box1(:)) + precorrImg_threshold * std(box1(:));
    tempmask = box1>level;    ;
    box1 = box1 .* tempmask ;
    tempmask = box2>level;
    box2 = box2 .* tempmask ;
    
    %% Covariance (or Correlation)
    layers(:,:,i)  = normxcorr2(box1,box2) ;
    size(layers(:,:,i)), size(bckgcorr)
    layers(:,:,i) = fliplr(flipud(layers(:,:,i) - bckfcorrflag*bckgcorr));
    s=size(layers(:,:,i)) ;
    layersraw(:,:,i) = layers(:,:,i) ;
    
    %% VALID SHIFTS CROP
    temp = layers(:,:,i);
    center = round(size(temp)/2);
    pivmeanvalid(:,:,i) = temp (center(1)-semiwindow : center(1)+semiwindow , ...
        center(2)-semiwindow : center(2)+semiwindow );
    pivlayers(:,:,i) = pivmeanvalid(:,:,i);
    %% POST: Constraints at the single-layer level after Correlation
    temp = pivlayers(:,:,i);
    
    % CASE 1: x-sigma level threshold
    level = mean(temp(:)) + pairwisecorr_threshold * std(temp(:));
    tempmask = temp>level;    
    mean(temp(:));
    pivlayers(:,:,i) = temp.*tempmask ./ max(temp(:));
    
      
    h = figure(10), clear gcf, subplot(ceil(sqrt(nframes)),ceil(sqrt(nframes)),i), imagesc(pivlayers(:,:,i)); axis off, axis equal,axis tight, hold on, plot(round([1 s(2)]), round([s(1)/2 s(1)/2])) ;
    set(h,'Visib','off') ;
    
end

pivmeanvalid = sum(pivlayers,3);
pivmeanvalid = pivmeanvalid(end:-1:1,:);

% Mask preps
[aa,bb] = size(pivmeanvalid) ;
[X, Y] = meshgrid( (1:aa)-(aa+1)/2, (1:bb)-(bb+1)/2) ;

R = sqrt (X.^2 + Y.^2) ;
theta = atan2d (Y,X) ;
%%%%%%%%%%%%%%%%%%

% Mask: remove correlation map corners (because vel. > max vel. designated (inputdlg))
pivmeanvalid = pivmeanvalid .* (R < floor(min(aa,bb)/2)) ;

% Mask: remove velocities at an angle higher than designated (inputdlg)
pivmeanvalid = pivmeanvalid .* ( abs(theta)<=181  );

% Subtract background: estimated as the average of the central vertical
%line (with margin 1+2*semiwid), corresponding to particles moving at +-90º, which are assumed to be
%noise.
semiwid=2; %plus the central element
bckgwindow = pivmeanvalid( : , ((bb+1)/2-semiwid):((bb+1)/2+semiwid));
%pivmeanvalid = pivmeanvalid - mean(bckgwindow(:))  ;

% Mask2: norm threshold below which values are ignored
sigmalevel = mean(pivmeanvalid(:)) + fullcorr_threshold*std(pivmeanvalid(:));
mask2 = pivmeanvalid > sigmalevel; 


% Plot pivmeanvalid

figure(2)
helpcrosscorr=uicontrol('Style','text','String',sprintf('HELP\nhover with mouse'),'Position',[.9 .5 .1 .1],'Un','Norm')
set(helpcrosscorr,'Tooltip',sprintf('Map on the right is a filtered version ogf that on the left.\n Filtering criteria are defined in the beginning of the script (LAPSOcrosscorr.m; section ''filters'').\n\n On the map on the right the centroid of the map is calculated, defining the VelocityWhite and the FlowAngle_white.\nDispersion angle is determined by intersection of the polar curve with 1/2 radial coordinate around the centroid angle.\n\nWith the dispersion angle at hand, a refined value for flow velocity and angle is retrieved by \n calculating a new centroid based on a correlation map which is truncated at the 1/2 criterium edge angles.  '))
subplot (2,2,1), hold on
imagesc( linspace(-maxvelocity, +maxvelocity, 2*semiwindow+1), ...
    linspace(-maxvelocity, +maxvelocity, 2*semiwindow+1), ...
    + interp2(pivmeanvalid,0) ) ,
plot([-maxvelocity +maxvelocity],[0 0],'w-','LineW',1)
plot([0 0],[-maxvelocity +maxvelocity],'w-','LineW',1)
impixelinfo, axis equal, axis tight
ylabel('y velocity (um/min)')
xlabel('x velocity (um/min)')

s=size(pivmeanvalid);
[ymax xmax] = find( pivmeanvalid == max(pivmeanvalid(:)) );
c1 = ymax*pixelvelocity-maxvelocity;
c2 = xmax*pixelvelocity-maxvelocity;


% Plot {pivmeanvalid.*mask2} (i.e. with after-math threshold)
PIVFINAL = pivmeanvalid .* mask2;

subplot(2,2,2), hold on
imagesc(linspace(-maxvelocity, +maxvelocity, 2*semiwindow+1), ...
    linspace(-maxvelocity, +maxvelocity, 2*semiwindow+1), ...
    PIVFINAL)  % 

plot([-maxvelocity +maxvelocity],[0 0],'w-','LineW',1)
plot([0 0],[-maxvelocity +maxvelocity],'w-','LineW',1)
impixelinfo, axis equal, axis tight
ylabel('y velocity (um/min)')
xlabel('x velocity (um/min)')

%%%%%%%%%%%%%%
% RADON-BASED AZIMUTHAL SCANNING

for nn=[3 4]
    if nn==4
        temppiv=pivmeanvalid;
        pivmeanvalid = PIVFINAL ;
    end
    
    pivmeanvalid = pivmeanvalid(:,end:-1:1);
    
    sizecorr = size(pivmeanvalid,1); % Correlation map size
    pivmeanvalid((sizecorr+1)/2,(sizecorr+1)/2) = 0;
    
    
    topcorr = pivmeanvalid ;
    mask= theta>=-180 & theta<=-0; topcorr = topcorr.*abs(1-mask);
    RRtop = radon( topcorr , (45:134)-90 ) ;
    
    leftcorr = pivmeanvalid;
    mask= theta>=-90 & theta<=90; leftcorr = leftcorr.*abs(1-mask);
    RRleft = radon( leftcorr , (135:224)-90 ) ;
    
    bottomcorr = pivmeanvalid;
    mask= theta>=0 & theta<=180; bottomcorr = bottomcorr.*abs(1-mask);
    RRbottom = radon( bottomcorr, (225:315)-90 ) ;
    
    rightcorr = pivmeanvalid;
    mask= theta>=90 | theta<=-90; rightcorr = rightcorr.*abs(1-mask);
    RRright = radon( rightcorr , (-45:44)-90 ) ;
    
    sizeradon = size(RRtop); % Radon map size
    
    RRtopp = RRtop( round(sizeradon(1)/2+0.5) , : ) ;
    RRrightt = RRright( round(sizeradon(1)/2+0.5) , : ) ;
    RRbottomm = RRbottom( round(sizeradon(1)/2+0.5), : ) ;
    RRleftt = RRleft( round(sizeradon(1)/2+0.5) , : ) ;
    
    
    figure,    
    imagesc([topcorr leftcorr bottomcorr rightcorr]), axis equal, grid on
    
    figure(2)
    subplot(2,2,nn) ,
    
    % rotscan is the full (2pi) radial (along radius) projection of the correlation map
    rotscan = [RRtopp RRrightt  RRbottomm RRleftt ];
    rotscan = circshift( rotscan,[1,45] );
    
    rotscanNorm = rotscan / max(rotscan(:)) ;    
    polar( 0: pi/180: 2*pi , rotscanNorm ) ;    
    figure(3), subplot(1,2,nn-2), plot( 180/pi *(0: pi/180: 2*pi) , rotscanNorm  ) ;
    
    if nn==4
        pivmeanvalid = temppiv ;
        
    end
    
end
%END OF RADON %%%%%%%%


% CALCULATIONS AND CLIPBOARD
% Calculate center of mass on pivfinal
figure, imagesc(PIVFINAL)
xx = 1 : size(PIVFINAL,2) ;
yy = 1 : size(PIVFINAL,1) ;

xproj = sum( PIVFINAL,1 ) ;
xcm = sum(xproj.*xx) ./ sum(xproj);
v_x = -maxvelocity + (xcm-1)*pixelvelocity;

yproj = sum( PIVFINAL,2 );
ycm = sum(yproj'.*yy) ./ sum(yproj);
v_y = -maxvelocity + (ycm-1)*pixelvelocity;

figure(2)
subplot(2,2,2)
plot( v_x , v_y , 'w+','MarkerSize',8)

plot( v_x , v_y , 'w+','MarkerSize',8)
vxy = sqrt( v_x^2 + v_y^2 ) ;
thetaV = atan2d( v_y , v_x ) ;

% ANGULAR RANGE DETERMINATION

% Option 1
% Find crossings of rotscan through 50% ('~FWHM);
% ... susceptible to multipeak interference
insidecone = find(rotscanNorm > conewidthcriterium)

if insidecone(1) > 1
    semicone = abs( (insidecone(1)-insidecone(end)) / 2 )
else
    findstep= find(max(diff(insidecone))) ;
    semicone = abs((insidecone(findstep) - insidecone(findstep+1)) / 2);
end

% Option 2

rotscanNormtripled = [rotscanNorm rotscanNorm rotscanNorm] ;

rotscanNormtripled( mod(thetaV,360)+360 : +1 : mod(thetaV,360)+360+90 )
aftercone = find( rotscanNormtripled( mod(thetaV,360)+360 : +1 : mod(thetaV,360)+360+90 )  < conewidthcriterium );
aftercone = aftercone(1)
beforecone = find( rotscanNormtripled( mod(thetaV,360)+360  : -1 : mod(thetaV,360)+360-90 ) < conewidthcriterium );
beforecone = beforecone(1)

maskCone = theta>thetaV-beforecone & theta<thetaV+aftercone;
PIVFINALcone = PIVFINAL .* maskCone;

xx = 1 : size(PIVFINALcone,2) ; 
yy = 1 : size(PIVFINALcone,1) ; 

xproj = sum( PIVFINALcone,1 ) ; 
xcm2 = sum(xproj.*xx) ./ sum(xproj) ; 
v_x2 = -maxvelocity + (xcm2-1) * pixelvelocity ; 

yproj = sum( PIVFINALcone,2 ) ; 
ycm2 = sum(yproj'.*yy) ./ sum(yproj) ; 
v_y2 = -maxvelocity + (ycm2-1) * pixelvelocity ; 

figure(2) 
subplot(2,2,2) 
plot( v_x2 , v_y2 , 'k+','MarkerSize',12) 

vxy_black = sqrt( v_y2^2+v_x2^2 ) ; %constrained to define center of mass within the cone between aftercone and beforecone
thetaVblack = atan2d( v_y2 , v_x2 ) ;

annotation('textbox', [0.0, 0.95, 1, 0.05], 'String', sprintf(['|v white| = ' , num2str(vxy),'  um/min;    Flow Angle ', num2str(thetaV),'º  +/- ', num2str(0.5*(-1+beforecone+aftercone)),'º     ;  (1 pixel is ',num2str(pixelvelocity),' um/min )']),'FontS',10,'EdgeColor','none')
annotation('textbox', [0.0, 0.92, 1, 0.05], 'String', sprintf([ '|v black| = ' , num2str(vxy_black),'  um/min;    Flow Angle ', num2str(thetaVblack) ] ),'EdgeColor','none') ;

annotation('textbox', [0.0, 0.05, 1, 0.05], 'String', ' values to clipboard are [ avg.velocity ,  flow direction  ,   +/- flow direction  ' ,'FontS',8 )

toclipboard( [vxy_black  thetaVblack  0.5*(-1+beforecone+aftercone)] ) 

set(h,'Visible','on')

end


