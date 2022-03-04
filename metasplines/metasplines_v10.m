function [fo, vo, ROICROP, ROICROPdummy, outdata] = metasplines_v10(background, threeDrepresentationsFlag)

% DESCRIPTION_
% Inter-polar versus kinetochore-attached microtubules quantification in a
% mitotic spindle. Also generates patch visualization of the microtubule
% density and of the generated spline plates.

% Metasplines_v10 reads an excel sheet containing kinetochore 3D coordinates
% previously obtained (e.g., in ImageJ).
% Kinetchore plates (spline surfaces connecting kinetochores)
% are calculated, dividing the cell volume in 3 regions:
% 'outer 1', 'outer 2', and 'intra-plates'.
%
% Based on parameters below ('polshift' and 'thickness'),
% one quantification volume is defined in each of the 3 regions:
%
% - Centromeric plate: this is a centered sub-slice in between the two kinetochores plates;
% this thin and curved 3D volume travels along the centromeric regions of the cell.
%
% - Two Poleward plates: each is a thin poleward-shifted replica of the kinetochore
% plate of the corresponding side. Slices thickness and poleward shift are
% user-defined. To strengthen the assumptions below, these plates should be defined
% very close (e.g. 200nm offset) to the kinetochore plates.
%
% Background-subtracted fluorescence signal is integrated in each of these volumes.
% Assuming that:
% i) only inter-polar microtubules cross the centromeric plate and
% ii) all the above inter-polar microtubules cross also the poleward plates
% (which are also crossed by the kinetochore-attached microtubules,
% then the inter-polar to kinetochore-attached microtubules ratio can be
% calculated. It is output in 'outdata.ipMTproportion'.
%
% INPUTS______
% background: previously estimated photon count, typically in a unpopulated
% sub-region of the 'intra-plates' region. For the quantitative part, it is
% very important to have a good estimate for background.
%
% threeDrepresentationsFlag (0-1): 1 to skip graphical representations
%
% OUTPUTS_____
%
% fo (faces) and vo (vertices) can be used for patch visualization
%
% ROICROP is the 3D matrix with shifts applied to compensate for
% kinetochore plate deformations; sum(sum(ROICROP,1),3) integrates signal in planes orthogonal to (warped) spindle axis
% ROICROPdummy is the equivalent 3D matrix without the correction applied.

%
% 'Out' structure fields
% outdata.filename
%        .background
%        .thickness
%        .polshift
%        .ipMTproportion
%
% NOTE________
% Each sheetname in the 'data.xls' file should have the data of a corresponding
% tiff (single tiff file): sheetname = tiff name @ same folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global zrange
global yrange
global xrange
global nelements;

global xyres
global ROICROP
global ROICROPdummy



%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toggle = -1 ; % +1 or -1, whichever aligns the spindle axis with the z-axis
xyres = 0.04 ; % [um] xy pixel size (@ microscope acquisition)
zres = 0.38 ; % [um] z pixel size (@ microscope acquisition)

% % Visualization-only parameters
undensify = 1e-1 ; % patches reduction factor
kernel = 1 ; % [pixels] odd number; for gaussian filter smoothing
nmax = 100 ; disp( ['max # of (ImageJ) points per plate is set to ', num2str(nmax)] )
isovalues = [50]; % [photon count] threshold: a vector of 1 to 3 elements (3 layers max to be patched)


% % For quantifications of integrated density
bckg = background ; % [photon count] for background subtraction
zrange = 6 %semi [micron] Spindle Axis, 2nd dim

thickness = 0.2 ; % [micron] the thickness of the layers where signal will be integrated for the very final calculation
% These 3 layerss are i) +/- thickness/2 around the virtual
% equator surface spline and ii) the two surfaces
% dislocated poleward (by ''polshift'') from the plate
% surfaces (left and right), with the same thickness
polshift = thickness/2 + 0.1 ; % [micron] see description above
%%%%% end of PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
[ fname, path ] = uigetfile('*.tif','Select TIFF ') ;
info = imfinfo(fname);
pathxls=[path,'\data.xlsx'];
prompt={'Isovalue (Visualization only) [photons]  Higher is much faster' , 'Background level [photons]'};
def = {'15','5'};
dlgTitle = fname;
lineNo = [ones(2,1),ones(2,1)*40];
dialog = inputdlg(prompt,dlgTitle,lineNo,def);
isovalues = str2num(dialog{1});
bckg = str2num(dialog{2});


%% Image original size [pixels]
width = info.Width;
height = info.Height;
depth = numel(info);

dimx = width*xyres;
dimy = height*xyres;

%% assign imported data (1st excel column disregarded; it's just
% the xyz points' sequence order (in ImageJ)); convert slice to space

in = xlsread( pathxls , fname(1:end-4), ['B2:H' , num2str(nmax)] ) ;

%% pole and PLATE '1'
x1pole = in(1,1);
y1pole = in(1,2);
z1pole = in(1,3) * zres;  % Calibrate 'z'
pole1 = [x1pole, y1pole, z1pole]; % all in um

x1 = in(2:end,1);
y1 = in(2:end,2);
z1 = in(2:end,3) * zres;  % Calibrate 'z'
kts1 = [x1 y1 z1] ; % [um] position vector relative to image corner !...


%% pole '2' and PLATE '2'
x2pole = in(1,5);
y2pole = in(1,6);
z2pole = in(1,7) * zres;  % Calibrate 'z'
pole2 = [x2pole, y2pole, z2pole];

x2 = in(2:end,5);
y2 = in(2:end,6);
z2 = in(2:end,7) * zres;  % Calibrate 'z'
kts2 = [x2 y2 z2] ; % [um] position vector relative to image corner !...

if z1pole-z2pole ~= 0
    disp('poles cannot be in different slices')
end

% coordinates [um] of a point inside the virtual metaphase plane.
midpoint = [ (pole1(:,1)+pole2(:,1))./2 , (pole1(:,2)+pole2(:,2))./2  , (pole1(:,3)+pole2(:,3))./2 ] ;

% a vector normal to the plane (pole to pole vector divided by its norm).
midplanenormal = (pole1-pole2) ./ sqrt( sum((pole1-pole2).^2, 2) ) ;

%% angles, distances
%interpolar angle (deg) to image x-axis and distance
theta = 180/pi * acos( dot([1 0 0], midplanenormal) ) / norm([1 0 0])
interpolardistance = 2 *dot ( pole1-midpoint , midplanenormal ,2 )

%% Rotation of the point cloud defined in ImageJ
theta = toggle*theta ; %* TOGGLE ; % rotation angle
rotmatrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)] ;

% Offsets and 'UNoffsets'; to mantain coherence with 'imrotate' behaviour
% below (axis of rotation centered)
off1 = repmat( [dimx/2 dimy/2] , numel(kts1(:,1)),1 );
off2 = repmat( [dimx/2 dimy/2] , numel(kts1(:,2)),1 );

xy1new = ( (kts1(:,1:2)-off1) * rotmatrix' ) + off1 ; % xy matrix rotation ;
xy2new = ( (kts2(:,1:2)-off2) * rotmatrix' ) + off2 ; %

pole1new = ( ( pole1(1:2)-off1(1,:) ) * rotmatrix' ) + off1(1,:) ;
pole2new = ( ( pole2(1:2)-off2(1,:) ) * rotmatrix' ) + off2(1,:) ;

%% FINAL PLATES (points cloud) with axes reordering !!
% following the convention:
% x, optical axis (was z slices @ microscope) ;
% z, spindle axis (was xnew ;i.e. the 'xy' interpolar axis @microscope)
% y, the remaining dimension

plate1 = [ z1 xy1new(:,2) xy1new(:,1) ] ;
plate1pole = [ pole1(3), pole1new(2), pole1new(1) ] ;

plate2 = [ z2 xy2new(:,2) xy2new(:,1) ];
plate2pole = [ pole2(3), pole2new(2), pole2new(1) ] ;


%% Query for disregarding (margin=0) the first two and the last two datapoints (the 'margin domain') in both plates
margin=1;
if margin == 0
    plate1 = plate1(3:end-2,:) ;
    plate2 = plate2(3:end-2,:) ;
end
%% Prepare fit and fit splines
warning off
[xData1, yData1, zData1] = prepareSurfaceData( plate1(:,1), plate1(:,2), plate1(:,3) );
[xData2, yData2, zData2] = prepareSurfaceData( plate2(:,1), plate2(:,2), plate2(:,3) );

[fitresult1, gof1, output1] = fit( [xData1, yData1], zData1, 'thinplateinterp'  );
[fitresult2, gof2, output2] = fit( [xData2, yData2], zData2, 'thinplateinterp'  );


%% create fine grid for plate visualization.
res = 50;

[xout1 yout1] = meshgrid( linspace(min(plate1(:,1)),max(plate1(:,1)),res) , linspace(min(plate1(:,2)),max(plate1(:,2)),res) );
zout1 = fitresult1 ( xout1, yout1 ) ;
zout2 = fitresult2 ( xout1, yout1 ) ;

figure,
for n=1:2
    %% figure( 'Name', 'splines' );
    subplot(1,3,n)
    h1 = surf(xout1, yout1, zout1, zout1); shading interp, alpha .7
    hold on
    h2 = surf(xout1, yout1, zout2, zout2); shading interp, alpha .7
    
    view(90,0)
    
    colormap(hsv)
    box on
    
    camlight right
    lighting phong
    material shiny
    plot3 ( plate1pole(1),plate1pole(2),plate1pole(3),'ko' , 'MarkerSize', 5), plot3 ( plate2pole(1),plate2pole(2),plate2pole(3),'ko' ),
    hold on, axis equal
    plot3 ( plate1pole(1),plate1pole(2),plate1pole(3),'ko'), plot3 ( plate2pole(1),plate2pole(2),plate2pole(3),'ko' ),
    xlabel( 'x' );    ylabel( 'y' );
    grid on
end

%%                 -- thinplates (alone) were fit and displayed --
%% 2. Add the 3D spindle patch representation
% note: spindle will be i) interpolated along the 'x' axis (i.e. along microscope
% slices) just to have cubic voxels and ii) smoothened by the gaussian
% kernel with width defined in the 1st lines of this code


aspect = [ 1 1 xyres/zres ];
pixeldims = [ width height depth ];
microndims = pixeldims.*[ xyres xyres zres ];
volume = [ 1:microndims(1) 1:microndims(2) 1:microndims(3) ] ;

A = NaN( height, width, depth);

for k = 1 : depth
    A(:,:,k) = imrotate(imread(fname, k, 'Info', info), -theta, 'bilinear','crop');
end

A = A-2^15; % [photon counts] just to revert it to signed 16-bit pseudoformat
Asmooth = smooth3(A,'gaussian',kernel);

[X Y Z] = meshgrid(1:width,1:height,1:depth) ;
[Xquery Yquery Zquery] = meshgrid( 1:width, 1:height , linspace(1,depth,depth/aspect(3)) ) ;
Asmoothinterp = interp3( X, Y, Z, Asmooth, Xquery, Yquery, Zquery ) ;

XX = Xquery * xyres;
YY = Yquery * xyres;
ZZ = Zquery * zres;

limits = [NaN NaN NaN NaN NaN NaN] ;
%[x, y, z, subA] = subvolume(Asmoothinterp, limits);
subA = Asmoothinterp;

if threeDrepresentationsFlag==1
    
    for isovalue = 1:numel(isovalues)
        [fo,vo] = isosurface( ZZ, YY, XX, subA, isovalues(isovalue) );
        for m=[2 3]
            disp('iso')
            subplot(1,3,m)
            p1 = patch('Faces', fo, 'Vertices', vo);
            axis([xout1(1) xout1(end) yout1(1) yout1(end) ])
            axis equal
            reducepatch(p1,undensify); % if PC complains
            switch isovalue
                case 1
                    set(p1,'FaceColor',[.6 .6 .6]); set(p1,'FaceAlpha',.5);
                case 2
                    set(p1,'FaceColor',[.3 .9 .3]); set(p1,'FaceAlpha',.5);
                case 3
                    set(p1,'FaceColor',[1 .6 .2]); set(p1,'FaceAlpha',.5);
            end
            set(p1,'EdgeColor','none');
            
            xlabel( 'x (um)' );
            ylabel( 'y (um)' );
            zlabel( 'z (um)' );
            
            view(90,0)
            
            colormap(hsv)
            box on
            
            
            camlight(45,45)
            lighting phong
            material shiny
        end
    end
    title(fname(1:end-4))
elseif threeDrepresentationsFlag==0
    fo=[];
    vo=[];
end
%% 3. INTENSITY PROFILES (Kt plate 1 is the master plate; Kt plate 2 is cropped or extrapolated to match the master)

% definition of subregion XY (i.e. equator-plane!) for quantifications
% By convention,PLATE ONE's xy domain was chosen to define (constrain or
% extrapolate) the domain of plate2

[Xm Ym] = meshgrid ( min(plate1(:,1)) : xyres : max(plate1(:,1)) , min(plate1(:,2)) : xyres : max(plate1(:,2)) ) ;
s = size(Xm)

Zm1 = fitresult1(Xm, Ym) ; % plate1 domain; all micron
Zm2 = fitresult2 (Xm, Ym) ; % plate1 (ONE) domain; all micron
ZmIntra = 1/2 * (Zm1 + Zm2) ;
absplatedist = (abs(Zm1-Zm2)) ;

distZm12 = Zm1-Zm2;
disp(['avg (Zm1-Zm2) = ',num2str( mean(distZm12(:)) )  ]);
disp(['min (Zm1-Zm2) = ',num2str( min(distZm12(:)) )  ]);


% Check for plate intersections and, if not, check if Kt plates (Zm1 and
% Zm2) should be switched. This is to warrant that the calculated poleward
% volumes (below) are indeed poleward and not shifted towards the equator.
if min(distZm12(:))<0
    if max(distZm12(:))>0
        disp(' Computation aborted (Kinetochore plates are intersecting).')
        return
        
    elseif max(distZm12(:))<0
        Zmtemp=Zm1;
        Zm1=Zm2;
        Zm2=Zmtemp;
    end
end

subplot(1,3,3), hold on
h.a = surf(Xm, Ym, Zm1, 'EdgeColor',[.3 .3 0],'EdgeAlpha',.2,'FaceColor',[0 1 0], 'FaceAlpha',.5), hold on,
h.b = surf(Xm, Ym, Zm2, 'EdgeColor',[.3 .3 0],'EdgeAlpha',.2,'FaceColor',[0 1 0], 'FaceAlpha',.5)

fcolor = [0 0.3 0.5];
surf(Xm, Ym, ZmIntra+thickness/2,'EdgeColor',[.3 .0 0],'EdgeAlpha',.2,'FaceColor',fcolor, 'FaceAlpha',.5),
surf(Xm, Ym, ZmIntra-thickness/2,'EdgeColor',[.3 .0 0],'EdgeAlpha',.2,'FaceColor',fcolor, 'FaceAlpha',.5),
surf(Xm, Ym, Zm1+polshift,'EdgeColor',[.3 .2 0],'EdgeAlpha',.2,'FaceColor',fcolor, 'FaceAlpha',.5),
surf(Xm, Ym, Zm1+polshift+thickness,'EdgeColor',[.3 .2 0],'EdgeAlpha',.2,'FaceColor',fcolor, 'FaceAlpha',.5),
surf(Xm, Ym, Zm2-polshift,'EdgeColor',[.3 .2 0],'EdgeAlpha',.2,'FaceColor',fcolor, 'FaceAlpha',.5),
surf(Xm, Ym, Zm2-polshift-thickness,'EdgeColor',[.3 .2 0],'EdgeAlpha',.2,'FaceColor',fcolor, 'FaceAlpha',.5),

ROI = subA ;
ROI = ROI-bckg ;
ROI ( find(ROI<0) ) = 0 ;

ROI = ROI ( round((Ym(1):xyres:Ym(end))/xyres) , : , round((Xm(1):xyres:Xm(end))/xyres) ) ;

sROI = size( ROI ) ;
disp(['ROI size (plate ''1'' domain + 2x the defined z half-range)  is ',mat2str(sROI),' pixels = ',mat2str(sROI*xyres),' micron'])

ROIxycrop=ROI;
s = size( ROIxycrop ) ;
disp(['final cropped matrix size was ',mat2str(s)])

ROICROP = zeros( size(ROIxycrop,1), 2*zrange/xyres+1, size(ROIxycrop,3) ) ;
ROICROPdummy = ROICROP ;
%leftvolume = zeros( size(ROIxycrop,1), 2*thickness/xyres+1, size(ROIxycrop,3) ) ;
%rightvolume = leftvolume;
%centervolume = leftvolume;



%% Calculation
figure, hold on
for xpix = 1:1:s(3)
    for ypix = 1:1:s(1)
        
        shift(ypix,xpix) = ZmIntra( ypix  , xpix )  ;    % [pixel!]
        ROICROP(ypix,:,xpix) = ...
            ROIxycrop( ypix , ...
            (ceil((shift(ypix,xpix)-zrange)/xyres ...
            :(shift(ypix,xpix)+zrange)/xyres)) , ...
            xpix);
        
        shiftdummy = round(s(2)/2*xyres); ROICROPdummy(ypix,:,xpix) = ROIxycrop( ypix , (ceil((shiftdummy-zrange)/xyres:(shiftdummy+zrange)/xyres))' , xpix);
        
        leftvolume(ypix,:,xpix) = ...
            ROIxycrop( ypix , ...
            ceil( (Zm2(ypix,xpix) - polshift - thickness) /xyres) ...
            : ceil( (Zm2(ypix,xpix) - polshift ) /xyres)  ...
            , xpix );
        
        rightvolume(ypix,:,xpix) = ...
            ROIxycrop( ypix , ...
            ceil( (Zm1(ypix,xpix) + polshift ) /xyres) ...
            : ceil( (Zm1(ypix,xpix) + polshift + thickness) /xyres)  ...
            , xpix );
        
        centervolume(ypix,:,xpix) = ...
            ROIxycrop( ypix , ...
            ceil( (ZmIntra(ypix,xpix)-thickness/2)/xyres) ...
            : ceil( (ZmIntra(ypix,xpix)+thickness/2)/xyres)  ...
            , xpix);
        
        
        
        %             subplot(2,1,1)
        %             plot( -zrange:xyres:zrange, ROICROP(ypix,:,xpix) )
        %             subplot(2,1,2)
        %             plot( -zrange:xyres:zrange, ROICROPdummy(ypix,:,xpix),'r-')
        
    end
end

%%



figure,
subplot(2,2,1), hold on, title('Corrected'), hold on
yrange = xyres*(size(ROICROP,1)-1)/2;
imagesc(-zrange:xyres:zrange, -yrange:xyres:yrange, (sum(ROICROP,3)))
impixelinfo, axis equal, axis xy, axis tight
xlabel('[um]'), ylabel('[um]')

subplot(2,2,3), hold on
% Display max and min inter-plate dist (red-shaded areas)
edges = max(absplatedist(:)) ;
edgesmin = min(absplatedist(:)) ;
area([-edgesmin/2 -edges/2],[1 1],'FaceColor',[1 .5 .5],'EdgeC','none') ;
area([edgesmin/2 edges/2],[1 1],'FaceColor',[1 .5 .5],'EdgeC','none') ;

intprofile = sum(sum(ROICROP,1),3) ;
intprofilenorm = intprofile/max(intprofile(:)) ;
plot(-zrange:xyres:zrange , intprofilenorm , 'Color',[.3 .3 .3], 'Linewidth',3);
xlabel( 'along spindle axis (um)' )
ylabel( 'intensity (norm.)' )

plot( [-thickness/2 -thickness/2],[0 1],'k'), plot([+thickness/2 +thickness/2],[0 1],'k')
[nelements, centers]=hist(distZm12(:)/2,100);
area(centers,.3*nelements/max(nelements(:)))
[nelements, centers]=hist(-distZm12(:)/2,100);
area(centers,.3*nelements/max(nelements(:)))

% Dummy (just to show the uncorrected spindle)
subplot(2,2,2), title('Uncorrected'), hold on
imagesc(-zrange:xyres:zrange, -yrange:xyres:yrange, (mean(ROICROPdummy,3)));
impixelinfo
axis equal, axis tight
xlabel('[um]')
ylabel('[um]')

subplot(2,2,4),
intprofiledummy = sum(sum(ROICROPdummy,1),3);
intprofiledummynorm = intprofiledummy / max(intprofiledummy(:));
plot(-zrange:xyres:zrange,intprofiledummynorm,'Color',[.3 .3 .3],'Linewidth',3); ylim([0 1])
xlabel('along spindle axis (um)')
ylabel('intensity (norm.)')
%% Proportion of interpolar MTs relative to all MTs
%Corrected
% ktplusipMTs = sum( intprofile( zrange/xyres + round((-edges/2-polshift-thickness/2 : -edges/2-polshift+thickness/2) / xyres )) ) ;
% ipMTs = sum( intprofile( zrange/xyres + round((-thickness/2 : +thickness/2) / xyres )) ) ;
% ipProportion_corrected = ipMTs / ktplusipMTs ;
%
% %Uncorrected
% ktplusipMTs = sum( intprofiledummy( zrange/xyres + round((-edges/2-polshift-thickness/2 : -edges/2-polshift+thickness/2) / xyres )) ) ;
% ipMTs = sum( intprofiledummy( zrange/xyres + round((-thickness/2 : +thickness/2) / xyres )) ) ;
% ipProportion_uncorrected = ipMTs / ktplusipMTs;
% ipProportion_uncorrected_fromMin = min(intprofiledummynorm( zrange/xyres + round((-edges/2-polshift-thickness/2 : +edges/2+polshift+thickness/2) / xyres )) );

prop = 2* sum(sum(centervolume,1),3) / sum(sum(leftvolume+rightvolume,1),3);

outdata.filename = fname;
outdata.background = bckg;
oudata.thickness = thickness;
oudata.thickness = thickness;
outdata.ipMTproportion = prop;

%% Display spline heatmaps/contours
figure,
subplot(1,2,1)
imagesc(Xm(1,:),Ym(:,1),ZmIntra), axis equal, impixelinfo, title ('Virtual equator (ZmIntra)')
xlabel( 'x (um) [along slices]' );
ylabel( 'y (um)' );

subplot(1,2,2), hold on
imagesc( Xm(1,:), Ym(:,1), distZm12) %distZm12 was defined above as Zm1-Zm2
contourlevels = [0.1*ceil(10*min(distZm12(:))), 0:.5:5, 0.1*ceil(10*max(distZm12(:)))];
contour( Xm(1,:), Ym(:,1), distZm12, contourlevels ,'showtext','on' ), axis equal, impixelinfo, title (['Interplate distance ( Zm1-Zm2 ) ; max is ', num2str(max(absplatedist(:))), 'um']);

xlabel( 'x (um) [along slices]' );
ylabel( 'y (um)' );

ss=size(ROICROP);
y = 1:ss(1);
try
    for n=1:ss(3)
        imwrite(uint16([ROICROP(:,:,n) ; ROICROPdummy(:,:,n)]),['ROICROP', sprintf('%03d',n), '.tif'],'tiff','Compression','none')        
    end
catch
    msg = sprintf('Couldn''t save TIFFs. ');
end
