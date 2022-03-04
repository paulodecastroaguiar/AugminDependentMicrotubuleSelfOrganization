% Read line profiles across kt pairs
% Calculate inter/kt distance and kt width
path='C:\ .... \test.xlsx'; % fill with folder

ctrl = xlsread(path,2,'A3:AL92');
haus6 = xlsread(path,3,'A3:AL92');
ndc80 = xlsread(path,4,'A3:AL92');

% create a sufficently large mosaic for the number of graphs
plotdim1=5;
plotdim2=6;

xpos = ctrl(:,1) ;
xstep = xpos(2)-xpos(1) ;
npos = numel( xpos ) ;

nctrl = size(ctrl,2)-1 ;
nhaus6 = size(haus6,2)-1 ;
nndc80 = size(ndc80,2)-1 ;

norma=repmat( max(ctrl(:,2:end)), npos, 1 );
ctrl(:,2:end) = ctrl(:,2:end)./norma;
norma=repmat( max(haus6(:,2:end)), npos, 1 ) ;
haus6(:,2:end) = haus6(:,2:end)./norma ;
norma=repmat( max(ndc80(:,2:end)), npos, 1 ) ;
ndc80(:,2:end) = ndc80(:,2:end)./norma ;

kernelwidth = npos/25 ;
kernelcenter = npos/2 ;
peakrange = 15 ; %pixels, for 2nd order polyfit 

% Scan controldata
locs = [];
cv = [];
mark_1=[];
mark_2=[];

list = { 'ctrl' , 'haus6' } ;


for ind = 1 : 2
    figure
    data = eval(list{ind});
    list{ind};
    for n = 1 : eval(strcat('n',list{ind}))
          
        subplot(plotdim1,plotdim2,2*n)
        plot( data(:,1+n) )        
        
        cv(:,n+1) = conv( data(:,1+n) , gaussmf(1:npos,[kernelwidth, kernelcenter]),'same' );
        
        [pks , locs(:,n+1)] = findpeaks(cv(:,n+1)./max(cv(:,n+1)),'NPEAKS',2,'MINPEAKHEIGHT',.5, 'MINPEAKDIST',15);
        
        
        %subplot(plotdim1,plotdim2,2*n+1)
        %plot(cv(:,n+1))
        %hold on
        
        % Parabolic fits on both peaks
        xpeak1 = round( locs(1,n+1)-peakrange/2 : locs(1,n+1)+peakrange/2 );
        xpeak1 = xpeak1( xpeak1>0 );
        
        xpeak2 = round( locs(2,n+1)-peakrange/2 : locs(2,n+1)+peakrange/2 );
        xpeak2 = xpeak2( xpeak2>0 );
        
        fit_1=polyfit(xpeak1',data(xpeak1,n+1),2);
        fit_2=polyfit(xpeak2',data(xpeak2,n+1),2);
        
        mark_1(ind,n+1) = -fit_1(2)/(2*fit_1(1));
        conc_1(ind,n+1) = -1./(2*fit_1(1));
        
        yfit_1=polyval(fit_1,xpeak1);
        hold on
        subplot(plotdim1,plotdim2,2*n)
        hold on
        plot(xpeak1,yfit_1,'r')
        
        mark_2(ind,n+1) = -fit_2(2)/(2*fit_2(1));
        conc_2(ind,n+1) = -1./(2*fit_2(1));
         
        yfit_2 = polyval(fit_2,xpeak2);
        hold on
        subplot(plotdim1,plotdim2,2*n)
        hold on
        plot(xpeak2,yfit_2,'r')
        
        dists=[mark_2-mark_1]'*xstep;
        dists(dists==0)=NaN; 
        conc_1(conc_1==0) = NaN; 
        conc_2(conc_2==0) = NaN; 
    end
        annotation('textbox', [0.2, 0.05, .6, 0.05], 'String', 'peaks below 50% of twin peak are dismissed')

end

    toclipboard([dists conc_1' conc_2'])
    
