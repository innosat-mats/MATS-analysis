function fig = MLTplotfunc(type,stripsNH,stripsSH,peaksNH,peaksSH)
    % Arguments: type = 1 or 2, strips = matlab structure allstrips file
    % peaks = matlab structure peak file
    if type == 1
        sgtitle({'\bf MLT plot of peak points correlated with intensity';'\rm Period: February 8th to 28th'},fontsize=16)
    
        %%%% North Hemisphere polar plot
        subplot(1,2,1);
       
        %plot settings
        p1 = polarscatter(peaksNH.MLT*(pi/12),peaksNH.Mlat,[],peaksNH.maxI,'filled');
        set(gca,"CLim",[2*10^3 3*10^4])  %Altitude limit = 90 115,  kplimit = 0-9 intensity limit = 1*10^3 3*10^4
        cb = colorbar ;
        colormap jet
        cb.Label.String = '10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)'; %'10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)' 'Altitude (km)' 'Kp-index'
        
        ax1 = gca;
        % MLT (theta settings)
        ax1.ThetaDir = 'counterclockwise' ;
        ax1.ThetaZeroLocation = 'bottom' ;
        ax1.ThetaTick = [0,45,90,135,180,225,270,315,360];
        ax1.ThetaTickLabel = {'00';'03';'06';'09';'12';'15';'18';'21'} ;
        % Latitude settings
        ax1.RLim = [40, 90] ;
        ax1.RDir = "reverse" ;
        rtickformat(ax1,"degrees") 
        rticks(40:10:90);
        
        %Settings for the peak points data tips
        p1.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",peaksNH.MLT);
        p1.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",peaksNH.Mlat);
        p1.DataTipTemplate.DataTipRows(3) = dataTipTextRow("Altitude",peaksNH.alt);
        p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Kp",peaksNH.kp);
        p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Time",peaksNH.time);
        p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Intensity",peaksSH.maxI);
        
        subtitle('Northern Hemisphere','FontWeight','bold',FontSize=12);
        
        %%%% South Hemisphere polar plot
        subplot(1,2,2) ;
        p2 = polarscatter(peaksSH.MLT*(pi/12),peaksSH.Mlat,[],peaksSH.maxI,'filled'); %change the last attribute as see fit
        set(gca,"CLim",[2*10^3 3*10^4])  %[0 9] , [2*10^3 4*10^4] , [90 115]
        cb = colorbar ;
        cb.Label.String = '10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)';  %'10^{13} photons / (nm \cdot m^2 \cdot sr \cdot s)'; 'Altitude (km)';'Kp-index' ; 

        ax2 = gca;
        ax2.ThetaDir = 'counterclockwise' ;
        ax2.ThetaZeroLocation = 'bottom' ;
        ax2.ThetaTick = [0,45,90,135,180,225,270,315,360];
        ax2.ThetaTickLabel = {'00';'03';'06';'09';'12';'15';'18';'21'} ;
        ax2.RLim = [-90, -40] ;
        ax2.RDir = "normal" ;
        rtickformat(ax2,"degrees")
        rticks(ax2,(-90:10:-40));
        
        %Settings for the peak points data tips
        p2.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",peaksSH.MLT);
        p2.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",peaksSH.Mlat);
        p2.DataTipTemplate.DataTipRows(3) = dataTipTextRow("Altitude",peaksSH.alt);
        p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Kp",peaksSH.kp);
        p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Time",peaksSH.time);
        p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Intensity",peaksSH.maxI);

        subtitle('Southern Hemisphere','FontWeight','bold',FontSize=12);

    end

    if type == 2
        sgtitle({'\bf Peak points correlated with Kp-level and TP-path';'\rm Period: May 1st to 8th'},fontsize=16)
        
        %%%% North Hemisphere polar plot
        subplot(1,2,1);
        
        %plot settings
        s1 = polarscatter(stripsNH.MLT*(pi/12),stripsNH.Mlat,0.8,'filled','black','MarkerFaceAlpha',.4);
        hold on
        p1 = polarscatter(peaksNH.MLT*(pi/12),peaksNH.Mlat,[],peaksNH.kp,'filled');
        
        set(gca,"CLim",[0 9])  
        cb = colorbar ;
        cb.Label.String = 'Kp-index' ;
        ax1 = gca;

        % MLT (theta settings)
        ax1.ThetaDir = 'counterclockwise' ;
        ax1.ThetaZeroLocation = 'bottom' ;
        ax1.ThetaTick = [0,45,90,135,180,225,270,315,360];
        ax1.ThetaTickLabel = {'00';'03';'06';'09';'12';'15';'18';'21'} ;
        % Latitude settings
        ax1.RLim = [40, 90] ;
        ax1.RDir = "reverse" ;
        rtickformat(ax1,"degrees") 
        rticks(40:10:90);
        
        %Settings for the peak points data tips
        p1.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",peaksNH.MLT);
        p1.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",peaksNH.Mlat);
        p1.DataTipTemplate.DataTipRows(3) = dataTipTextRow("Altitude",peaksNH.alt);
        p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Kp",peaksNH.kp);
        p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Time",peaksNH.time);
        
        %Settings for the path data tips
        s1.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",stripsNH.MLT);
        s1.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",stripsNH.Mlat);
        
        subtitle('Northern Hemisphere','FontWeight','bold',FontSize=12);
        
        %%%% South Hemisphere polar plot
        subplot(1,2,2) ;
        
        s2 = polarscatter(stripsSH.MLT*(pi/12),stripsSH.Mlat,0.8,'filled','black','MarkerFaceAlpha',.3);
        hold on
        p2 = polarscatter(peaksSH.MLT*(pi/12),peaksSH.Mlat,[],peaksSH.kp,'filled');
        
        set(gca,"CLim",[0 9])  
        cb = colorbar ;    
        cb.Label.String = 'Kp-index' ;
        ax2 = gca;

        ax2.ThetaDir = 'counterclockwise' ;
        ax2.ThetaZeroLocation = 'bottom' ;
        ax2.ThetaTick = [0,45,90,135,180,225,270,315,360];
        ax2.ThetaTickLabel = {'00';'03';'06';'09';'12';'15';'18';'21'} ;
        ax2.RLim = [-90, -40] ;
        ax2.RDir = "normal" ;
        rtickformat(ax2,"degrees")
        rticks(ax2,(-90:10:-40));
        
        %Settings for the peak points data tips
        p2.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",peaksSH.MLT);
        p2.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",peaksSH.Mlat);
        p2.DataTipTemplate.DataTipRows(3) = dataTipTextRow("Altitude",peaksSH.alt);
        p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Kp",peaksSH.kp);
        p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Time",peaksSH.time);
        
        %Settings for the path data tips
        s2.DataTipTemplate.DataTipRows(1) = dataTipTextRow("MLT",stripsSH.MLT);
        s2.DataTipTemplate.DataTipRows(2) = dataTipTextRow("MLat",stripsSH.Mlat);
        
        subtitle('Southern Hemisphere','FontWeight','bold',FontSize=12);
    end
end