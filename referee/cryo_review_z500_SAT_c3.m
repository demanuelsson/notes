%%%%%%%%%%%%%%%%%%%%%%%%%
% cryo_review_z500_SAT_c3.m
% Daniel Emanuelsson
% Matlab 2017a
%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * corrcoef_df.m   UoW Steig
% * cosweight.m     UoW online archive Atmospheric Science
% * annave.m        UoW online archive Atmospheric Science
% * star_coord_WA.m DE
% * fig.m           fileexchange
% * export_fig.m    fileexchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close

site='TA';

cor_nr=1;  % (14) d18O (2) d-excess %%%%%%%%%%%%%%%%%
if cor_nr==1
      iso='d18O';
elseif cor_nr==2     
      iso='d_excess';
end

%%%%%%%%%%%%%%%%%%%%%%%
crop_nr=0; % (1/0) crop?

lock_scalebar=1; % (1/0)
%rcontour_p=0;% (1/0) SAT reg contour

% corr_time=1; % (1) annual, (2) monthly
 show_max=0;   % (1/0) show max label
 max45=0;

coast_nr=1; % (1/0) On/Off, slow turn off these if you just want to check something, but have them on when you generate final figure 

%%%%%%%%%%Figure format
figure_format=2; % (1) - EPS, (2) - PNG 
%%%%%%%%%%
corr_label=1; %(1/0) %%%%%%%%%%% Show colorbar and label

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_value=1; % sign (1)  positive (2) negative r-values in the ABS/Ross Sea region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era_name_nr=1; % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. z500 *
% 2. u850
% 3. v850
% 5. 2mT *

if era_name_nr==1
    name='z500';
elseif era_name_nr==2
    name='u850';
elseif era_name_nr==3
    name='v850';
elseif era_name_nr==5
    name='2mT'; 
end


age_scale_nr=1;%%%%%(1) new (2) old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if era_name_nr==1

                if age_scale_nr==1
                    %yr_e=2012;
                    yr_e=2014;
                    yr_s=1997;
                    letter='a';
                    age_scale_str='New';
                elseif age_scale_nr==2
                    letter='c';
                    yr_e=2014;
                    yr_s=1998;
                    age_scale_str='Original';
                end

            elseif era_name_nr==3 || era_name_nr==2
                
                if age_scale_nr==1
                  %  yr_e=2012;
                    yr_e=2014;
                    yr_s=1997;
                    if era_name_nr==3
                        letter='a';
                    elseif era_name_nr==2
                        letter='b';
                    end
                    age_scale_str='New';
                elseif age_scale_nr==2
                    
                    yr_e=2014;
                    yr_s=1998;
                    
                    if era_name_nr==3
                        letter='c';
                    elseif era_name_nr==2
                        letter='d';
                    end                
                    
                    age_scale_str='Orig.';
                end                


            elseif era_name_nr==5
%                 letter='d';
                if age_scale_nr==1
                  %  yr_e=2012;
                    yr_e=2014;
                    yr_s=1997;
                    letter='b';
                    age_scale_str='New';
                elseif age_scale_nr==2
                    letter='d';
                    yr_e=2014;
                    yr_s=1998;
                    age_scale_str='Original';
                end
            end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sea_nr=1;                  Season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    season='annual';
    
proj='stereo';
 



%       Size:       480x241x421
%        Dimensions: longitude,latitude,time
 addpath C:\PHD\ERA_interim

% ncdisp('ERA_int_monthly_z500.nc') 
% ncdisp('ERA_int_monthly_2m_T.nc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 0 Read in ERA-interim data
addpath C:\PHD\ERA_interim

name_c='ERA_int_monthly_z500_2.nc';

era_time=ncread(name_c,'time');
era_long=ncread(name_c,'longitude');
era_lat=ncread(name_c,'latitude');

% time
dayssince111=era_time/24;
datevalue=dayssince111+datenum(1900,1,1);
 date_vec=datevec(double(datevalue)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 era_year_num=yyyy+(mm/12-1/12); 
%

if strcmp(name,'2mT')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % surface 2mT 
    name_c='ERA_int_monthly_2m_T_c.nc';
    era_T=ncread(name_c,'t2m'); % 2m temp
    era_T=era_T- 273.15;
%     letter='b';
    lat1=-90; 
    lat2=-30;  
 
 
 elseif strcmp(name,'z500')==1 && strcmp(season,'annual')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % geopotentail 500 hPa  
     
     
    name_c='ERA_int_monthly_z500_2.nc'; 
    era_z500=ncread(name_c,'z'); 
    era_z500=era_z500/9.80665;
  
        if strcmp( proj,'stereo')==1
        lat1=-90;
        lat2=-10;
        end
        
elseif strcmp(name,'v850')==1 || strcmp(name,'u850')==1        
%  era_v850=ncread('ERA_int_monthly_v850.nc','v'); % v850 'm s**-1'
 

    if strcmp(name,'v850')==1 
        par_str='v';
    elseif strcmp(name,'u850')==1
        par_str='u';
    end


 
%============================================
% using newer files, each file is one year

            era_str=name;

 
            ERA_I_X_store=[];
            era_time_store=[];
            for i=1979:2016

            ERA_I_X=ncread(['ERA_interim_',era_str,'_monthly_means_of_daily_means_',num2str(i),'_c.nc'],par_str); % 1979
            

            era_time=ncread(['ERA_interim_',era_str,'_monthly_means_of_daily_means_',num2str(i),'_c.nc'],'time'); 
            
            % help cat
            ERA_I_X_store=cat(3,ERA_I_X_store,ERA_I_X);
            era_time_store=vertcat(era_time_store,era_time);
            
            end     
        
        era_v850=ERA_I_X_store;
 
    if strcmp(name,'v850')==1 
        era_v850=ERA_I_X_store;
    elseif strcmp(name,'u850')==1
        era_u850=ERA_I_X_store;
    end
        
        
 
    lat1=-90;
    lat2=-40; 

        
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%  elseif strcmp(name,'u850')==1  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  era_u850=ncread('ERA_int_monthly_u850.nc','u'); 
%  
   
 
 
 end




%
% era_count=396; % that overlaps with RICE record 1979-2011
era_count=find(era_year_num==yr_e)+11;  % data in months
%era_start=find(era_year_num==1980-1);

% annual 
era_start=find(era_year_num==yr_s);
div=ceil((era_count-era_start)/12);


    if strcmp(name,'z500')==1
        M=era_z500;
    elseif strcmp(name,'2mT')==1
        M=era_T;    
    elseif strcmp(name,'v850')==1
        M=era_v850;
    elseif strcmp(name,'u850')==1
        M=era_u850;
    end

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% ===================================================
 
 era_lat_c=era_lat;

 [nlon, nlat, ntim]= size(M);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of these three steps doesnt matter
%%% -=(1)=-
%  z = cosweight(z, lat);
% weighted to account for the decrease in area towards the pole.

%  era_z500_c = reshape(era_z500, ntim, nlat, nlon);  %  z500_eof

era_z500_c = permute(M, [3 2 1]);  %  z500_eof

%%%%%%%%%%%%%%%
%     B = permute(A,ORDER) rearranges the dimensions of A so that they
%     are in the order specified by the vector ORDER.  The array produced
%     has the same values as A but the order of the subscripts needed to 
%     access any particular element are rearranged as specified by ORDER.
%%%%%%%%%%%%%%%%%

era_z500_c2=cosweight(era_z500_c,era_lat_c); % Original UoW function    (time x lat x lon) 

%%%%%%%%%%%era_z500_c2=cosweight_c(era_z500,era_lat_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 %%% -=(2)=-
% Remove monthly climatology for all cells 

 era_z500_c3 = double(reshape(era_z500_c2, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33       until 33X240=7920
 %  
 %1
 %time
 %
 %444
 
 
 [era_z500_c4,clim_z500] = annave(era_z500_c3);   % checked, Removes
%  seasonal cycle

% Back to old format again
era_z500_c6= reshape(era_z500_c4,  ntim, nlat, nlon );

M = permute(era_z500_c6, [3 2 1]);  %  z500_eof
 
 %==========================================================
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Since the grid size decreases as you move towards the pole, 
% weight each grid box (i.e., multiply the time series at each grid box) 
% by the square root of the cosine of latitude (the weights are based on 
% the square root of the cosine so that the covariance matrix is weighted by the cosine of latitude).

% S_lon=266;S_lat=227;%check
% M ( S_lon, S_lat-20,2) 
%  M=cosweight_c(M,era_lat); % UoW function
%  M ( S_lon, S_lat-20,2) % check should changes value before and after

%            Size:       480x241x421
%            Dimensions: longitude,latitude,time
%

[m n t]=size(M(:,:,1:era_count)); %
for i= 1:m
    for j= 1:n

        

        dummy1=reshape(squeeze(M(i,j,era_start:era_count)), 12, div);          
%           
     % Annual    
    if  strcmp(season,'annual')==1
         
%          if corr_time==1
          ERA_M_annual=mean(dummy1,1); % annual ERA-values
          ERA_M_annual_reg(i,j,:)=ERA_M_annual(1:end); % save annual values for regression
          ERA_M_annual_dummy=detrend(ERA_M_annual(1:end)); % change to 1 to get 1979 annual
          ERA_M_annual_detrend(i,j,:)=ERA_M_annual_dummy;            
%          end
          
    end  
    
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isotopes record

if strcmp(iso,'d18O')==1 
% if corr_time==1    % Annual 

% load('C:\g_drive\review\Cryosphere\Data001_d18O_cryo_review.mat')
load('C:\g_drive\review\Cryosphere\MA_TA_d18O.mat')

% date_annual=MA_TA_d18O(:,1);

    if age_scale_nr==1 % New
        date_annual=MA_TA_d18O(:,1);
        stacked_record_annual_Ma=MA_TA_d18O(:,2);
    elseif age_scale_nr==2 % Original  
        date_annual=MA_TA_d18O(:,3);
        stacked_record_annual_Ma=MA_TA_d18O(:,4);
    end

start_t=find(date_annual==yr_s);
end_t=find(date_annual==yr_e); 

    
    X=detrend(stacked_record_annual_Ma((start_t:end_t),1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(iso,'d_excess')==1 
    
load('C:\g_drive\review\Cryosphere\MA_TA_d_excess.mat')

% date_annual=MA_TA_d_excess(:,1);    

    if age_scale_nr==1 % New
        date_annual=MA_TA_d_excess(:,1);
        stacked_record_annual_Ma=MA_TA_d_excess(:,2);
    elseif age_scale_nr==2 % Original  
        date_annual=MA_TA_d_excess(:,3);
        stacked_record_annual_Ma=MA_TA_d_excess(:,4);
    end

start_t=find(date_annual==yr_s);
end_t=find(date_annual==yr_e); 

    
 X=detrend(stacked_record_annual_Ma((start_t:end_t),1));


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERA-data choose annual or season 

if strcmp(season,'annual')==1
%    season='Annual';
   Y=ERA_M_annual_detrend;

end
  
[A B C]=size(Y);
R=zeros(2,2,A,B);P=zeros(2,2,A,B);
for i=1:A
    for j=1:B
        Q=squeeze(Y(i,j,:));
      [R(:,:,i,j) P(:,:,i,j)]=corrcoef_df(X,Q, 'rows','pairwise');
     % moving average
     %   [R(:,:,i,j) P(:,:,i,j)]=corrcoef_df(moving_average(X,5),moving_average(Q,5), 'rows','pairwise');
    end
end

% for RICE
R(:,:,find(floor(era_long)==199),find(floor(era_lat)==-79))
P(:,:,find(floor(era_long)==199),find(floor(era_lat)==-79))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure
% axesm define projection
f1=fig('units','inches','width',10.5,'height',10.5,'font','Helvetica','fontsize',18,'border','on'); 

 proj='stereo';
%proj='mercator';

   axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...
       'mlabellocation',[0:30:179,0:-30:-180]); 
   %Frame - around figure
   
  set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes  
  
 gridm('GLineStyle','--','Gcolor',[.6 .6 .6],'GLineWidth',1.5,...
    'MLineLimit',[lat1 lat2],'Galtitude', .02)

hold on

%%%%%%%%%%%%%%%%%%%%
c_min=squeeze(R(1,:,:));
c_min=c_min(2,:);
c_min=min(c_min);

c_max=squeeze(R(1,:,:)); 
c_max=c_max(2,:);
c_max=max(c_max);

c_max_c=c_max;

% Masking out areas that are not statistically significant
Rc=R;

% Rc2=R;

% find lat long max 
% c_max=max(c_max,[],2);
% find(c_max==c_max)   (54-lat,
%c_max=X_r(2,:);


X_r=squeeze(Rc(1,2,:,:))';
% X_r=X_r((181:241),:); % from -45S to -90 max



if max45==1

        X_r((1:181),:)=NaN;  %%% MAX below -45 degrees, and east ant E 0-150
        X_r(:,(1:201))=NaN; 
else
        X_r=X_r;  %%% MAX below -45 degrees, and east ant E 0-150
%         X_r(:,(1:201))=NaN; 
    
   
end


%X_r(:,(465:480))=NaN; 

if r_value==1
X_r_max=nanmax(nanmax(X_r));
X_r_long=nanmax(X_r); % 1 long
X_r_lat=nanmax(X_r,[],2); %2

elseif r_value==2
X_r_max=nanmin(nanmin(X_r));
X_r_long=nanmin(X_r); % 1 long
X_r_lat=nanmin(X_r,[],2); %2
end
    
xr=find(X_r_long==X_r_max); 
xr_lon=era_long(xr);

xr=find(X_r_lat==X_r_max); 
xr_lat=era_lat(xr);
 


c_limit=max(abs(c_min),abs(c_max));
c_limit_c=c_limit;


if lock_scalebar==1 && strcmp(name,'2mT') % && (strcmp(season,'annual')) && (strcmp(iso,'d18O'))
  %  c_limit=0.7484;
    c_limit=0.8;
%    c_limit=1;
elseif lock_scalebar==1 && strcmp(season,'annual') && strcmp(iso,'d18O')    
%     c_limit=0.7351;
        c_limit=0.8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rc(P>=1)=NaN;
% Rc(P>=0.1)=NaN; % maksing comment out if you dont want it

%%%%%%%%%%%%% fill in seam
Rc_c=squeeze(Rc(1,2,:,:))';
era_long_c=era_long;
era_lat_c=era_lat;
era_long_c(481,:)=359.99;
Rc_c(:,481)=Rc_c(:,1);
%%%%%%%%%%%%%%%%%%

hSurf=surfm(double(era_lat),double(era_long_c),squeeze(Rc_c));

colormap(b2r(-c_limit,c_limit));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(name,'2mT')
    color_alt=1;
elseif strcmp(name,'v850') 
    color_alt=3;
elseif strcmp(name,'u850') 
    color_alt=5;
else    
    color_alt=4;    
end


if color_alt==1
    colormap(brewermap(256,'*RdBu')); % mainly used
elseif color_alt==2
    colormap(flipud(cbrewer('div','Spectral',10)));
  
elseif color_alt==3
    colormap(flipud(cbrewer('div','PiYG',20))); 
    
elseif color_alt==4
    colormap(flipud(cbrewer('div','RdYlGn',20))); 
    
elseif color_alt==5
    colormap(flipud(cbrewer('div','PuOr',20)));     
    
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hSurfAx=(gca);
cRange= caxis;

if corr_label==1
    h=colorbar('EastOutside');     
 
    set(h, 'Position', [.76 .145 .015 .75], 'FontSize',18, 'FontWeight', 'bold');  
 

 %%%%%%%%%%%% Move Colorbar  %%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%% move colobar and corr title
 
        if strcmp(name,'2mT') % && (strcmp(season,'annual')) && (strcmp(iso,'d18O'))
            x_move_colorbar=0.065;  % 0.16; 0.11
         
        elseif  strcmp(name,'v850') || strcmp(name,'u850')  
            x_move_colorbar=0.135; 
        else   
          %  x_move_colorbar=0.137;  % 0.16; 0.11
            x_move_colorbar=0.10;
        end
 
 %         if  strcmp(iso,'Accumulation')==1 || strcmp(iso,'SST_1')==1 || strcmp(iso,'SST_2')==1
            pos_c = get(h,'position'); 
            pos_c(1,1) = pos_c(1,1)+x_move_colorbar;    
            pos_c(2)=pos_c(2)+0.14;  
%             pos_c(3)=pos_c(3)- 0.0052;  % widthcolorbar

            pos_c(3)=pos_c(3)+0.014; 

            pos_c(4)=pos_c(4)-0.30;  % height colorbar
            set(h,'pos',pos_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if era_name_nr==1  
        txt_h=axestext_c(1.05,0.580,'Correlation','rotation',-90,'FontSize',20, 'FontWeight', 'bold'); % z500
        
    elseif era_name_nr==3 || era_name_nr==2
         txt_h=axestext_c(1.09,0.580,'Correlation','rotation',-90,'FontSize',20, 'FontWeight', 'bold');
    elseif era_name_nr==5
        txt_h=axestext_c(1.008,0.580,'Correlation','rotation',-90,'FontSize',18, 'FontWeight', 'bold'); % 2mT     
         
        
    else
        txt_h=axestext_c(1.03,0.58,'Correlation','rotation',-90,'FontSize',18, 'FontWeight', 'bold'); 
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contourm function keeps projection defined above
    h1= contourm( double(era_lat),double(era_long),squeeze(P(1,2,:,:))', [0.05],'--k','ShowText','off',...
    'Linecolor',[.1 .1 .1],'LineWidth', 2);

caxis(cRange);  
alpha(1)

%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(iso,'d18O')==1
    iso_str='{\delta}^1^8O';
elseif cor_nr==2     
    iso_str='d-excess';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coast_nr=0; % (1/0) On/Off

if coast_nr==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load coast
%%%%%%%%%%%%%
% to be able to use coastline for continents in combination with ant
% coastline and grounding line from bedmap
in_c=find(lat<-60);
lat_cr=lat;
lon_cr=long;
lat_cr(in_c)=NaN;
lon_cr(in_c)=NaN;

%%%%%%%%%%%%%
color_code=[.4 .4 .4];

 plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
%plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);

addpath C:\PHD\matlab_mapping
bedmap2('gl','LineWidth', 2, 'color',color_code );
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%

%%%%%%%%%%%%%%%%%%
show_title=1; %(1/0/2/3/4) 2- just iso string (3) use(4) season

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_title==1

    hp1=title([site,' ',age_scale_str,' r(',iso_str,', ',name,') ',season,' ERA-Interim ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],...
    'FontSize',20, 'FontWeight', 'bold');    

elseif show_title==2
  
    hp1=title([iso_str],...
    'FontSize',22, 'FontWeight', 'bold'); 

elseif show_title==3

    hp1=title([iso_str,', ',name],...
    'FontSize',30, 'FontWeight', 'bold');   

elseif show_title==4

    hp1=title([iso_str,' ',name,' ',season],...
    'FontSize',24, 'FontWeight', 'bold');   


end


    if strcmp(proj,'stereo')==1

    pos=get(hp1,'Position');
    
%         if strcmp(iso,'Accumulation')==1 || cor_nr==15
%             pos(2)=pos(2)-0.00;
            
        if era_name_nr==1    % z500
            %pos(2)=pos(2)-0.05;%0.08
            pos(2)=pos(2)-0.25;%0.08
        elseif era_name_nr==3 || era_name_nr==2   
           pos(2)=pos(2)-0.08; 
        else
            pos(2)=pos(2)-0.27;%0.08
        end
        
    set(hp1,'Position',pos)

    elseif strcmp(proj,'mec')==1
    pos(2)=pos(2)-0.10;
    set(hp1,'Position',pos)  
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Max corr point  %%%%%%%%%%%%%%%

if xr_lat<=-45 %&&  (strcmp(iso,'dD_1213B')==1 ||  strcmp(iso,'dD')==1)    

 if strcmp(iso,'Accumulation')==1 &&  strcmp(name,'tp')==1 
 max_position_str= [330 175 210 50];
 else
 max_position_str= [160 115 190 50];
 end
   
 if show_max==1
         plotm(double(xr_lat),double(xr_lon),'.','MarkerSize',30,'color',[.4,.6,.1]);  % max corr point

  
         TextBox = uicontrol('style','text');
         set(TextBox,'String',[' r = ',num2str(round(X_r_max*1000)/1000)],...
             'position',max_position_str,'FontWeight','bold','FontSize',22 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    site_coor=star_coord_WA(site); % Site Marker


%         plotm(site_coor(1),site_coor(2),'o','MarkerSize',12,'MarkerFaceColor','none','MarkerEdgeColor',[.4,.6,.1],'LineWidth',3.3); 
%       % plotm(site_coor(1),site_coor(2),'o','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.3,.9],'LineWidth',3.2);    
        % plotm(site_coor(1),site_coor(2),'o','MarkerSize',22,'MarkerFaceColor','none','MarkerEdgeColor',[.4,.6,.1],'LineWidth',2.5); 
        plotm(-66.78, 139.56 ,'.','MarkerSize',36,'MarkerFaceColor','none','MarkerEdgeColor',[.0,.0,.0]); % [.9,.1,.9]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% font size for colorbar 
% placed here because it gets over writen if place earlier 

if corr_label==1
    set(h, 'FontSize',18, 'FontWeight', 'bold');
end


if strcmp(name,'z500')==1 
    plabel on;
    mlabel('MLabelParallel',-20,'PLabelLocation',[-75 -45 -15 0],'PLabelMeridian',100) 

elseif strcmp(name,'2mT')==1 
    plabel on;
    mlabel('MLabelParallel',-37 ,'PLabelLocation',[-75 -60 -45 -30],'PLabelMeridian',100)
    
elseif strcmp(name,'v850')==1 || strcmp(name,'u850')==1
    plabel on;
    mlabel('MLabelParallel',-45 ,'PLabelLocation',[-75 -60 -45 -30],'PLabelMeridian',100)
    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold off

if strcmp(name,'2mT')==1
    letter_pos=[200 760 60 60]; 
elseif strcmp(name,'v850')==1  || strcmp(name,'u850')==1  
    letter_pos=[200 790 60 60];     
else
    letter_pos=[200 760 60 60];  
end

         TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
    
%============================================================
% save figure
if corr_label==0
    filename=[site,'_',iso,'_',name,'_',season,'_ERA_interim_correlation _',...
        num2str(round(date_annual(start_t))),'_',num2str(round(date_annual(end_t))),proj];
elseif corr_label==1 % w. colorbar
    filename=[site,'_',iso,'_',age_scale_str,'_',name,'_',season,'_ERA_interim_correlation _',...
        num2str(round(date_annual(start_t))),'_',num2str(round(date_annual(end_t))),proj,'_crop'];
end

% filedir ='rice_isotopes\annual\';
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);

% save as png 
orient landscape
% increase quality of picture by increase to -r500

if figure_format==1
  export_fig('-eps','-nocrop','-painters', '-depsc','-opengl', '-r100',savefilename_c); % EPS works  func needs ghostscript 
elseif figure_format==2
  
   if crop_nr==1
   export_fig('-png','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG 110
   else
       
  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG
  
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 