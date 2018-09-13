%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cryo_review_SIC_c2.m
% Daniel Emanuelsson 2018
% Matlab 2017a
%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * corrcoef_df.m   UoW Steig
% * cosweight.m     UoW online archive Atmospheric Science
% * annave.m        UoW online archive Atmospheric Science
% * star_coord_WA.m DE
% * fig.m           fileexchange
% * export_fig.m    fileexchange
% * grabit.m        fileexchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Digitazing 
% images are saved to this folder 
% dir='C:\g_drive\review\Cryosphere\';
grabit
%%%%%%%%%%%%%%%%
%%

figure
plot(Data002(:,1),Data002(:,2),'.-r')
hold on
plot(Data001_cryo_review(:,1),Data001_cryo_review(:,2))
grid
[r, p]=corrcoef_df(Data002(:,2),Data001_cryo_review((1:16),2));

figure
plot(Data002_d_excess_cryo_review(:,1),Data002_d_excess_cryo_review(:,2))
hold on
plot(Data001(:,1),Data001(:,2),'.-r')
grid


[r, p]=corrcoef_df(Data001(:,2),Data002_d_excess_cryo_review((1:16),2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all; close all

dir_cryo='C:\g_drive\review\Cryosphere\';


iso_nr=1; % (1) 18O (2) d_excess
save_data=0; % (1/0) on/off


if iso_nr==1
    load([dir_cryo,'Data001_d18O_cryo_review.mat'])
    yrs_c=round(Data001_cryo_review(:,1)); % original time-scale
    iso_data_d18O=Data001_cryo_review(:,2);

    % time=[1998:2014]';
    time=[1997:2014]';


    load([dir_cryo,'d18O_digi_high_res_2013_c.mat']); % nssSO4 and iso peak (timing guidance by iso)
    d18O_digi_2013_new_mean=nanmean(d18O_digi_high_res_2013_c(:,2));
    load([dir_cryo,'d18O_digi_high_res_2014_c.mat']);
    d18O_digi_2014_new_mean=nanmean(d18O_digi_high_res_2014_c(:,2));

    % load([dir_cryo,'d18O_digi_high_res_2011.mat']);
    % d18O_digi_2011_new_mean=nanmean(d18O_digi_high_res_2011(:,2));
    % load([dir_cryo,'d18O_digi_high_res_2010.mat']);
    % d18O_digi_2010_new_mean=nanmean(d18O_digi_high_res_2010(:,2));


    %  check_d18O=(d18O_digi_2013_new_mean+d18O_digi_2014_new_mean)/2;
    % Data001_cryo_review(17,2)-check_d18O;

    % matrix
    MA_TA_d18O=NaN(18,4);
    MA_TA_d18O(:,1)=[1997:2014]';
    MA_TA_d18O((1:16),2)=Data001_cryo_review((1:16),2);

    MA_TA_d18O(17,2)=d18O_digi_2013_new_mean;
    MA_TA_d18O(18,2)=d18O_digi_2014_new_mean;

    % %%%%%%%%%
    % MA_TA_d18O(14,2)=d18O_digi_2010_new_mean;
    % MA_TA_d18O(15,2)=d18O_digi_2011_new_mean;
    % MA_TA_d18O((16:18),2)=Data001_cryo_review((15:17),2);
    %%%%%%%%

    MA_TA_d18O((1:17),3)=yrs_c; % Original age-scale
    MA_TA_d18O((1:17),4)=Data001_cryo_review((1:17),2); % Original data
    % save modified record

    name_c='MA_TA_d18O';

    if save_data==1
        savefilename =[dir_cryo, name_c, '.mat'];
        save(savefilename,'MA_TA_d18O'); 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iso_nr==2 % d-excess
    
    load([dir_cryo,'Data002_d_excess_cryo_review.mat'])   
    yrs_c=round(Data002_d_excess_cryo_review(:,1)); % original time-scale
    iso_data_d_excess=Data002_d_excess_cryo_review(:,2);   
    
    time=[1997:2014]';

    % load([dir_cryo,'d_excess_digi_high_res_2013.mat']);
    % d_excess_digi_2013_new_mean=nanmean(d_excess_digi_high_res_2013(:,2));
    % load([dir_cryo,'d_excess_digi_high_res_2014.mat']);
    % d_excess_digi_2014_new_mean=nanmean(d_excess_digi_high_res_2014(:,2));

    load([dir_cryo,'d_excess_digi_high_res_2013_c.mat']);
    d_excess_digi_2013_new_mean=nanmean(d_excess_digi_high_res_2013_c(:,2));
    load([dir_cryo,'d_excess_digi_high_res_2014_c.mat']);
    d_excess_digi_2014_new_mean=nanmean(d_excess_digi_high_res_2014_c(:,2));


    % check_d18O=(d18O_digi_2013_new_mean+d18O_digi_2014_new_mean)/2;
    % Data001_cryo_review(17,2)-check_d18O;

    % matrix
    MA_TA_d_excess=NaN(18,4);
    MA_TA_d_excess(:,1)=[1997:2014]'; % new time scale

    MA_TA_d_excess((1:16),2)=Data002_d_excess_cryo_review((1:16),2);
    MA_TA_d_excess(17,2)=d_excess_digi_2013_new_mean;
    MA_TA_d_excess(18,2)=d_excess_digi_2014_new_mean;

    MA_TA_d_excess((1:17),3)=yrs_c(1:17); % Original age-scale
    MA_TA_d_excess((1:17),4)=Data002_d_excess_cryo_review((1:17),2); % Original data
    % save modified record


    name_c='MA_TA_d_excess';


    if save_data==1
        savefilename =[dir_cryo, name_c, '.mat'];
        save(savefilename,'MA_TA_d_excess'); 
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HadISST SIC corr
 clear all
 close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
lag_fact=0;%%%%%%%-= !! =-%%%%%%%%%%% shift climate data if not zero

param_nr=2;% (1-SST/ 2 SIC)

age_scale_nr=2; % (1) New (2) Original

% SST_dataset=1; % (1-HadISST, 2-ERSST) 

rcontour_psa=0; %(0)(1) show where PSA pattern overlap interferance (2) shows where PSA2 is significant
coast_nr=1;  % On/Off
    
if param_nr==2

name='SIC';

end
 
%  site='RICE';
 site='TA';
  
%  yr_e=2012; %%%%%%%%% 
  
 str_yr_nr=1; % (1-1979  / 2-1950)%%% 
  
  
%%%%%%%%%%%%%%%%%%%%%
sea_nr=1;
season='annual'; 
%%%%%%%%%%%%%%%%%%%%%%
 
 
  %%%%%%%%%%%%%%%%
  show_title=1; % (1)-long (2)-short (3) use for Interval SIC figure
  

        if  str_yr_nr==1
  
           
        
            if param_nr==1
                letter='a';
%                 RSAS_SIC_box=0;
            elseif param_nr==2
%                 letter='d';
                if age_scale_nr==1
                    %yr_e=2012;
                    yr_e=2014;
                    yr_s=1997;
                    letter='a';
                    age_scale_str='New';
                elseif age_scale_nr==2
                    letter='b';
                    yr_e=2014;
                    yr_s=1998;
                    age_scale_str='Original';
                end
%                 RSAS_SIC_box=1;
            end
        
        
        area_2_box=0 ; % 0/1 for RICE SST corr Box and text for Area 1, 2, and 3
        
 
        
        elseif  str_yr_nr==2
        yr_s=1950;
        %letter='a';
        letter='b';
        area_2_box=0 ; % 0/1 for RICE SST corr Box and text for Area 2
        end
        
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Int_in=0; %%% turn (1/0) on/off, on for Fig 10a, b 
 
 interval_nr=2; %(1/2)%%%%%%%%%%%%%%%%%%% 2,3 in use for both SST and SIC
 

 
 show_colorbar=1;
 figure_format=2; %(1)  EPS (2) PNG
 
 iso_nr=2; % (1) d18O (2) d_excess
 %%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%
 if  iso_nr==1
    iso='d18O';
 elseif iso_nr==2
    iso='d_excess';
 end
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
% 
 
show_maximum_point=0;
 
p_level_rc=0.05; % p level choose 0.1 for figures 0.05 for SST time series
  
if  strcmp(name,'SIC')==1
  proj='stereo';
  lat1=-90;
  lat2=-50;
  box_use=0;  % (0/1)  
  
end
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HadISST data
% HadISST data
% version 1.1 Rayner et al. 2003
%
%           standard_name = 'sea_ice_area_fraction'
%           long_name     = 'Monthly 1 degree resolution sea ice concentration'

addpath C:\PHD\HadISST_data\ % Ienovo


% if SST_dataset==1

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    name_c='HadISST_ice_c.nc';
    HadISST_time=ncread(name_c,'time'); %units         = 'days since 1870-1-1 0:0:0'
    % missing values  missing_value = -1e+30 ( looks more like -1000)\
    %monthly values
 
    HadISST_time_c=HadISST_time+datenum(1870,1,1);
 
    date_vec=datevec(double(HadISST_time_c)); 
    yyyy=date_vec(:,1); 
    mm=date_vec(:,2);
    HadISST_year_num=yyyy+(mm/12-1/12); 

    HadISST_time_bnds =ncread(name_c,'time_bnds');
    HadISST_lat=ncread(name_c,'latitude');
    HadISST_lon=ncread(name_c,'longitude');
    
    
         %   name='ice'; 
        HadISST_ice=ncread(name_c,'sic'); % 360x180x1728 lon x lat x time
        M=HadISST_ice;   
    
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HadISST_count=1705; % that overlaps with RICE record 1979-2011
HadISST_count=find(HadISST_year_num==yr_e+lag_fact)+11; % values in months
% HadISST_start=289; %1894
HadISST_start=find(HadISST_year_num==yr_s+lag_fact); % 1950 hardly any data before this
% HadISST_start=1309; %1979 satelite
% (HadISST_count-HadISST_start)/12

% era_count=1693; 
% div=118;
div=ceil((HadISST_count-HadISST_start)/12); %1979

%%%%%%%%%%%%%%%%%%%%%%%
% % cosweight
% S_lon=17;S_lat=170;%check
% M ( S_lon, S_lat-5,2) 
%  M=cosweight_c(M,HadISST_lat); % UoW function
%  M ( S_lon, S_lat-5,2)

HadISST_year_num(HadISST_start);
HadISST_year_num(HadISST_count);

% HadISST_year_num()-HadISST_year_num()
% HadISST_start:HadISST_count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% ===================================================
 % Added but does not make a differnce for annual corr at least
 
  % seasonality
 
 % help annave
 
 HadISST_lat_c=HadISST_lat;

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

era_z500_c2=cosweight(era_z500_c,HadISST_lat_c); % Original UoW function    (time x lat x lon) 

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

% run from start if year is changed above
% (360x180x1728) (Long,Lat,month)
[m n t]=size(M(:,:,HadISST_start:HadISST_count));% time(months) % 1894-2011
for i= 1:m
    for j= 1:n

        dummy1=reshape(squeeze(M(i,j,HadISST_start:HadISST_count)), 12, div);

   
if strcmp(season,'annual')==1
         HadISST_M_annual_temp=nanmean(dummy1,1); % annual HadISST-values

         
%          HadISST_M_annual_temp=anomaly(HadISST_M_annual_temp); % try
%          normalizing
         HadISST_M_annual_dummy=detrend(HadISST_M_annual_temp(1:end));
         
         HadISST_M_annual(i,j,:)=HadISST_M_annual_temp(1:end);
         HadISST_M_annual_detrend(i,j,:)= HadISST_M_annual_dummy; % detrended
%          
end

%   clear       dummy1 HadISST_M_annual_dummy HadISST_M_annual_temp

    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(iso,'d18O')==1
    
% load('C:\g_drive\review\Cryosphere\Data001_d18O_cryo_review.mat')
load('C:\g_drive\review\Cryosphere\MA_TA_d18O.mat')    
    
% isotopes record
% date_annual=[1900:2011]';
%date_annual=[1998:2014]';
% date_annual=[1997:2014]';

if age_scale_nr==1 % New
    date_annual=MA_TA_d18O(:,1);
    stacked_record_annual_Ma=MA_TA_d18O(:,2);
elseif age_scale_nr==2 % Original  
    date_annual=MA_TA_d18O(:,3);
    stacked_record_annual_Ma=MA_TA_d18O(:,4);
end


start_t=find(date_annual==yr_s);
end_t=find(date_annual==yr_e); 


%     stacked_record_annual_Ma=MA_TA_d18O;
    X=detrend(stacked_record_annual_Ma((start_t:end_t),1));


elseif strcmp(iso,'d_excess')==1
    
%load('C:\g_drive\review\Cryosphere\Data002_d_excess_cryo_review.mat')
load('C:\g_drive\review\Cryosphere\MA_TA_d_excess.mat')
% date_annual=[1997:2014]';
if age_scale_nr==1 % New
    date_annual=MA_TA_d_excess(:,1);
    stacked_record_annual_Ma=MA_TA_d_excess(:,2);
elseif age_scale_nr==2 % Original  
    date_annual=MA_TA_d_excess(:,3);
    stacked_record_annual_Ma=MA_TA_d_excess(:,4);
end

start_t=find(date_annual==yr_s);
end_t=find(date_annual==yr_e);     
    
    

%     date_annual=MA_save(:,1);
    
    X=detrend(stacked_record_annual_Ma((start_t:end_t),1));

       
    
end 
   


if strcmp(season,'annual')==1
    %season='Annual';
    Y=HadISST_M_annual_detrend; 
end


[A B C]=size(Y);
R=zeros(2,2,A,B);P=zeros(2,2,A,B);
for i=1:A
    for j=1:B
        Q=squeeze(Y(i,j,:));
            
        [R(:,:,i,j) P(:,:,i,j)]=corrcoef_df(X,Q, 'rows','pairwise');
        
    end
end
%for the WAIS site:
% % for RICE
% R(:,:,find(floor(HadISST_lon)==199),find(floor(HadISST_lat)==-80))
% P(:,:,find(floor(HadISST_lon)==199),find(floor(HadISST_lat)==-80))
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Figure  
 % axesm defined projection
 
  lock_scalebar=1;%%%%%%%%%%%%%%%%%%%%%%%
    
    RSAS_SIC_box=0; % (1) area defined using RICE corr pattern (2) ADP following Yuan 2004
    alt_pc=0;
 
 %%%%%%%
% if strcmp(proj,'stereo')==1

%     if strcmp(name,'SIC')
    f2=fig('units','inches','width',10.5,'height',10.5,'font','Helvetica','fontsize',16,'border','on'); 


  axesm( 'stereo','MapLatLimit',[lat1 lat2],'Grid','on',...
      'Frame','on','ParallelLabel','on',...
       'MeridianLabel','on','FontWeight','bold','FontSize',18, 'mlabellocation',[0:-30:-360, 0:30:180]); 
 
  set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes 
  
 
 hold on

                set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );
      
 gridm('GLineStyle','--','Gcolor',[.5 .5 .5],'GLineWidth',1.5,... %%%%%%%%%% Grid
    'MLineLimit',[lat1 lat2],'Galtitude', .02)



Rc=R;
Rc5=R;
Rc2=R;

Rc(P>=p_level_rc)=NaN; 
Rc2(P>=0.01)=NaN; % Mask out above or equal R's

Rc5(P>=0.1)=NaN; 




% surfm(double(HadISST_lat),double(HadISST_lon),HadISST_ice(:,:,2)')


c_min=squeeze(R(1,:,:));
c_min=c_min(2,:);
c_min=min(c_min);

c_max=squeeze(R(1,:,:)); 
c_max=c_max(2,:);
c_max=max(c_max);

c_limit=max(abs(c_min),abs(c_max));
c_limit_c=c_limit;


if strcmp(name,'SIC')==1

 X_r=squeeze(Rc(1,2,:,:))';

X_r(:,(100:360))=NaN;  %!!!!!!!!!! Look ONLY in RS ABS region,  nan out around east antarcica 
X_r((1:90),:)=NaN; % southern Hem
 
%  X_r_max=nanmax(nanmax(X_r));   
  X_r_max=nanmin(nanmin(X_r));  
 
X_r_long=nanmin(X_r); % 1 long
xr=find(X_r_long==X_r_max); 
xr_lon=HadISST_lon(xr);

X_r_lat=nanmin(X_r,[],2); %lat
xr=find(X_r_lat==X_r_max); 
xr_lat=HadISST_lat(xr);
 
   
end


 if lock_scalebar==1       
      c_limit=0.8;
     
 end
        mlabel('MLabelParallel',-55  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100)   


%%%%%%%%%%%%%%%% Fill in seam  % only one side
% HadISST_lat_c=HadISST_lat;
HadISST_lon_c=HadISST_lon;
HadISST_lon_c(361,:)=179.99;

Rc6=squeeze(Rc5(1,2,:,:))';
Rc6(:,361)=Rc6(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SST
hSurf=surfm(double(HadISST_lat),double(HadISST_lon_c),squeeze(Rc6));
colormap(b2r(-c_limit,c_limit));
%colormap(jet);
% colormap(brewermap(256,'*RdBu'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

if Int_in==1
   color_alt=1;   
else
    color_alt=6;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end


if color_alt==1
    colormap(brewermap(256,'*RdBu'));
elseif color_alt==2
    colormap(flipud(cbrewer('div','Spectral',10)));
elseif color_alt==3
    colormap(flipud(cbrewer('div','PiYG',20)));      
elseif color_alt==4
    colormap(flipud(cbrewer('div','RdGy',20)));      
elseif color_alt==5
    colormap(flipud(cbrewer('div','RdBu',20)));    
elseif color_alt==6
    colormap(cbrewer('div','BrBG',20));
elseif color_alt==7
    colormap(flipud(cbrewer('div','RdYlBu',20)));    
end

%%%%%%%%%%%%%%%%%%%%%%


%  hSurf=surfm(double(HadISST_lat),double(HadISST_lon),squeeze(Rc(1,2,:,:))'); colormap(flipud(jet)); % for sea ice

hSurfAx=(gca);
cRange= caxis;

% contourm function keeps projection defined above
    h1= contourm( double(HadISST_lat),double(HadISST_lon),squeeze(P(1,2,:,:))', [0.05],'-k','ShowText','off',...
   'LineWidth', 2);


% Then reset the color axis according to the range determined above:
caxis(cRange);
%  
alpha(1)

 
%%%%%%%%%%%%%%%%%%%%%%%%
if RSAS_SIC_box==1    
        
        lon_b_c=[-120 -120 -150 -150 -120];
        lat_b_c=[ -73 -60 -60 -73 -73];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .7 .1]) 
        
 elseif RSAS_SIC_box==2         % ADP
        
        lon_b_c=[-130 -130 -150 -150 -130];
        lat_b_c=[ -70 -60 -60 -70 -70];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .2 .2])  % 1 .7 .1
        
        lon_b_c=[-40 -40 -60 -60 -40];
        lat_b_c=[ -70 -60 -60 -70 -70];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .2 .2])         
        
end        
  
%%%%%%%%%%%%%%%%%%%%%%%%

%c_limit_c
contour_lim=1.3;
max_contour=1; % (1/0)
polygon_nr=2; % (1/2)
contour_p_level=0.05;

if max_contour==1
contour_lat_lon=contourm( double(HadISST_lat),double(HadISST_lon),squeeze(P(1,2,:,:))',...
       [contour_p_level contour_p_level],'-k','ShowText','off',...
    'Linecolor',[.1 .1 .1],'LineWidth', 2);%,'LineColor','none');  %%%%%%%%%%%%%%%%%


end
%%%%%%%%%%%%%%%%%%%%%%%%%


% ind_lon = find(HadISST_lon<=-93.5 & HadISST_lon>=-104);
% ind_lat =find(HadISST_lat<=-33 & HadISST_lat>=-38);

hSurfAx=(gca);
cRange= caxis;


if show_colorbar==1
  %h=colorbar('SouthOutside'); 
  h=colorbar('EastOutside');  
 % h=colorbar2('vert'); % from UoW atmo page 
% h=colorbar;


  
end

if  param_nr==1
       a1=axestext_c(1.093,+0.67, ['Correlation'] ,'rotation',-90,'FontSize',18, 'FontWeight', 'bold');
elseif    param_nr==2
    
    a1=axestext_c(1.065,+0.57, ['Correlation'] ,'rotation',-90,'FontSize',20, 'FontWeight', 'bold');    %050
    
        pos=get(h,'Position');
%         pos(1)=pos(1)+0.024;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
        pos(1)=pos(1)+0.028; % 0.035
        pos(2)=pos(2)+0.14;  
        pos(3)=pos(3)- 0.0052;  % widthcolorbar
        pos(4)=pos(4)-0.31;  % height colorbar
        
        set(h,'Position',pos)
end



if strcmp(iso,'dD')==1 || iso_nr==1
    iso_str='{\delta}^1^8O';
elseif strcmp(iso,'d_excess')==1 || iso_nr==2   
    iso_str='d-excess';
end

 

TextBox = uicontrol('style','text');
if  strcmp(proj,'merc')==1          
            
  letter_position=  [180 500 40 50];% [200 635 50 50];
  

    %set(TextBox,'String',letter,'position',[130 820 50 50],'FontWeight','bold','FontSize',28 ); % x position y position xsize ysize  
    set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',28 ); % x position y position xsize ysize            
            
elseif  strcmp(proj,'stereo')==1
    
       if strcmp(name,'SIC')==1 
       letter_position=  [200 820 60 60];%
       else
       letter_position=   [330 840 50 50];
        end
    
            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize 
            %285
end
            set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coast_nr=0;  % On/Off

if coast_nr==1
    
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

addpath C:\PHD\matlab_mapping

if strcmp(name,'SIC')==1
color_code=[.4 .4 .4];    
    plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
    bedmap2('gl','LineWidth', 2,'color',color_code);
elseif strcmp(name,'SST')==1
    
    if rcontour_psa==0
    color_code=[.4 .4 .4];
    elseif rcontour_psa==1
     color_code=[.1 .1 .1];   
    end
    
    plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
%     plot3m(lat_cr,lon_cr,'-k','LineWidth', 1.5);   
    bedmap2('gl','k','LineWidth', 1.5);
end 

addpath C:\PHD\matlab_mapping
bedmap2('gl','LineWidth', 2, 'color',color_code );
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if strcmp(name,'SIC')==1 && show_maximum_point==1     
 
 plotm(double(xr_lat),double(xr_lon),'.','color',[.4,.8,.2],'MarkerSize',36);  % max corr point
 

         TextBox = uicontrol('style','text');
         set(TextBox,'String',[' r = ',num2str(round(X_r_max*1000)/1000)],...
             'position',[310 130 170 50],'FontWeight','bold','FontSize',22 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);    
end


if  strcmp(site,'TA')==1
%         site_coor=star_coord_WA(site); % RICE marker
        % Circle Alt
        plotm(-66.78,139.56,'.','MarkerSize',30,'MarkerFaceColor','none','MarkerEdgeColor',[.0,.0,.0]); %,'LineWidth',2.5); % [.4,.6,.1] [.9,.1,.9]
        %%%% Dot Alt
        %plotm(site_coor(1),site_coor(2),'.','MarkerSize',40,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.1,.9],'LineWidth',2.5); 
        
   
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_title==1 || show_title==2 || show_title==3
      
  SST_dataset_str='HadISST';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if show_title==1

                        hp1=  title([site,' ',age_scale_str,' r(',iso_str,', ',name,') ',SST_dataset_str,' ',season,...
                            ' ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],'FontSize',20, 'FontWeight', 'bold');
        elseif show_title==2
                        hp1=  title([iso_str,', ',name],'FontSize',28, 'FontWeight', 'bold');
                                               
        elseif show_title==3
                        hp1=  title([iso_str,' ',name,' ',SST_dataset_str,' ',...
                            ' ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],'FontSize',20, 'FontWeight', 'bold');
             
        end

        if strcmp(proj,'stereo')==1  
            pos=get(hp1,'Position');
            pos(2)=pos(2)-0.05;
            set(hp1,'Position',pos,'FontSize',28, 'FontWeight', 'bold');    
  
        elseif strcmp(proj,'mec')==1

            pos=get(hp1,'Position');
            pos(2)=pos(2)-0.02;
        	set(hp1,'Position',pos,'FontSize',18, 'FontWeight', 'bold');    

        end
       
end

    if show_colorbar==1
    set(h, 'FontSize',16, 'FontWeight', 'bold'); 
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off
%%%%%%%%%%%%%%%%%%%%
% save figure

 filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
 filename=strcat(filedir,'HadISST_',site,'_',age_scale_str,'_',name,...
     '_',iso,'_',season,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e));

% Save as Matlab fig
% saveas(f2,[filename,iso,'_HadISST_corr_',season,'_',num2str(round(date_annual(start_t))),...
%     '_',num2str(round(date_annual(end_t)))],'fig')

savefilename_c=filename;

 % save as png 
orient landscape
% increase quality of picture by increase to -r500
if figure_format==1
  export_fig('-eps','-nocrop','-painters', '-depsc','-opengl', '-r100',savefilename_c); % EPS works  func needs ghostscript 
elseif figure_format==2
  export_fig('-png','-painters', '-depsc','-nocrop','-opengl', '-r170',savefilename_c); % PNG  110
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%