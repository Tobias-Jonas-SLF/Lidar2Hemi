function Lidar2Hemi

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% DESCRIPTION
  %   This scripts creates synthethic hemispherical images from high
  %   resolution LiDAR point cloud data of forest canopy. The script
  %   relates to research published as
  %
  %     D. Moeser, J. Roubinek, P. Schleppi, F. Morsdorf, T. Jonas;
  %     Canopy closure, LAI and radiation transfer from airborne LiDAR
  %     synthetic images; 2014; Agricultural and Forest Meteorology;
  %     doi: 10.1016/j.agrformet.2014.06.008
  %
  %   The script comes along with a demonstration dataset which is provided
  %   for the purpose of testing the script's functionality. Adaptations
  %   might become necessary as you apply this script to other datasets.
  %   Comments in the script may guide you through possible adaptations.
  % 
  %   Note further that this script has been developed for aerial LiDAR 
  %   data flown with a point density of approx 36pt/m2. Applications to
  %   LiDAR data from terrestrial scanning and/or with significantly
  %   different point densities may require other methodology.
  %
  % SETUP REQUIRED TO RUN THIS SCRIPT
  %   Matlab base version 7.0 or higher. As far as we are aware of no 
  %   additional toolbox is needed. Add parent folder Lidar2Hemi
  %   to the Matlab path, or run script from the very folder.
  % 
  % IMPLEMENTATION
  %   by Tobias Jonas with input from David Moeser*, Clare Webster*,  
  %   Gulia Mazzotti*, Johanna Malle*, and Felix Morsdorf**
  %   * WSL Institute for Snow and Avalanche Research SLF
  %     Davos, Switzerland 
  %   **University of Zurich, Department of Geography
  %     Zurich, Switzerland 
  %
  % VERSION / LAST CHANGES
  %   v2.6 / 2019-06-24 / by TJ
  %
  % DATA REQUIREMENTS
  %   DSM/DTM/DEM. For this script we need all three; the DSM representing
  %     the surfsce including trees and the DTM representing the underlying
  %     terrain, both LiDAR-based point-cloud data of the area of interest;
  %     and the DEM describing the surrounding topography in a gridded 
  %     representation.
  %   DSM data. Expected in the standard LAS format, or the Matlab 
  %     equivalent of las files, i.e. the LASM format. The LAS format is a 
  %     public file format for the interchange of 3-dimensional point cloud
  %     data. The LAS format is defined by the American Society for 
  %     Photogrammetry and Remote Sensing (ASPRS). Please refer to their 
  %     webpage for more information http://www.asprs.org/Committee-General
  %     .../LASer-LAS-File-Format-Exchange-Activities.html. 
  %     All data parsed through this script as well as the demonstration
  %     dataset use LAS 1.4. Compatibility with other versions of LAS
  %     cannot be guaranteed. The LAS reader and affiliated files
  %     utilized here are from Dr. Felix Morsdorf within the Remote Sensing
  %     Laboratories of the University of Zürich. Of note, it is possible
  %     to use other LAS readers as long as the output remains uniform to
  %     the below example and there is access to the X, Y, Z data. Note
  %     that this scripts can also handle nDSM data instead of DSM data,
  %     which is a differential representation of DSM where nDSM = DSM-DTM
  %   DTM data. Is either expected in LAS format, in LASM format, or can 
  %     alternatively be provided as a list of point data. This list must
  %     consist of 3 columns for easting, northing, elevation and be saved 
  %     in a text format. Note, the DTM data should not contain 
  %     discontinous data (e.g. canopy data misclassified as ground 
  %     returns). If necessary, filter/smooth DTM data prior to using it 
  %     here as input.
  %   DEM data. Is expected in standard ascii grid format, see the example
  %     in the demonstration dataset for more information.
  %   Evaluation points
  %     Evaluation points are coordinates for which synthethic
  %     hemispherical images are to be calculated. These points are defined
  %     from a list of point data. This list must consist of 2 columns
  %     for easting and northing and be saved in a text format.
  %   Note: all above datasets must feature data in the same coordinate 
  %     system (true coordinates). The demonstration dataset features data 
  %     in the SwissGrid coordinate syste CH1903/LV03.
  %
  % OUTPUT
  %   this script will generate a synthethic hemispherical image per
  %   coordinate. The images are saved according to type/format/color
  %   settings specified in the user settings below. 
  %
  % USAGE
  %   adapt user settings as necessary
  %   run Lidar2Hemi from the command line of Matlab
  %
  % ADVISE
  %   this script uses screenshots to save images, be advised to not 
  %   interfere with the computer (run other applications, use the mouse, 
  %   click to other windows, etc.) while excecuting this script, even
  %   though faulty screenshots occur only rarely
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% USER SETTINGS
  basefolder          = '...\Lidar2Hemi';
  % path to folder which contains this script as well as all subfolders 
  % to auxiliary functions / data required to run this script 
  
  DSM                 = 'DSM_Data\dsm_laret_low.las';
  % relative path\file to data file that contains DSM/nDSM data. The
  % absolute path to the DSM data file is derived as 
  % fullfile(basefolder,DSM)
  
  DTM                 = 'DTM_Data\dtm_laret_low.las';
  % relative path\file to data file that contains DTM data. The absolute
  % path to the DTM data file is derived as fullfile(basefolder,DTM)
  
  DEM                 = 'DEM_Data\dem_100m.txt';
  % relative path\file to data file that contains DEM data. The absolute
  % path to the DEM data file is derived as fullfile(basefolder,DEM)
  
  PTS                 = 'PTS_Data\pts_laret_low_5.txt';
  % relative path\file to text file that contains coordinates for which
  % synthethic hemispherical images are to be calculated. The absolute path
  % to the text file is derived as fullfile(basefolder,points)
  
  output_path         = 'Output_Images';
  % relative path to folder into which output images will be saved. The
  % absolute path to the output folder is derived as 
  % fullfile(basefolder,output_path)
  
  eval_peri           = 100;
  % this parameter defines a perimeter around every evaluation point within
  % which DTM/DSM is considered when calculating a synthethic hemispherical
  % image. The evaluation perimeter is defined as a diameter around a point
  % where the diameter is to be given in the same units as the DTM/DSM data
  % (m in case of the demonstration dataset). Small evaluation perimeter
  % will allow for quick calculation, but decrease the accuracy of the 
  % hemispherical images at high zenith angles (close to to the horizon). 
  % Large evaluation perimeter on the other hand will allow more accurate 
  % represenation of far distance canopy elements in the hemispherical 
  % images but require more CPU time.
  
  topo_peri           = 30000;
  % this parameter defines a perimeter around every evaluation point within
  % which DEM is considered when calculating topographic shading for a 
  % synthethic hemispherical image. The topography perimeter is defined as 
  % a radius around a point where the radius is to be given in the same
  % units as the DTM/DSM data (m in case of the demonstration dataset). 
  % Small topography perimeters will allow for quick calculation, but may
  % decrease the accuracy of terrain shading from far distance topography.
  % Large topography perimeter on the other hand will allow more accurate 
  % represenation of far distance topographic elements in the hemispherical 
  % images but require more CPU time.
  
  hg_cutoff           = 1.25;
  % this parameter defines a height of canopy elements above the forest
  % floor underneath which all LiDAR points are neglected. A height cutoff
  % is necessary to a) simulate a hemispherical image taken from a certain
  % height above the forst floor, and b) to avoid unrealiable output images
  % due to insufficient point cloud densities as you approach the forest
  % floor (assuming LiDAR data to be taken from above the canopy). 
  
  marker_size         = [7 0.5];
  % this parameter is an import metric that defines the plotting size of  
  % LiDAR points as a function of distance to the camera position. In the 
  % synthetic hemispherical image, a far canopy element should be plotted
  % as a small item, whereas a near canopy element should be represented as
  % a large item. The appearance of the resulting image is significantly 
  % influenced by this parameter. The first number is the plotting size of
  % a single LiDAR point at zero distance to the camera, the second number 
  % is the plotting size of a single LiDAR point at maximum distance to the
  % camera (i.e. eval_peri). Between these two values the plotting size is 
  % a linear function of the distance to the camera. Note: an optimal value
  % for this parameter may depend on the point density of the LiDAR data
  % used as well as the output image size. We suggest to adapt this 
  % parameter to allow calculation of realistic LAI values when inputing
  % the resulting synthethic hemispherical images into respective software
  % packages like Hemisfer (cf. figure 4 of Moeser et al., 2014)
  
  radius              = 350;
  % this parameter relates to the output size of the resulting synthethic
  % hemispherical images. The value defines the radius of the 180° image 
  % circle in pixels. Make sure your display has a height of at least 2.5 
  % times the radius.
  
  camera_height       = 1.5;
  % this parameter sets the virtual camera above ground by camera_height
  
  hidegrid            = 0;                                         
  % default value is 1. Set to 0 in order to display a polar coordinate
  % grid on top of the resulting hemispherical image. Note, hidegrid must 
  % be enabled (i.e. hidegrid = 1) for further processing of the images,
  % this includes both, subsequent external processing with hemisfer and
  % internal calculation of SVF and LAI 
  
  image_type          = 1;                                         
  % set to 1 for output image to include DSM/DTM/DEM
  % set to 2 for output image to include DSM/DTM
  % set to 3 for output image to include DSM
  
  image_color         = 1;                                         
  % set to 1 for output color image (DEM = gray; DTM = brown; DSM = green)
  % set to 0 for output black/white image
  
  image_format        = 'png';                                       
  % format for output of hemispherical images
  % available options are: 'gif','png','tif','jpg'
  % note, gif images are available in black/white format only
  % note, jpg is a lossy format, png and tif are preferred formats
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CALCULATIONS
  %% > preparatory steps
  
  %%% add path to aux functions 
  addpath([basefolder '\Aux_Functions']);
  
  %%% gather information about computer screen used
  scrsz               = get(0,'ScreenSize');
  screenratio         = scrsz(3) / scrsz(4);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > read DTM/DSM/DEM/point data
  
  %%% read coordinate of points for which synthethic hemispherical images are to be calculated
  ev_coors            = dlmread(fullfile(basefolder,PTS),'');
 
  %%% read DSM/nDSM data
  [path,name,ext]     = fileparts(DSM);
  if strcmpi(ext,'.las')                                                   % read DSM/nDSM data assuming las format
    [dsm,hdr]           = readlas(fullfile(basefolder,DSM));               
  elseif strcmpi(ext,'.lasm')                                              % read DSM/nDSM data assuming lasm format
    dsm = load(fullfile(basefolder,DSM),'-mat');
    fnames = fieldnames(dsm);
    dsm = dsm.(fnames{1});
  else
    error('unknown input file format for dsm data')
  end;
  
  %%% read DTM data                                                        
  [path,name,ext]     = fileparts(DTM);
  if strcmpi(ext,'.las')                                                   % read DTM data assuming las format
    [dtm,hdr]         = readlas(fullfile(basefolder,DTM));                 % reads point cloud data (results in a structure array)
  elseif strcmpi(ext,'.lasm')                                              % read DTM data assuming lasm format
    dtm = load(fullfile(basefolder,DSM),'-mat');                           % reads point cloud data (saved as a structure array)
    fnames = fieldnames(dtm);
    dtm = dtm.(fnames{1});
  else                                                                     % read DTM data assuming text format
    xyz               = dlmread(fullfile(basefolder,DTM));                 % reads xyz point list data and convert to structure array 
    dtm.x             = xyz(:,1);
    dtm.y             = xyz(:,2);
    dtm.z             = xyz(:,3);
  end

  %%% read DEM data
  dem                 = load_ascii_grid(fullfile(basefolder,DEM));         % read DEM data assuming ascii/text format
  xdem                = [dem.xllcorner:dem.cellsize:dem.xllcorner+dem.cellsize*(dem.ncols-1)]+dem.cellsize/2 ; % convert reference from lower left corner of each cell to its center of cell
  ydem                = [dem.yllcorner:dem.cellsize:dem.yllcorner+dem.cellsize*(dem.nrows-1)]+dem.cellsize/2'; % convert reference from lower left corner of each cell to its center of cell
  [dem.x dem.y]       = meshgrid(xdem,ydem);
  dem.z               = dem.data;                                          % unnecessary,  but copied into field z for consistency reasons
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > preprocess DSM/DTM data

  %%% interpolate DTM to DSM coordinates (point to point)
  tsi = TriScatteredInterp(dtm.x,dtm.y,dtm.z);
  dsm.e = tsi(dsm.x,dsm.y);
  
  %%% convert DSM to nDSM if not done already
  if mode(dsm.z) > 100
    dsm.z = dsm.z - dsm.e;                                                 % height about terrain
  end

  %%% ignore values below height cutoff hth
  dsm.z(dsm.z <= hg_cutoff) = NaN;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > create data for synthetic hemispherical image
  
  for coorix = 1:length(ev_coors(:,1));                                    % loop through evaluation coordinates  
    pts.x(coorix)     = ev_coors(coorix,1);
    pts.y(coorix)     = ev_coors(coorix,2);
    
    %%% create empty output figure which is the size of your screen
    fh = figure('position', [0 0 scrsz(3)  scrsz(4)]);
    figpos = get(fh,'position');
    
    %%% create empty axis with correct size to results in a hemispherical image of radius as specified 
    ah = axes;
    midpoint = [((((scrsz(4)-figpos(4))/2)*screenratio)+figpos(3))/2 (((scrsz(4)-figpos(4))/2)+figpos(4))/2];
    set(ah,'units','pixel','position',[midpoint(1)-radius, midpoint(2)-radius*1.15, 2*radius, 2*radius*1.15]);
        
    %%% transfer the global coordinate system into a project ccoordinate system relative to the point in question   
    dsm.xpc           = dsm.x - pts.x(coorix);
    dsm.ypc           = dsm.y - pts.y(coorix);
    dsm.zpc           = dsm.z + dsm.e - tsi(pts.x(coorix),pts.y(coorix)) - camera_height;
    
    %%% convert data to a spherical system
        % where theta (.tht) is the counterclockwise angle in the xy plane measured from the positive x axis
        % where phi (.phi) is the elevation angle from the xy plane
        % where radius (.rad) is the 3d distance from the center coordinate
    [dsm.tht,dsm.phi,dsm.rad] = cart2sph(dsm.xpc(:),dsm.ypc(:),dsm.zpc(:));
    
    %%% flip theta (convert viewing perspective from top-down to bottom-up) 
    dsm.tht           = - dsm.tht;
    
    %%% convert elevation angle phi into zenith angle
    dsm.phi           = (pi/2) - dsm.phi;
      
    %%% convert unit of zenith angle to degree
    dsm.phi           = dsm.phi * (180/pi);
    
    %%% eliminate all zenith angle values greater than viewing angle (90 degrees)
    dsm.phi(dsm.phi >= 90) = NaN;
      
    %%% crop DSM data outside of a horizontal radius as specified (eval_peri)
    keepix            = find(dsm.rad.*sin(dsm.phi) < eval_peri);
    dsm.tht           = dsm.tht(keepix);
    dsm.phi           = dsm.phi(keepix);
    dsm.rad           = dsm.rad(keepix);
    
    if image_type < 3 % to further include DTM data
    
      %%% transfer DTM to polar coordinates in an equivalent manner
      dtm.xpc           = dtm.x - pts.x(coorix);
      dtm.ypc           = dtm.y - pts.y(coorix);
      dtm.zpc           = dtm.z - tsi(pts.x(coorix),pts.y(coorix)) - camera_height;
      [dtm.tht,dtm.phi,dtm.rad] = cart2sph(dtm.xpc(:),dtm.ypc(:),dtm.zpc(:));
      dtm.tht           = - dtm.tht;
      dtm.phi           = (pi/2) - dtm.phi;
      dtm.phi           = dtm.phi * (180/pi);
      dtm.phi(dtm.phi >= 90,1) = NaN;

      %%% crop DTM data outside of a horizontal radius as specified (eval_peri)
      keepix            = find(dtm.rad.*sin(dtm.phi) < eval_peri);             
      dtm.tht           = dtm.tht(keepix);
      dtm.phi           = dtm.phi(keepix);
      dtm.rad           = dtm.rad(keepix);
      
      if image_type < 2 % to further include DEM data
    
        %%% transfer DEM to polar coordinates in an equivalent manner
        dem.xpc           = dem.x - pts.x(coorix);
        dem.ypc           = dem.y - pts.y(coorix);
        dem.zpc           = dem.z - interp2(dem.x,dem.y,dem.z,pts.x(coorix),pts.y(coorix)) - camera_height;
        [dem.tht,dem.phi,dem.rad] = cart2sph(dem.xpc(:),dem.ypc(:),dem.zpc(:));
        dem.tht           = -dem.tht;    
        dem.phi           = (pi/2) - dem.phi;
        dem.phi           = dem.phi * (180/pi);
        dem.phi(dem.phi >= 90,1) = 90; % do NOT delete those values but set them to phi = 90 to allow consistent interpolation of horizont line in line 363

        %%% crop DEM data outside of a horizontal radius as specified (topo_peri)
        keepix            = find(dem.rad.*sin(dem.phi) < topo_peri);       % crop data further away than 30km
        dem.tht           = dem.tht(keepix);
        dem.phi           = dem.phi(keepix);
        dem.rad           = dem.rad(keepix);

        %%% calculate horizon lines
        rbins             = [eval_peri:sqrt(2)*dem.cellsize:topo_peri];
        tbins             = [-pi+pi/360:pi/720:pi-pi/360]';
        minphi            = ones(length(tbins),1).*90;                     % initilization of horizont line
        for rbix = length(rbins)-1:-1:1
          fix1 = find(dem.rad >= rbins(rbix) & dem.rad < rbins(rbix+1));
          [spp,fix2] = sort(dem.tht(fix1));
          minphi = min(minphi,interp1q(dem.tht(fix1(fix2)),dem.phi(fix1(fix2)),tbins));
        end;
      end;
    end;
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > plot synthetic hemispherical image
      
    if image_type < 2 % to include DEM data
      
      %%% plot terrain first
      ph = polar_special(ah,[tbins; tbins(1)],[min(minphi,90); min(minphi(1),90)],1,[0 0 0],hidegrid);
      hold on;
      set(ph,'marker','none','linestyle','-','linewidth',2,'color',[0.60 0.60 0.60]);
      for lix = 0.25:0.25:min(minphi)
        ph = polar_special(ah,[tbins; tbins(1)],[min(minphi+lix,90); min(minphi(1)+lix,90)]);
        set(ph,'marker','none','linestyle','-','linewidth',2,'color',[0.60 0.60 0.60]);
      end
    end;

    if image_type < 3 % to include DTM data
        
      %%% plot local terrain second
      ph = polar_special(ah,dtm.tht,dtm.phi,5,[0.50 0.27 0.07],hidegrid);
      hold on;
    end;
    
    %%% plot DSM data last in case a constant plotting size of DSM points is selected (marker_size)
    if length(unique(marker_size)) == 1 
      
      %%% plot all points at once
      ph = polar_special(ah,dsm.tht,dsm.phi,marker_size(1),[0.00 0.30 0.00],hidegrid);
      hold on;
        
    %%% plot DSM data last in case a distance-dependant plotting size of DSM points is selected (marker_size)
    else
      
      %%% allow 10 categories of marker sizes 
      no_categories = 10;                                                  % 10 provides a good balance between plotting accuracy / time, however this value can be increased if needed
      sbins = 0:eval_peri/no_categories:eval_peri;

      %%% loop through marker size category
      for sbix = 1:length(sbins)-1
        
        %%% identify data in distance categories according to marker size categories
        fix = find(dsm.rad >= sbins(sbix) & dsm.rad < sbins(sbix+1));      % select according to 3-d radius (i.e. spherical polar coordinates)
    
        %%% plot only points identifies above
        ph = polar_special(ah,dsm.tht(fix),dsm.phi(fix),(marker_size(1)+(diff(marker_size))*(sbix-1)/(no_categories-1)),[0.00 0.30 0.00],hidegrid);
        hold on;
        
      end;
    end;
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > save synthetic hemispherical image
    
    %%% capture resulting image by recording a screenshot
    set(fh,'color',[1 1 1]);
    if hidegrid
      scrshot = getframe(fh,[midpoint(1)-radius-1, midpoint(2)-radius-1, 2*radius+1, 2*radius+1]);
    else
      scrshot = getframe(fh,[midpoint(1)-radius-50, midpoint(2)-radius-50, 2*radius+100, 2*radius+100]);
    end;
    
    %%% convert to black/white for image output if requested
    if ~image_color
      scrshot.cdata(:,:,1) = double(scrshot.cdata(:,:,3) == 255).*255;
      scrshot.cdata(:,:,2) = double(scrshot.cdata(:,:,3) == 255).*255;
      scrshot.cdata(:,:,3) = double(scrshot.cdata(:,:,3) == 255).*255;
    end;
      
    %%% save resulting image
    switch image_format
    case 'gif'
      imwrite(scrshot.cdata(:,:,3),fullfile(basefolder,output_path,['SHI_' sprintf('%6.0f',pts.x(coorix)) '_' sprintf('%6.0f',pts.y(coorix)) '.gif']),'gif');
    case 'png'
      imwrite(scrshot.cdata(:,:,:),fullfile(basefolder,output_path,['SHI_' sprintf('%6.0f',pts.x(coorix)) '_' sprintf('%6.0f',pts.y(coorix)) '.png']),'png');
    case 'tif'
      imwrite(scrshot.cdata(:,:,:),fullfile(basefolder,output_path,['SHI_' sprintf('%6.0f',pts.x(coorix)) '_' sprintf('%6.0f',pts.y(coorix)) '.tif']),'tif');
    case 'jpg'
      imwrite(scrshot.cdata(:,:,:),fullfile(basefolder,output_path,['SHI_' sprintf('%6.0f',pts.x(coorix)) '_' sprintf('%6.0f',pts.y(coorix)) '.jpg']),'jpg');
    otherwise
      error('unknown output option')
    end
    
    %%% close figure
    close(fh)                                                              % remove this statement if matlab figures shall remain open 
      
  end
  
end