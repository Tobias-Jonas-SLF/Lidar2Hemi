# Lidar2Hemi
  Creating synthetic hemispherical images from LiDAR point cloud data

% DESCRIPTION

  This Matlab scripts creates synthethic hemispherical images from 
  high-resolution LiDAR point cloud data of forest canopy. The script
  relates to research published as
  
    D. Moeser, J. Roubinek, P. Schleppi, F. Morsdorf, T. Jonas;
    Canopy closure, LAI and radiation transfer from airborne LiDAR
    synthetic images; 2014; Agricultural and Forest Meteorology;
    doi: 10.1016/j.agrformet.2014.06.008
  
  The script comes along with a demonstration dataset which is provided
  for the purpose of testing the script's functionality. Adaptations
  might become necessary as you apply this script to other datasets.
  Comments in the script may guide you through possible adaptations.
   
  Note further that this script has been developed for aerial LiDAR 
  data flown with a point density of approx 36pt/m2. Applications to
  LiDAR data from terrestrial scanning and/or with significantly
  different point densities may require other methodology.

% SETUP REQUIRED TO RUN THIS SCRIPT

  Matlab base version 7.0 or higher. As far as we are aware of no 
  additional toolbox is needed. Add parent folder Lidar2Hemi
  to the Matlab path, or run script from the very folder.
  
% GETTING STARTED

  1) copy entire folder/file system to a local directory 'Lidar2Hemi'
  2) open Lidar2Hemi.m in your Matlab editor, specify the full path 
     of your local directory as the base folder in line 98, then save
     Lidar2Hemi.m
  3) set Matlab path to include the base folder and all subfolders
  4) unzip demo las files within folder DSM_Data, adding the two
     files 'dsm_laret_low.las' and 'ndsm_laret_low.las'
  5) run Lidar2Hemi from the command line 

% IMPLEMENTATION

  by Tobias Jonas with input from David Moeser*, Clare Webster*,  
  Gulia Mazzotti*, Johanna Malle*, and Felix Morsdorf**
  
  *WSL Institute for Snow and Avalanche Research SLF
    Davos, Switzerland
     
  **University of Zurich, Department of Geography
    Zurich, Switzerland 
  
% VERSION / LAST CHANGES

  v2.6 / 2019-06-24 / by TJ
  
% DATA REQUIREMENTS

  DSM/DTM/DEM. For this script we need all three; the DSM representing
    the surfsce including trees and the DTM representing the underlying
    terrain, both LiDAR-based point-cloud data of the area of interest;
    and the DEM describing the surrounding topography in a gridded 
    representation.
    
  DSM data. Expected in the standard LAS format, or the Matlab 
    equivalent of las files, i.e. the LASM format. The LAS format is a 
    public file format for the interchange of 3-dimensional point cloud
    data. The LAS format is defined by the American Society for 
    Photogrammetry and Remote Sensing (ASPRS). Please refer to their 
    webpage for more information http://www.asprs.org/Committee-General
    .../LASer-LAS-File-Format-Exchange-Activities.html. 
    All data parsed through this script as well as the demonstration
    dataset use LAS 1.4. Compatibility with other versions of LAS
    cannot be guaranteed. The LAS reader and affiliated files
    utilized here are from Dr. Felix Morsdorf within the Remote Sensing
    Laboratories of the University of Zürich. Of note, it is possible
    to use other LAS readers as long as the output remains uniform to
    the below example and there is access to the X, Y, Z data. Note
    that this scripts can also handle nDSM data instead of DSM data,
    which is a differential representation of DSM where nDSM = DSM-DTM
    
  DTM data. Is either expected in LAS format, in LASM format, or can 
    alternatively be provided as a list of point data. This list must
    consist of 3 columns for easting, northing, elevation and be saved 
    in a text format. Note, the DTM data should not contain 
    discontinous data (e.g. canopy data misclassified as ground 
    returns). If necessary, filter/smooth DTM data prior to using it 
    here as input.
    
  DEM data. Is expected in standard ascii grid format, see the example
    in the demonstration dataset for more information.
    
  Evaluation points
    Evaluation points are coordinates for which synthethic
    hemispherical images are to be calculated. These points are defined
    from a list of point data. This list must consist of 2 columns
    for easting and northing and be saved in a text format.
    
  Note: all above datasets must feature data in the same coordinate 
    system (true coordinates). The demonstration dataset features data 
    in the SwissGrid coordinate syste CH1903/LV03.
    
% OUTPUT

  this script will generate a synthethic hemispherical image per
  coordinate. The images are saved according to type/format/color
  settings specified in the user settings below. 
  
% USAGE

  adapt user settings as necessary
  run Lidar2Hemi from the command line of Matlab
  
% ADVISE

  this script uses screenshots to save images, be advised to not 
  interfere with the computer (run other applications, use the mouse, 
  click to other windows, etc.) while excecuting this script, even
  though faulty screenshots occur only rarely
