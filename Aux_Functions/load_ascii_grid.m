function answer = load_ascii_grid(varargin)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DESCRIPTION:
  %   copy of function loadgrid within MetDataWizard for use outside of MDW
  %
  % IMPLEMENTATION:
  %   by TJ in November-2014 @ SLF Switzerland
  %   last changes 04.11.2014
  %
  % INPUT:
  %   prompts the user to select input file in case of nargin = 0
  %   else varargin{1} = path/file to grid to open
  %   supported formats: .mat / ASCII GIS / Massimiliano's binary format
  % 
  % OUTPUT: grid structure (internal MDW format) or error message (char)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  answer = [];
  if nargin == 0
    [file,path] = uigetfile('*.*','Load MetDataWizard grid file');
    if isequal(file,0) || isequal(path,0)
      return;
    else
      filepath = fullfile(path,file);
    end;
  else
    filepath = varargin{1};
  end;
  [path,name,ext] = fileparts(filepath);
  
  fid = fopen(filepath,'r');
  if fid == -1
    answer = 'File inaccessible';
    return;
  end;
  try %reading grid according to ASCII GIS format
    answer = 'Error reading grid.ncols';  
    grid.ncols = fgets(fid);
    fix = strfind(lower(grid.ncols),'ncols');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.ncols(1:fix-1) grid.ncols(fix+5:end)];
    grid.ncols = str2num(hlpstr);
    answer = 'Error reading grid.nrows';  
    grid.nrows = fgets(fid);
    fix = strfind(lower(grid.nrows),'nrows');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.nrows(1:fix-1) grid.nrows(fix+5:end)];
    grid.nrows = str2num(hlpstr);
    answer = 'Error reading grid.xllcorner';  
    grid.xllcorner = fgets(fid);
    fix = strfind(lower(grid.xllcorner),'xllcorner');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.xllcorner(1:fix-1) grid.xllcorner(fix+9:end)];
    grid.xllcorner = str2num(hlpstr);
    answer = 'Error reading grid.yllcorner';  
    grid.yllcorner = fgets(fid);
    fix = strfind(lower(grid.yllcorner),'yllcorner');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.yllcorner(1:fix-1) grid.yllcorner(fix+9:end)];
    grid.yllcorner = str2num(hlpstr);
    answer = 'Error reading grid.cellsize';  
    grid.cellsize = fgets(fid);
    fix = strfind(lower(grid.cellsize),'cellsize');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.cellsize(1:fix-1) grid.cellsize(fix+8:end)];
    grid.cellsize = str2num(hlpstr);
    answer = 'Error reading grid.NODATA_value';  
    grid.NODATA_value = fgets(fid);
    fix = strfind(lower(grid.NODATA_value),'nodata_value');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.NODATA_value(1:fix-1) grid.NODATA_value(fix+12:end)];
    grid.NODATA_value = str2num(hlpstr);
    answer = 'Error reading grid.data';  
    formatstr = '';
    for cix = 1:grid.ncols
      formatstr = [formatstr '%f'];
    end;
    data = textscan(fid,formatstr,grid.nrows);
    for cix = 1:grid.ncols
      grid.data(:,cix) = data{cix};
    end;
    grid.data = flipud(grid.data); %conversion necessary in order to comply with own asciigrid standards
    clear data;
    fclose(fid);
    answer = grid;
  catch
    try
      fclose(fid);
    end;
    answer = 'File is not a standard ASCII grid';  
  end;
end