function [raw,hdr] = readlas(fname);
%  function [raw,hdr] = readlas(fname);
%  load .las data
    
hdr = readlas_hdr(fname);

[fid] = fopen(fname,'r','l');
dum = fread(fid,hdr.OffsetToPointData,'char');
[raw] = rekpoints(fid,hdr);

fclose(fid);


%-------------------------------------------------------------------------
function [raw] = rekpoints(fid,hdr); 
% function to reconstruct point data from point record  
  
  num = hdr.NumberOfPointRecords;
  x = ones(1,num);raw.y = x;raw.z = x;raw.int = x;
  raw.x = x;
  
  if 1
    raw.rnnr = x;raw.nrrt =x;
    raw.ReturnNumber = x; raw.NumberOfReturns = x;
    raw.ScanDirection = x;raw.EdgeFlightLine = x;
    raw.Classification = x; raw.ScanAngle = x;
    raw.UserData = x; raw.PointSourceID = x;
    raw.GPSTime = x;
  end
  
  clear x;
  
  skip = hdr.PointDataRecordLength;
  fseek(fid,hdr.OffsetToPointData,-1);
  
  raw.x = fread(fid,num,'long',skip-4);
  fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,4,0);
  
  raw.y = fread(fid,num,'long',skip-4);
  fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,8,0);
  
  raw.z = fread(fid,num,'long',skip-4);
  fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,12,0);
  
  raw.int = fread(fid,num,'ushort',skip-2);
  
  if 1
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,14,0);
    % Return Number and Number of Returns are coded with 3 bits each!
    raw.rnnr = fread(fid,num,'*ubit3',((skip-1)*8)+5);

    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,14,0);
    %move three bits
    fread(fid,1,'*ubit3');
    raw.rnnr = fread(fid,num,'*ubit3',((skip-1)*8)+5)*10 + raw.rnnr;
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,15,0);
    raw.Classification = fread(fid,num,'*ubit4',((skip-1)*8)+4);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,16,0);
    raw.ScanAngle = fread(fid,num,'char=>uchar',skip-1);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,17,0);
    raw.UserData = fread(fid,num,'char',skip-1);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,18,0);
    raw.PointSourceID = fread(fid,num,'ushort',skip-2);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,20,0);
    raw.GPSTime= fread(fid,num,'double',skip-8);
  end
  
  raw.x = (raw.x * hdr.XScaleFactor) + hdr.XOffset;
  raw.y = (raw.y * hdr.YScaleFactor) + hdr.YOffset;
  raw.z = (raw.z * hdr.ZScaleFactor) + hdr.ZOffset;
  
    
  return
  % old code
  %h = waitbar(0,'Loading single echos from .las file...');
  for i = 1:num
    %if hdr.PointDataFormatID == 1
      
    x(i) = fread(fid,1,'long');
    y(i) = fread(fid,1,'long');
    z(i) = fread(fid,1,'long');
    Int(i) = fread(fid,1,'ushort');
    
    if nargout == 4
      ReturnNumber = fread(fid,1,'bit3');
      NumberOfReturns = fread(fid,1,'bit3');
      ScanDirection = fread(fid,1,'bit1');
      EdgeFlightLine = fread(fid,1,'bit1');
      Classification = fread(fid,1,'uchar');
      ScanAngle = fread(fid,1,'char');
      UserData = fread(fid,1,'uchar');
      PointSourceID = fread(fid,1,'ushort');
    else
      ReturnNumber(i) = fread(fid,1,'bit3');
      NumberOfReturns(i) = fread(fid,1,'bit3');
      ScanDirection(i) = fread(fid,1,'bit1');
      EdgeFlightLine(i) = fread(fid,1,'bit1');
      Classification(i) = fread(fid,1,'uchar');
      ScanAngle(i) = fread(fid,1,'char');
      UserData(i) = fread(fid,1,'uchar');
      PointSourceID(i) = fread(fid,1,'ushort');
    end
    if hdr.PointDataFormatID == 1
      GPSTime(i) = fread(fid,1,'double');
    end
    if fix(i/10000)==i/10000
     waitbar(i/num,h)
    end
  end
  close(h);
  