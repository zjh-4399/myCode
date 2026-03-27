function rayData = read_ray_file_bellhop2d(rayFile)
%READ_RAY_FILE_BELLHOP2D Read a Bellhop-style 2D ASCII .ray file.

fid = fopen(rayFile, 'rt');
if fid < 0
    error('Cannot open ray file: %s', rayFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

rayData = struct();
rayData.file = rayFile;
rayData.title = local_strip_quotes(strtrim(fgetl(fid)));
rayData.freq_hz = sscanf(fgetl(fid), '%f', 1);
rayData.source_counts = sscanf(fgetl(fid), '%d').';
rayData.angle_counts = sscanf(fgetl(fid), '%d').';
rayData.top_depth = sscanf(fgetl(fid), '%f', 1);
rayData.bottom_depth = sscanf(fgetl(fid), '%f', 1);
rayData.coord = local_strip_quotes(strtrim(fgetl(fid)));

if numel(rayData.angle_counts) >= 1 && isfinite(rayData.angle_counts(1))
    nAlpha = max(rayData.angle_counts(1), 1);
else
    nAlpha = 1;
end

beams = cell(nAlpha, 1);
ib = 0;
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    line = strtrim(line);
    if isempty(line)
        continue;
    end

    alpha = sscanf(line, '%f', 1);
    metaLine = fgetl(fid);
    if ~ischar(metaLine)
        error('Unexpected EOF while reading beam metadata in %s.', rayFile);
    end
    meta = sscanf(metaLine, '%d').';
    if numel(meta) < 3
        error('Invalid beam metadata line in %s: %s', rayFile, metaLine);
    end

    nPts = meta(1);
    rz = nan(nPts, 2);
    for ip = 1:nPts
        dataLine = fgetl(fid);
        if ~ischar(dataLine)
            error('Unexpected EOF while reading beam points in %s.', rayFile);
        end
        vals = sscanf(dataLine, '%f');
        if numel(vals) < 2
            error('Invalid beam point line in %s: %s', rayFile, dataLine);
        end
        rz(ip,1:2) = vals(1:2).';
    end

    ib = ib + 1;
    beam = struct();
    beam.alpha_deg = alpha;
    beam.nPoints = nPts;
    beam.surface_bounces = meta(2);
    beam.bottom_bounces = meta(3);
    beam.r = rz(:,1);
    beam.z = rz(:,2);
    beams{ib} = beam;
end

rayData.beams = beams(1:ib);
end


function s = local_strip_quotes(s)
if isempty(s)
    return;
end
if s(1) == ''''
    s = s(2:end);
end
if ~isempty(s) && s(end) == ''''
    s = s(1:end-1);
end
end
