p = jsondecode(fileread('/private/tmp/claude-32353/-Users-sprinkle-work-data-mvt-nature/b36285a7-76cc-4f89-ab9c-6bf640615298/scratchpad/probe_values.json'));
v = p.values(:); r = round(v, 4);
out = struct('input', v, 'matlab_round', r, 'released', p.released(:));
fid = fopen('/private/tmp/claude-32353/-Users-sprinkle-work-data-mvt-nature/b36285a7-76cc-4f89-ab9c-6bf640615298/scratchpad/round_probe_out.json','w');
fwrite(fid, jsonencode(out), 'char'); fclose(fid);
agree = sum(r == p.released(:));
fprintf('MATLAB round(v,4) agrees with released file on %d of %d probe values\n', agree, numel(v));
