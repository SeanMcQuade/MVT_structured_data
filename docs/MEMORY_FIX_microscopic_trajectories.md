# Memory fix: `plot_microscopic_trajectories.m` batched patches

**File:** `Scripts/plot_microscopic_trajectories.m`
**Commit:** `388d890` (on branch `tooling`); pre-change checkpoint at `ee1c2cc`, original code at `0742185`.
**Symptom fixed:** a full-day westbound run exhausted RAM and hard-rebooted the machine.

---

## What went wrong

Running `plot_microscopic_trajectories(18)` rendered the two low-res PNGs, then
the machine rebooted mid-zoom-stage (confirmed: `..._lowres.png` and
`..._zoomwin_lowres.png` were written, but `..._zoom_lowres.png` never was).

Two things combined to blow up memory, **neither of which `skip_t_plot` can
help with** — that knob only trims the number of *vertices per polygon*, not the
number of *objects*:

1. **One `patch()` object per trajectory.** The inner loop called `patch(...)`
   once per trajectory. The released `slim` files have **no `direction` field**
   (slim is already westbound-only), so the `isfield(data,'direction')` guard
   fell through to the `else` branch and selected *every* trajectory —
   ~15,438 per file. The loop reads every other file (`1:2:end`, 12 files), so a
   full day accumulated **~185,000 persistent HG2 patch objects in one axes**.
   Each HG2 object carries kilobytes of fixed overhead (handle, property structs,
   listeners, transform) regardless of how few vertices it holds. That fixed
   overhead × 185k objects is the RAM.

2. **The zoom stage walked all of them again.** To "lighten" the zoom render it
   ran `findobj(gca,'type','patch')` (returning all ~185k handles) and then a
   MATLAB `for` loop reading `h_traj(j).Vertices` on each — a second
   O(#trajectories) pass that materialized every object's data. This is the step
   that was executing when the machine rebooted.

A per-file `drawnow limitrate` inside the build loop made it worse: it forced a
full-scene re-render on every iteration while the object count kept growing, so
cost grew super-linearly (the full-day render took ~34 minutes, i.e. swap
thrashing).

## The fix

Draw **all of a file's trajectories as a single `patch` object** using the
`Faces`/`Vertices`/`FaceVertexCData` form, instead of one `patch()` call per
trajectory. Each trajectory's polygon vertices are collected into shared arrays;
a NaN-padded face-index matrix (`F`) tells the one patch which vertices belong to
which trajectory, so they render as separate filled faces. Object count per full
day drops from **~185,000 to ~12** (one per loaded file).

Consequently:
- The zoom-stage `findobj`/`.Vertices` deletion loop is **removed** — with few
  objects it buys nothing, and axes clipping already discards off-window faces
  once `xlim`/`ylim` are set (that `xlim(zoom_t), ylim(zoom_x/1000)` line is
  unchanged and still performs the zoom).
- The per-file `drawnow limitrate` is **removed**; the single `drawnow` before
  saving remains.

### Output is unchanged

Pixels should be identical: same `skip_t_plot` subsampling, same `int32`
rounding when `flag_use_int32_vars` is set, same `colormap`/`caxis`, and the
per-vertex speed coloring is preserved via `'FaceColor','interp'` with
per-vertex `FaceVertexCData` (the old high-level `patch(X,Y,C)` with a per-vertex
`C` also uses interpolated, scaled coloring). The only subtlety: the low-level
`Vertices`/`FaceVertexCData` properties reject integer types, so the stacked
arrays are cast to `double` **at the patch call** — after the `int32` rounding
has already been applied, so the visible quantization is retained.

`flag_use_int32_vars` is kept as-is. It no longer buys graphics memory (there is
now one object per file), but it still changes the rendered values via rounding,
so it is left in place pending a separate decision.

---

## Exact diff

### 1. Inner loop — batch into one patch per file; drop per-file `drawnow`

```diff
-    % Process trajectories
-    fprintf('Adding %d trajectories to plot ...',length(ind)), tic
-    % iterate over all trajectories
-    for j = 1:length(data) % loop over used trajectories
+    % Process trajectories: collect each trajectory's polygon into shared
+    % vertex/face/color arrays, then draw them ALL as one patch object per
+    % file. Calling patch() once per trajectory creates ~10^5 persistent HG2
+    % objects for a full day (each with kByte-scale fixed overhead regardless
+    % of vertex count) and exhausts RAM; one patch/file keeps memory bounded.
+    fprintf('Adding %d trajectories to plot ...',length(ind)), tic
+    nTraj = length(data);
+    vertsC = cell(nTraj,1);    % per-trajectory polygon vertices [2*ns x 2]
+    cdataC = cell(nTraj,1);    % per-trajectory vertex colors     [2*ns x 1]
+    faceLens = zeros(nTraj,1); % vertex count of each polygon
+    for j = 1:nTraj % loop over used trajectories
         traj_len = data(j).length*ft2meterFactor; % length of vehicle [m]
-        traj_t = data(j).timestamp; 
+        traj_t = data(j).timestamp;
         traj_x = data(j).(subfield_name_x); % vehicle position [m]
         % Calculate velocity
         traj_v = (traj_x([2:end,end])-traj_x([1,1:end-1]))./...
             (traj_t([2:end,end])-traj_t([1,1:end-1]))*direction;
         % Subsample trajectories for plotting
         traj_n = length(traj_t); % number of trajectory points
         traj_ns = max(ceil(traj_n/skip_t_plot),2); % number of segments
         traj_ind = round(linspace(1,traj_n,traj_ns)); % subsample vector
         traj_t = traj_t(traj_ind);
         traj_x = traj_x(traj_ind);
         traj_v = traj_v(traj_ind);
-        patch_t = [traj_t;traj_t(end:-1:1)] - ax_t(1) ; 
+        patch_t = [traj_t;traj_t(end:-1:1)] - ax_t(1) ;
         patch_x = [traj_x;traj_x(end:-1:1)+traj_len*direction]/1000;
         patch_v = [traj_v;traj_v(end:-1:1)];
         if flag_use_int32_vars
             patch_t = int32(patch_t*toInt32Factor);
             patch_x = int32(patch_x*toInt32Factor);
             patch_v = int32(patch_v*toInt32Factor);
         end
-         patch(patch_t, patch_x, patch_v,'EdgeColor','None')
-    end
-    drawnow limitrate nocallbacks
+        vertsC{j} = [patch_t(:), patch_x(:)];
+        cdataC{j} = patch_v(:);
+        faceLens(j) = numel(patch_t);
+    end
+    % Stack all polygons and build a NaN-padded face-index matrix: row j
+    % lists the vertex indices of trajectory j, so one patch renders them all
+    % as separate filled faces. Cast to double at the boundary because the
+    % Vertices/FaceVertexCData properties reject integer types (the int32
+    % rounding above is preserved, so pixel output is unchanged).
+    V = double(vertcat(vertsC{:}));
+    C = double(vertcat(cdataC{:}));
+    maxLen = max(faceLens);
+    F = nan(nTraj,maxLen);
+    voff = 0;
+    for j = 1:nTraj
+        F(j,1:faceLens(j)) = voff + (1:faceLens(j));
+        voff = voff + faceLens(j);
+    end
+    patch('Faces',F,'Vertices',V,'FaceVertexCData',C,...
+        'FaceColor','interp','EdgeColor','none')
     fprintf(' Done (%0.0fsec).\n',toc)
 end
```

### 2. Zoom stage — remove the `findobj`/`.Vertices` deletion loop

```diff
-    % Remove all trajectories fully outside of zoom window
-    fprintf('Remove trajectories outside zoom window, and zoom in ...'), tic
-    h_traj = findobj(gca,'type','patch');
-    ind_delete = false(1,length(h_traj));
-    for j = 1:length(h_traj)
-        ind_delete(j) = all(h_traj(j).Vertices(:,1)<zoom_t(1))||...
-            all(h_traj(j).Vertices(:,1)>zoom_t(2))||...
-            all(h_traj(j).Vertices(:,2)<zoom_x(1)/1000)||...
-            all(h_traj(j).Vertices(:,2)>zoom_x(2)/1000);
-    end
-    delete(h_traj(ind_delete))
+    % Zoom in. (Previously this walked every patch object and deleted those
+    % outside the window to lighten the render -- a second O(#trajectories)
+    % pass that spiked RAM and rebooted the machine. With one patch/file the
+    % objects are few, and setting xlim/ylim below lets axes clipping discard
+    % off-window faces at draw time, so no deletion is needed.)
+    fprintf('Zoom in ...'), tic
```

The `xlim(zoom_t), ylim(zoom_x/1000)` line immediately below (unchanged) is what
performs the zoom.

---

## How to re-run

The `_reduced.mat` caches for a day are reused if all 24 are present, so a re-run
skips the multi-GB JSON decode and starts drawing immediately:

```bash
cd .../MVT_structured_data/Scripts
matlab -batch "plot_microscopic_trajectories(18)"
```

The edited code is newer than the existing day-18 figures, so the staleness check
marks them stale and re-runs without needing `Force`. To rebuild regardless, pass
`Force`:

```bash
matlab -batch "plot_microscopic_trajectories(18, 'Force', true)"
```

## Rollback

```bash
git checkout 0742185 -- Scripts/plot_microscopic_trajectories.m   # restore original
# or reset the branch to the pre-change checkpoint:
git reset --hard ee1c2cc
```
