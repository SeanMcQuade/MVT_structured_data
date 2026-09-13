function version = dataVersion()
% MVT.DATAVERSION  Version of the released MVT data set this code produces.
%
% Purpose
%   One place to state which version of the derived data the pipeline emits, so
%   the number cannot drift between scripts, logs, and documentation. This is
%   the version of *our* derived products (slim/full/gps), not of the upstream
%   I-24 MOTION data, which is versioned independently by the observatory and
%   is unchanged by anything here.
%
% Outputs
%   version  char, e.g. '2.1.1'
%
% Versioning scheme
%   MAJOR.MINOR.PATCH, where a consumer reads it as:
%     MAJOR  the data means something different; re-read the documentation.
%     MINOR  fields added, removed, renamed, or re-specified; code that parses
%            the data may need changing.
%     PATCH  same fields, same meaning, same format - only the last decimal of
%            a small number of values can differ. Analyses do not need redoing;
%            checksums do change, so two copies with different PATCH versions
%            are not expected to be byte-identical.
%
% History
%   2.1    The data set used for the Nature submission.
%   2.1.1  Deterministic fuel quadrature: mvt.neumaierDot (compensated
%          summation) replaces MATLAB's `dot` in the trapezoidal fuel integral.
%          `dot` calls the BLAS, whose accumulation order differs between
%          platforms, so total_fuel_consumed_grams depended on whether the run
%          happened on a Mac or a PC. No field, format, or meaning changes; the
%          affected values move by 0.0001 g on roughly one trajectory in ten
%          thousand. Rationale and measurements:
%          docs/REPRODUCIBLE_QUADRATURE.md.
%
% Notes
%   Bump PATCH whenever a change alters, or could alter, the bytes of the
%   released data even though nothing about its structure or interpretation
%   changes. The point of the version is to let someone holding two copies of
%   the data tell whether they should be identical.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

version = '2.1.1';
end
