function total = neumaierDot(a, b)
% MVT.NEUMAIERDOT  Platform-independent dot product for the fuel quadrature.
%
% Purpose
%   Replaces `dot(a,b)` in the trapezoidal fuel integral. `dot` dispatches to
%   the BLAS, whose accumulation order is implementation-defined: Accelerate on
%   Apple silicon and MKL on Windows disagree in the last bit, so the same code
%   on the same data produces different `total_fuel_consumed_grams` on a Mac and
%   on a PC. Because those totals are then rounded to 4 decimals, a value sitting
%   within a ULP of a rounding tie lands on different sides, and the released
%   JSON differs byte-for-byte between machines.
%
%   This computes the same sum with compensated (Kahan-Babuska-Neumaier)
%   summation: a fixed sequence of IEEE-754 double operations with no
%   library-, thread-, or platform-dependent freedom. Every machine and every
%   language implementing this loop gets identical bits.
%
% Inputs
%   a, b  numeric vectors of equal length
%
% Outputs
%   total  double, sum of a.*b
%
% Notes
%   It is also *more accurate* than the BLAS form, not merely more consistent:
%   on 60 real trajectories it matched the exactly-rounded sum on all 60, while
%   `dot` deviated from it on 19 by 1-2 ULP. See
%   docs/REPRODUCIBLE_QUADRATURE.md.
%
%   The running compensation `c` is added once at the end rather than folded in
%   each step; that is what makes the algorithm robust when a term is much
%   larger than the running total, which is the case here whenever a trajectory
%   starts with near-zero fuel rates.
%
% (C) 2026 CIRCLES Consortium. BSD-3-Clause.

if numel(a) ~= numel(b)
    error('mvt:neumaierDot:sizeMismatch', ...
        'Inputs must have the same number of elements (%d vs %d).', ...
        numel(a), numel(b));
end

total = 0;          % running sum
c = 0;              % running compensation for what the sum dropped
for i = 1:numel(a)
    x = a(i) * b(i);
    t = total + x;
    if abs(total) >= abs(x)
        c = c + ((total - t) + x);
    else
        c = c + ((x - t) + total);
    end
    total = t;
end
total = total + c;
end
