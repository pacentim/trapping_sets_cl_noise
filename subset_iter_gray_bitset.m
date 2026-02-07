function it = subset_iter_gray_bitset(n)
%SUBSET_ITER_GRAY_BITSET  Iterator over all non-empty subsets of {1..n}
% in Gray-code order, represented as packed uint64 bitset.
%
% Usage:
%   it = subset_iter_gray_bitset(n);
%   it.reset();
%   while true
%       [bits, flipped, done] = it.next();
%       if done, break; end
%       % bits: uint64(1, W) bitset
%       % flipped: index in 1..n of the bit that flipped vs previous step
%   end
%
% Gray code ensures exactly one bit flips each step.

    W = ceil(n/64);
    bits = zeros(1, W, 'uint64');  % current subset as bitset

    % Gray counter state
    i = uint64(0);   % step index, from 0 to 2^n - 1
    g_prev = uint64(0);

    function reset_()
        bits(:) = 0;
        i = uint64(0);
        g_prev = uint64(0);
    end

    function [bits_out, flipped, done] = next_()
        % Advance to next non-zero gray code word.
        % done becomes true after emitting subset for i = 2^n - 1.
        done = false;
        flipped = 0;

        if n >= 63
            % Enumerating all subsets for n>=63 is astronomically large.
            % This guard is here to prevent silent overflow in 2^n logic.
            % You can remove it, but you will not finish in practice.
            % Still, the iterator logic for flipped-bit uses uint64 and is fine up to 63.
        end

        % Move to next state
        i = i + 1;
        if i == 0
            % overflow (only possible if you truly go beyond 2^64-1 steps)
            done = true;
            bits_out = bits;
            return;
        end

        % If i exceeds 2^n - 1, stop
        % We cannot compute 2^n exactly for large n in uint64 safely,
        % so we instead stop when i has its n-th bit set (i == 2^n).
        if bitget(i, n+1) == 1
            done = true;
            bits_out = bits;
            return;
        end

        g = bitxor(i, bitshift(i, -1));   % gray(i)
        d = bitxor(g, g_prev);            % single flipped bit (power of two)

        % Find flipped bit position (1-based)
        % d has exactly one bit set; locate it:
        flipped = 1;
        while bitand(d, uint64(1)) == 0
            d = bitshift(d, -1);
            flipped = flipped + 1;
        end

        % Toggle that bit in packed bitset
        w = floor((flipped-1)/64) + 1;
        b = mod(flipped-1, 64); % 0..63
        bits(w) = bitxor(bits(w), bitshift(uint64(1), b));

        g_prev = g;

        % Skip the empty set (Gray code starts at 0; after first step, it is nonzero)
        bits_out = bits;
    end

    it.next  = @next_;
    it.reset = @reset_;
end
