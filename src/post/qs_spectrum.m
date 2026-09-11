function S = qs_spectrum(sig, Fs, band_frac)
%QS_SPECTRUM  Hann-windowed periodogram of one signal segment.
%
%   S = QS_SPECTRUM(sig, Fs, band_frac)
%
%   Reproduces exactly the spectrum the original drivers computed inline:
%     sg = sig - mean(sig);  Hann window;  P = |fft|^2;
%     f  = (0:Nc-1)*Fs/Nc;   search bins  hf = 2:floor(Nc/2)
%   (bin 1, the DC bin, is excluded, which is why the mean is removed first).
%
%   OUTPUT
%     S.f, S.P    full frequency axis [Hz] and power
%     S.hf        the searched bin range
%     S.fpk       frequency of the largest bin in hf   [Hz]
%     S.ppk       its power
%     S.pfrac     fraction of the total power in hf that sits within
%                 +/- band_frac*fpk of the peak. This is the "is it one tone
%                 or is it broadband" number: a clean LCO gives ~0.6-0.95,
%                 chaotic/broadband response gives a small value.
%     S.df        frequency resolution [Hz]

    if nargin < 3 || isempty(band_frac), band_frac = 0.05; end

    sg = sig(:) - mean(sig(:));
    Nc = numel(sg);
    S = struct('f',[], 'P',[], 'hf',[], 'fpk',NaN, 'ppk',NaN, 'pfrac',NaN, 'df',NaN);
    if Nc < 8
        return
    end

    win = 0.5*(1 - cos(2*pi*(0:Nc-1).'/(Nc-1)));
    P   = abs(fft(sg.*win)).^2;
    f   = (0:Nc-1).' * Fs/Nc;
    hf  = 2:floor(Nc/2);

    [ppk, kk] = max(P(hf));
    fpk = f(hf(kk));

    tot = sum(P(hf));
    if tot > 0 && fpk > 0
        bw   = max(band_frac*fpk, 2*(Fs/Nc));     % at least +/- 2 bins wide
        inb  = abs(f(hf) - fpk) <= bw;
        pfrac = sum(P(hf(inb))) / tot;
    else
        pfrac = NaN;
    end

    S.f = f;  S.P = P;  S.hf = hf;
    S.fpk = fpk;  S.ppk = ppk;  S.pfrac = pfrac;  S.df = Fs/Nc;
end
