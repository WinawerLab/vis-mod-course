function oi = getOptics(whichObserver, verbose)

% Add optics

switch whichObserver
    case 'ideal'
        
        % (Question: how to add pupil size?)
        wvf = wvfCreate('calc wave', [400:10:700]);
        wvf = wvfSet(wvf,'measured pupil size', 1.5);
        wvf = wvfSet(wvf,'calc pupil size', 1.5);
        wvf = wvfComputePSF(wvf);
        
        oi = wvf2oi(wvf);
        oi = opticsSet(oi, 'model', 'diffraction limited');
        
    case 'human'
        oi = oiCreate('human');
        
end

if verbose
    
    % Plot ocular transmittance function
    figure; plot(oi.optics.lens.wave, oi.optics.lens.transmittance);
    title('Ocular transmittance function');
    
    % Request pupil diameter
    p = opticsGet(oi.optics,'pupil diameter','mm');
    
    if verbose
        % Plot OTF
        oiPlot(oi, 'OTF', 'wavelength', 550);
        
        % Plot Point spread function
        oiPlot(oi, 'PSF', 'wavelength', 550);
        
        % Plot line spread function for multiple wave lengths
        oiPlot(oi,'ls wavelength');
        title(sprintf('F/# = %.0d',opticsGet(oi.optics,'f number')));
    end
    
end
