NframesPerLoop = 40;
Nloops = 6;
Nframes = NframesPerLoop*Nloops;

rf_phase_0 = 117;
rf_count = 1;  

phases = zeros(Nframes,1);

for frame = 1:Nframes
    rf_phase = mod(0.5 * rf_phase_0 * rf_count^2, 360.0);
    phases(frame) = rf_phase;

    rf_count = rf_count + 1;
    if rf_count > Nframes
        rf_count = 1;
    end
end

