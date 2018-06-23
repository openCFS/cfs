% This is an example script for the creation of an interpolation filter for upsampling in order to avoid imaging.
% Similar filters could also be used for downsampling in order to avoid aliazing
% Some comments:
%    - In order to create a new filter with the desired properties one can play with the fdatool (firls is the Least-Squares Algorithm, Equiripple also seems to be promising
%    - The fvtool can be used to analyse the magnitude response of the filter
%    - The order of the filter is even, such it has an odd number of coefficients and an integer delay of time steps
%    - In the actual cfsdat implementation the filter coefficients need to be copied manually into the xml file

filename='int_4th58_forcfsdat';

factor = 4;

intcoeffs = firls(58,[0 0.2 0.3 1],[1 1 0 0],[1 1]);

%fvtool(intcoeffs,'Fs',100000);
[h,w]=freqz(intcoeffs);
plot(0:50000/(length(h)-1):50000,20*log10(abs(h)));
xlabel('Frequency (Hz)');
grid on;
ylabel('Magnitude (dB)');
ylim([-200 20]);

Fs = 100000;
Fn = Fs / 2;
Lf = length(h);
frequency = 0:Fn/(Lf-1):Fn;
magnitude = abs(h);
magnitudedb = 20*log10(magnitude);
magnitude4 = magnitude * 4;
magnitudedb4 = 20*log10(magnitude4);

%csvwrite('downupfilters_int_4th58.csv',[transpose(frequency) magnitude magnitudedb magnitude4 magnitudedb4]);


fileID = fopen(filename,'w');

%fprintf(fileID,'%s\n','Direct-Form FIR Polyphase Interpolator');

%fprintf(fileID,'%d\n',factor);

%fprintf(fileID,'%d\n',length(intcoeffs));

fprintf(fileID,'%.20f',iCoef);
for iCoef = intcoeffs(2:end)
    fprintf(fileID,',%.20f',iCoef);
end

fclose(fileID);




