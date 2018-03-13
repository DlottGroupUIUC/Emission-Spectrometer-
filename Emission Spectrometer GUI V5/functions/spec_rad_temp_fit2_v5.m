function [ model,Temperature,Temperature_unc,emissivity,emissivity_unc ] = spec_rad_temp_fit2_v4( wavelength,spectral_radiance,starting_temperature, idx_of_wavelengths_to_exclude)
%spec_rad_temp_fit Function to solve for the temperature given a
%temperature

%wavelenght in nm
%spectral radiance in W/m3/sr

wavelength=wavelength*1e-9; %convert to meters
spectral_radiance=spectral_radiance'; %Convert to column array

%Exclude wavelengths
orig_wavelength=wavelength;
spectral_radiance(idx_of_wavelengths_to_exclude)=[];
wavelength(idx_of_wavelengths_to_exclude)=[];

%Define fitting options for non-linear fit

s = fitoptions('Method','NonlinearLeastSquares',... %fit method
               'Lower',[0,1500],...   %lower bound for emissivity and temp
               'Upper',[1,1e4],...    %upper bound for emissivity and temp
               'Startpoint',[0.1 starting_temperature]);   %Starting point
           
 bb=fittype('E*2*6.63e-34*3e8^2/wavelength^5/(exp((6.63e-34*3e8)/(wavelength*1.38e-23*T))-1)',...
    'independent',{'wavelength'},'coefficients',{'E','T'},'options',s);  %blackbody fitting model

fitteddata=fit(wavelength,spectral_radiance,bb);  %fit results

%Get coeffecients
fit_coeffs=coeffvalues(fitteddata);
%Isolate slope and intercept
emissivity=fit_coeffs(1);
Temperature=fit_coeffs(2);

%Get 95% confidence interval
conf=confint(fitteddata);
emissivity_unc=conf(2,1)-emissivity;
Temperature_unc=conf(2,2)-Temperature;

%Generate model
model=emissivity*2*6.63e-34*3e8^2./orig_wavelength'.^5./(exp((6.63e-34*3e8)./(orig_wavelength'.*1.38e-23*Temperature))-1);

