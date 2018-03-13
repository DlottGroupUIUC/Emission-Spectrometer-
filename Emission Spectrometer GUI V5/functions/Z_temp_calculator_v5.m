function [temperature,temperature_error,slope,intercept,Z_full] = Z_temp_calculator_v4(spectral_radiance,wavelength,idx_of_wavelengths_to_exclude)
%Z_temp_calculator Calculate the temperature of greybody data

%Physical Constants
h=6.6260695729E-34;
c=299792458;
k=1.3806488E-23;
%Derived Constants
c1=2*h*c^2;
c2=h*c/k;

%Convert wavelength to meters
wavelength=wavelength*1e-9;

%Exclude wavelengths
orig_wavelength=wavelength;
orig_spec_rad=spectral_radiance;
spectral_radiance(idx_of_wavelengths_to_exclude)=[];
wavelength(idx_of_wavelengths_to_exclude)=[];

%Calculate Z
Z=log(c1./(wavelength.^5).*spectral_radiance'.^-1);
Z_full=log(c1./(orig_wavelength.^5).*orig_spec_rad'.^-1);
%Calculated inverse wavelength
inverse_wavelength=1./(wavelength);

%Convert imaginary data into Nan (remove negative spectral radiance)
Z(imag(Z)~=0)=NaN;
Z_full(imag(Z_full)~=0)=NaN;

%Isolate real portion after removing imaginary values
Z=real(Z);
Z_full=real(Z_full);

%check for excessive NaN values
if sum(isnan(Z))>20
    temperature=NaN;
    temperature_error=NaN;
    slope=NaN;
    intercept=NaN;
else

%Indexes for values that arn't NaN
k=~isnan(Z);

%Only fit the data that isn't NaN
lin_fit=fit(inverse_wavelength(k)',Z(k)','poly1');
%Get coeffecients
fit_coeffs=coeffvalues(lin_fit);
%Isolate slope and intercept
slope=fit_coeffs(1);
intercept=fit_coeffs(2);

%Get 95% confidence interval
conf=confint(lin_fit);
slope_er=conf(:,1);
int_er=conf(:,2);

%Calculate temperature and error
temperature=c2/slope;
Temp_conf=c2./slope_er;
temperature_error=mean([abs(temperature-Temp_conf(1)) abs(temperature-Temp_conf(2))]);

end
end

