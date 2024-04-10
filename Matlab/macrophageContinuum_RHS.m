function dy=macrophageContinuum_RHS(~,y,p,Afwd,Aback,f1,f2,R,params)

   g=y(1); c=y(2); m=y(3:end);
   
   g_integral=trapz(p,f1(p).*m);
   c_integral=trapz(p,f2(p).*m);
   M=trapz(p,m);
   
   fluxPos=(1-p).*m;
   fluxNeg=(1+p).*m;
      
   dy=zeros(size(y));
   
   dy(1)= g_integral - params.gammaG*g;
   dy(2)= params.kappaC*c_integral - c*g - c;
   dy(3:end)= -params.alpha1*c*Aback*fluxPos+params.alpha2*g*Afwd*fluxNeg+(c+params.cT)*R(p)*M*(1-M/params.mMax)-params.gammaM.*m;

end
