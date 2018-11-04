function [Smodel,Jconverted] = dtv1_fitmri(theta,Q,opts)
% dtv1_fitmri
% theta1: C1 on domain [0,1] (fiso=C1)
% theta2: C2 on domain [0,1] (f1=C2-C1*C2, f2=1-C1-C2+C1*C2)
% theta3: C3 on domain [0,1] (lambda_par = C3*3e-3)
% theta4: C4 on domain [0,1] (lambda_perp1 = C4*lambda_par)
% theta5: C5 on domain [0,1] (lambda_perp2 = C5*lambda_par)
% theta6: theta1
% theta7: phi1
% theta8: theta2
% theta9: phi2
% theta10: S0

if nargout==1
    Smodel = dt_all_fitmri(theta,Q,opts); 
else
    [Smodel,Jconverted] = dt_all_fitmri(theta,Q,opts);   
end

end