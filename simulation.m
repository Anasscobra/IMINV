% Simulation de l'image observee
figure();
imagesc(Z);
colormap gray
title('Visualisation de l image idéale');
% Simulation de l'image simulee apres correction par le gain

%% On trace l'erreur gth-1
[m,n] = size(Z);
vecteur_erreur= gth-ones(n,1);
plot(vecteur_erreur);
title('Allure du vecteur erreur');
xlabel('Indice du vecteur'); 
ylabel('Amplitude de l erreur'); 
err=norm(vecteur_erreur,2);
disp('Norme de l erreur gth-1 : ')
disp(err);
%% On prend g=gth
W_simul=Z*diag(gth);
figure();
imagesc(W_simul);
title('Visualisation de l image en utilisant le gain')
colormap gray
% Visualisation du gth et du log(gth)
figure(); 
plot(gth); 
title('Allure du vecteur gain gth'); 
xlabel('Indice du vecteur'); 
ylabel('Amplitude du vecteur gth'); 

figure(); 
plot(log(gth)); 
title('Allure du vecteur log(gth)'); 
xlabel('Indice du vecteur'); 
ylabel('Amplitude du vecteur log(gth)'); 

%%% methode empirique
vecteur_moyen=zeros(1,n);
for i=1:n
    vecteur_moyen(i)=mean(log(W_simul(:,i)));
end
vecteur_moyen_centre=vecteur_moyen-mean(vecteur_moyen);
figure(); 
plot(vecteur_moyen_centre);
title('Visualisation du vecteur moyen centré Méthode empirique');
xlabel('Indice du vecteur'); 
ylabel('Amplitude du vecteur moyen centré'); 
%%% signal erreur entre gth et methode empirique 
figure(); 
hold on 
plot(vecteur_moyen_centre);
hold on 
plot(log(gth)-mean(log(gth))); 
legend('Vecteur moyen centré','log(gth) centré');

%%% erreur d'estimation
figure();
vecteur_estim_erreur=vecteur_moyen_centre-log(transpose(gth));
plot(vecteur_estim_erreur);
title('Evaluation de l erreur entre moyenne empirique et log(gth)');
xlabel('Indice du vecteur'); 
ylabel('Amplitude de l erreur'); 


N=nextpow2(n);
fourier_erreur_estim=abs(fft(vecteur_estim_erreur,2^N));
figure();
plot(fourier_erreur_estim);
title('Transformée de Fourier de l erreur entre moyenne empirique et log(gth)')
xlabel('Fréquence'); 
ylabel('Amplitude de la transformée de Fourier'); 
disp('La norme de L erreur entre moyenne empirique et log(gth) : '); 
disp(norm(vecteur_estim_erreur)); 
%%% methode maximum de vraisemblance
mat_log=log(W_simul);
vect_delta=zeros(1,n);
for i=2:n
    vect=mat_log(:,i)-mat_log(:,i-1);
    delta_f_i=median(vect);
    vect_delta(i)=delta_f_i;
end
f_max_vraisemblance_norme_1=cumsum(vect_delta);
f_max_vraisemblance_norme_1_centre=f_max_vraisemblance_norme_1-mean(f_max_vraisemblance_norme_1);
figure(); 
hold on 
plot( f_max_vraisemblance_norme_1_centre);
hold on 
plot(log(gth)); 
legend('Maximum de vraissemblance','log(gth)'); 

figure(); 
plot(f_max_vraisemblance_norme_1_centre-log(transpose(gth)));
title('erreur max vraisemblance norme 1');
xlabel('Indice du vecteur'); 
ylabel('Amplitude de l erreur'); 
disp('La norme de l erreur entre le max de vraissemblance et log(gth)');
disp(norm(f_max_vraisemblance_norme_1_centre-log(transpose(gth)))); 


fourier_erreur_estim_max=abs(fft(f_max_vraisemblance_norme_1_centre-log(transpose(gth)),2^N));
figure();
plot(fourier_erreur_estim_max);
title('Transformée de Fourier de l erreur entre Max de vraissemblance et log(gth)')
xlabel('Fréquence'); 
ylabel('Amplitude de la transformée de Fourier'); 

%%% Méthode par maximum à posteriori 
e = ones(n-1,1);
D= spdiags([e -e],0:1,n-1,n);
full(D); 
M=transpose(D)*D;
s=vecteur_moyen;
mu=logspace(-1,1,100); 
erreur=zeros(1,100); 
fth=log(gth);
for i=1:100
   f_hat=mldivide(M+mu(i)*eye(n),M*transpose(s));
   moyenne_estimee=mean(f_hat);
   f_hat=f_hat-moyenne_estimee;
   erreur(i)=norm(f_hat-fth-mean(fth),2)^2/n; 
end 
figure(); 
loglog(mu,erreur); 
title ('Evolution de l erreur en fonction de mu par maximum à posteriori'); 
xlabel('Valeur de mu'); 
ylabel('Amplitude de L erreur'); 

disp('La valeur minimale de L erreur par maximum à posteriri est :'); 
disp(min(erreur)); 
disp('la valeur de mu correspondante à la valeur d erreur minimale est : '); 
disp(mu(minima(erreur)));

% 2ème Régularisation en norme 2 en utilisant la fonction MAPL1 
lbd=logspace(-5,5,100); 
D_mapl= spdiags([e -e],0:1,n-1,n);
erreur_mapl=zeros(1,100);
for i=1:100
   f_hat1=MAPL1(transpose(mat_log),D_mapl,lbd(i));
   moyenne_estimee=mean(f_hat1);
   f_hat1=f_hat1-moyenne_estimee;
   erreur_mapl(i)=norm(f_hat1-fth-mean(fth),2)^2/n; 
end 
figure(); 
loglog(mu,erreur_mapl); 
title ('Evolution de l erreur en fonction de mu par la fonction MAPL1'); 
xlabel('Valeur de mu'); 
ylabel('Amplitude de L erreur'); 

disp('La valeur minimale de L erreur par la fonction MAPL1 est :'); 
disp(min(erreur_mapl)); 
disp('la valeur de mu correspondante à la valeur d erreur minimale est : '); 
disp(mu(minima(erreur_mapl)));


