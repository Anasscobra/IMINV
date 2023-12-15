function indice=minima(liste)
indice=1; 
valeur_minimale=liste(1); 
for i=2:length(liste)
    if liste(i)<valeur_minimale
        valeur_minimale=liste(i); 
        indice=i; 
    end 
end
end 

