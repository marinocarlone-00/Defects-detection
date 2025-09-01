function bolle(isColor)

tic;
clc; 
close all;
if nargin<1

    isColor =true;

end

%parametri per valutare/stimare la grandezza delle bolle individuate'
immagine_percorso= '/Users/marinocarlone/Desktop/Codici Ad.tech/Rilevazione Bolle/primafoto.png';

min_DIAM_MM=2.0;% range fisico atteso (diametri bolla)
max_DIAM_MM=12.0; 
Soglia_LOG= 98.9;%soglia su LoG per isolare picchi 
ECC_MAX= 0.90;%eccentricità max accettata per bolla
CirC_MIN= 0.30;%circolarità minima accettata
I_SIG= 0.15;%intensità
S_SIG= 5;%spazio

%% CALIBRAZIONE FILE

%qui di carica il file precedentemente ottenuto dall'add on della
%calibrazione

LOAD=load('camera_calib.mat');
cam=[];
for fn=fieldnames(LOAD)'
    v=LOAD.(fn{1});
    if isa(v,'cameraParameters')
        cam= v;
        break;
    end
end

assert(~isempty(cam),'Nel .mat non c''è nessun cameraParameters');

unita = "mm";
if isprop(cam,'WorldUnits')&& ~isempty(cam.WorldUnits)
    unita=string(cam.WorldUnits); 
end

switch lower(unita)

    case {"mm","millimetro","millimetri"}

        scaleToMM=1;

    case {"cm","centimetro","centimetri"} 

        scaleToMM=10;

    case {"m","metro","metri"} 

        scaleToMM=1000;

    otherwise

        scaleToMM=1;

end

WP=cam.WorldPoints;  

if size(WP,2)==2

    WP3=[WP,zeros(size(WP,1),1)];

else

    WP3=WP; 

end

Dimen=squareform(pdist(WP));

L=min(Dimen(Dimen>0)); 
Lmm=L*scaleToMM;

[ii,jj] =find(Dimen>0 & abs(Dimen-L)<=2e-3*L &triu(true(size(Dimen)),1));

pixDist=[];

for k=1:size(cam.RotationVectors,1)

    rvec=cam.RotationVectors(k,:); 

    tvec =cam.TranslationVectors(k,:);

    try      
        uv =worldToImage(cam,rvec, tvec,WP3);

    catch   
        uv= worldToImage(cam,rotationVectorToMatrix(rvec),tvec,WP3);

    end

    if any(~isfinite(uv(:)))

        continue; 

    end

    pixDist=[pixDist; hypot(uv(ii,1)-uv(jj,1),uv(ii,2)-uv(jj,2))]; 

end

assert(~isempty(pixDist),'Nessuna distanza pixel valida');

pxPerMM =median(pixDist)/Lmm;
mmPerPx = 1/pxPerMM;

fprintf('Scala: %.4f px/mm (%.4f mm/px)\n', pxPerMM,mmPerPx);

%% PROCESSAMENTO DELLE IMMAGINI

immagine=imread(immagine_percorso);

if size(immagine,3)==3

    Ig=rgb2gray(immagine);

else

    Ig=immagine;

end

%Riadattamento a fattore di scala se necessario 
if all(size(Ig,1:2)==cam.ImageSize(1:2))

    IgU= undistortImage(Ig, cam);

else

    IgU =Ig;

end

Grigia =im2double(IgU);


%Flat-field+CLAHE 
bg=max(85,round(85*pxPerMM));

%questo filtro 'appiattisce' uniformando l'illuminazione, eliminando le discontinuità di luce ed eventuali
%vignettature

filtro_piatto=imflatfield(Grigia, bg);

%questo filtro aumenta il contrasto, anche in zone poco illuminate,
%dividendo l'immagine in una griglia 16X16 e applicando una equalizzazione
%locale in ogni tile
filtro_adattivo=adapthisteq(filtro_piatto, 'NumTiles',[10 10],'ClipLimit',0.01);

% questo filtro omogenizza l'intera texture del cartone, attenua la grana
% ma gonfia i bordi delle bolle, come soluzione vi è un processo di
% adattamento a scala (+ traslazione) finale
filtro_texture =imbilatfilt(filtro_adattivo, I_SIG, S_SIG);

%creazione di una maschera che integra la nuova texture
[H,W]=size(filtro_texture);
margine=max(10, round(0.02*min(H,W)));
maschera_roi=true(H,W); 
maschera_roi([1:margine,end-margine+1:end],:)=false;
maschera_roi(:,[1:margine,end-margine+1:end])=false;

%applico filtro LoG multiScala: (Laplacian of Gaussian), si applica una sfocatura 
%gaussiana di scala SIGMA e poi il laplaciano (seconda derivata).Le bolle non 
%hanno tutte lo stesso diametro. Variando SIGMA si ottiene una "mappa di blobness"
%alle diverse scale e per ogni pixel si prende il massimo su SIGMA
%(con normalizzazione in scala) così da essere quasi invariante alla
%dimensione. Questo filtro, aggrega il contrasto della bolla e sopprime molta 
%texture del cartone, creando una nuova maschera ideale che metta in
%evidenza solo le bolle

minDiam_px=max(6,round(min_DIAM_MM/mmPerPx));
maxDiam_px=round(max_DIAM_MM/mmPerPx);
minR_px=max(3,round(minDiam_px/2));
maxR_px=max(minR_px+1,round(maxDiam_px/2));

%Aumento della scala del filtro di riferimento
alpha_scala =0.55;
variabile=unique(round(linspace(minR_px, round(1.80*maxR_px), 10)));  
sigma=max(1.0,variabile/sqrt(2));
resp_luminosa=zeros(size(filtro_texture));
resp_scura=zeros(size(filtro_texture));

%È il cuore del LoG multiscala:per ogni SIGMA costruisce un filtro LoG e filtra l'immagine e accumula 
% – pixel per pixel – la risposta massima alle varie scale, separando "blob
% chiari" e "blob scuri''

for s=sigma

    var=2*ceil(3*s)+1;
    LoG=fspecial('log',[var var],s);
    Rf=(s^2)*imfilter(filtro_texture, LoG, 'replicate', 'conv');  % LoG norm. scala
    peso=s^alpha_scala;                                    
    resp_luminosa=max(resp_luminosa,peso*(-Rf));% cupole chiare
    resp_scura=max(resp_scura,peso*(Rf));% avvallamenti scuri

end

resp_luminosa=zeros(size(filtro_texture));%cupole chiare(LoG <0)
resp_scura=zeros(size(filtro_texture));%avvallamenti scuri(LoG> 0)


%Applica il filtro precedente all'immagine pre-filtrata:si ottengono due
%mappe di blobness multiscala (luminosa e scura), invarianti alla
%dimensione del difetto su cui poi si applicano filtri morfologici e/o di
%forma per segmentare le bolle

for s=sigma

    var=2*ceil(3*s)+1;
    LoG=fspecial('log',[var var],s);
    R=(s^2)*imfilter(filtro_texture,LoG,'replicate', 'conv');  % normalizzato per scala
    resp_luminosa=max(resp_luminosa,-R);
    resp_scura=max(resp_scura,R);

end


% Scelta della polarità dominante, e della mappa sulla quale effettuo la
% segmentazione
if max(resp_luminosa(:))>=max(resp_scura(:))
    score =resp_luminosa;
    pol='bright';
else
    score=resp_scura;  
    pol='dark';
end

%normalizzazione della mappa precedente, applicando un blur Gaussiano e rendendo
%l'istogramma più "pulito", stabilizzando la soglia (percentile/Otsu) utile per
%la morfologia successiva
score_nuovo= mat2gray(score);
score_nuovo_modificato= imgaussfilt(score_nuovo,1.0); 

%Soglia robusta sui picchi LoG
soglia_perc =prctile(score_nuovo_modificato(:), Soglia_LOG); %è una soglia robusta agli outlier, tiene solo i picchi più alti
soglia_otsu= graythresh(score_nuovo_modificato); %soglia OTSU massimizzando la separazione tra due classi fondo e oggetti
soglia_max=max(soglia_perc,soglia_otsu); %scegliere la soglia più severa
filtro_LoG=score_nuovo_modificato>= soglia_max;;

% pulizia+buchi+vincoli area (applicazione dei filtri morfologici - toolbox
% Matlab)
filtro_LoG =imopen(filtro_LoG, strel('disk',2));
filtro_LoG=imfill(filtro_LoG,'holes');
minAreaPx=round(pi*(minR_px)^2);
maxAreaPx=round(pi*(maxR_px)^2);


%Filtro del tool pre-processing: rimozione chiazze attraverso
% l'image processing tool di Matlab,
% morfologia di pulizia,
filtro_LoG= bwareaopen(filtro_LoG,minAreaPx);
filtro_LoG= filtro_LoG&~bwareaopen(filtro_LoG,maxAreaPx);

%Filtri di forma per ricercare le forme pseudo-circolari delle forme
%Pipeline morfologico dei filtri applicati
Perimetro=bwperim(filtro_LoG); %estrazione del perimetro, trasforma i blob pieni in anelli sottili
Perimetro =imdilate(Perimetro,strel('disk',1)); %Ispessisce gli anelli di 1–2 px e chiude micro-gap, quindi un bordo più continuo accumula più voto.
Perimetro=Perimetro | edge(score_nuovo,'Canny',[0.2 0.4], 1.0);  %estrazioni bordi

% Hough circolare sull'orlo LoG - estrazione delle forme circolari
[centL, radL, metL]=imfindcircles(Perimetro, [minR_px maxR_px],'ObjectPolarity','bright', 'Method','PhaseCode','Sensitivity',0.96,'EdgeThreshold',0.00);

% filtraggio della maschera
if ~isempty(centL)

    funzione=arrayfun(@(k) maschera_roi(round(centL(k,2)), round(centL(k,1))), 1:size(centL,1));

    centL=centL(funzione,:);  

    radL =radL(funzione);  

    metL= metL(funzione);

end


%Validazione anulare per identificare se un rigonfiamento sia davvero una
%bolla

filtro_anelli=false(H,W); 
variabile=false(1,numel(radL));
[X,Y]=meshgrid(1:W,1:H); %due matrici di coordinate della stessa taglia dell'immagine
tolleranza_anello_fill= 0.35;%frazione di 1 nella banda anulare
tolleranza_interna_esterna=0.06;%differenza media LoG interno vs esterno


%ciclo for per la validazione dei cerchi trovati, attraverso lo
%sfruttamento di un controllo metrico grazie all'utilizzo del tool 2d di
%matlab 

for i=1:numel(radL)
    centramentox =centL(i,1);
    centramentoy =centL(i,2); 
    R =radL(i);
    raggio_in= 0.70*R;% bordo interno della banda
    raggio_out =1.30*R; % bordo esterno della banda

    D2 =(X-centramentox).^2 +(Y-centramentoy).^2;
    ringBand =(D2>= raggio_in^2) & (D2<=raggio_out^2) & maschera_roi;
    fraz=nnz(filtro_LoG & ringBand) /max(nnz(ringBand),1);

    dentro  =(D2 <= raggio_in^2) & maschera_roi;
    fuori= (D2 >= (1.4*R)^2) & (D2 <= (1.9*R)^2) & maschera_roi;
    mu_in= mean(score_nuovo(dentro), 'omitnan');
    mu_out= mean(score_nuovo(fuori), 'omitnan');
    delta_io= mu_in - mu_out;

    if fraz >= tolleranza_anello_fill && delta_io >=tolleranza_interna_esterna
        fattore_di_crescita = 0.22;                   
        Rg=R *(1 +fattore_di_crescita);           
        diskMask =(D2 <=(Rg).^2) & maschera_roi;   

        CCv =bwconncomp(diskMask);
        Sv= regionprops(CCv,'Area','Perimeter','Eccentricity');
        if ~isempty(Sv)
            Area =Sv(1).Area; 
            Perimetro= Sv(1).Perimeter;
            Eccentricita= Sv(1).Eccentricity;
            Volume = 4*pi*Area / max(Perimetro^2, eps);
            if Eccentricita <= ECC_MAX && Volume>= CirC_MIN
                filtro_anelli =filtro_anelli | diskMask;
                variabile(i)= true;
            end
        end
    end
end

% pulizia finale per area fisica
filtro_anelli =bwareaopen(filtro_anelli, minAreaPx);
filtro_anelli= filtro_anelli & ~bwareaopen(filtro_anelli, maxAreaPx);

% maschera finale (se Hough non trova nulla, si usa solo LoG validato)
risultato_finale= filtro_anelli;
if ~any(risultato_finale(:))
    CC0 =bwconncomp(filtro_LoG);
    S0= regionprops(CC0,'Area','Perimeter','Eccentricity');
    if ~isempty(S0)
        A0 =[S0.Area];
        P0 =[S0.Perimeter]; E0=[S0.Eccentricity];
        circ0= 4*pi*A0 ./max(P0.^2, eps);
        vari0= (A0>=minAreaPx & A0<=maxAreaPx & E0<=ECC_MAX & circ0>=CirC_MIN);
        L0= labelmatrix(CC0);
        risultato_finale =ismember(L0,find(vari0));;
    end
end

Circon_finali =bwconncomp(risultato_finale);

if Circon_finali.NumObjects==0

    fprintf('Nessuna bolla riconosciuta');

else

    Sf=regionprops(Circon_finali,'Area','EquivDiameter','Eccentricity','Perimeter','Centroid','BoundingBox');
    Apx= [Sf.Area];  
    Dpx =[Sf.EquivDiameter]; 
    P=[Sf.Perimeter]; 
    Perimetro=[Sf.Eccentricity];
    D_mm= Dpx *mmPerPx;
    circf= 4*pi*Apx ./max(P.^2, eps);

    fprintf('Bolle riconosciute (LoG %s + anello/Hough): %d\n', pol, numel(D_mm));
    fprintf('Diametro medio: %.2f mm | mediana: %.2f mm\n', mean(D_mm), median(D_mm));
    fprintf('Ecc. media: %.3f | Circ. media: %.3f\n', mean(Perimetro), mean(circf));
end


figure('Name','Confronto filtri - Bolle circolari');
subplot(2,3,1);
imshow(Grigia,[]);   
title(' Originale');

subplot(2,3,2); 
imshow(filtro_adattivo,[]); 
title('Flat + CLAHE');

subplot(2,3,3); 
imshow(filtro_texture,[]); 
title(' Bilaterale anti-texture');

subplot(2,3,4); 
imshow(score_nuovo,[]); 
title(sprintf('4) LoG multiscala [%s]',pol));

subplot(2,3,5);
imshow(filtro_LoG);    
title('5) Maschera LoG (thr+area+holes)');

%adattamento a scala causa applicazione maschera
dx = -20;   
dy = 0;  

risultato_finale = imtranslate(risultato_finale, [dx dy],'OutputView','same', 'FillValues', 0);          
subplot(2,3,6);
imshow(risultato_finale);
title('Anelli riconosciuti (bolle)');;

figure('Name','Overlay finale');
imshow(Grigia,[]);
hold on;

visboundaries(risultato_finale,'Color','r','LineWidth',1.8);
title('Bolle (pseudo-circonferenze su LoG)'); 
hold off;

toc

end
