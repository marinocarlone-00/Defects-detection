
function graffi(isColor)

clc
close all
tic;

%% INSERISCO E RICHIAMO L'IMMAGINE

% Percorso immagine
imagine_percorso='/Users/marinocarlone/Desktop/primo esempio/primafoto.png';

% Calibrazione della fotocamera con aggiunta del file ottenuto dall'add_on
CalibrazioneFile= mfilename('fullpath');
CalibrazioneCartella= fileparts(CalibrazioneFile);
calib_percorso =fullfile(CalibrazioneCartella, 'camera_calib.mat');

%Controllo esistenza del file
if exist(calib_percorso, 'file')

    calib_per =calib_percorso;

else

  fprintf('NoN  si è ancora caricato il file matlab utile per la calibrazione!! ')

end


%% METTO A SCALA 
%Utilizzo del file per la valutazione del rapporto pixel/mm/cm
S=load('camera_calib.mat');
cam =[];
cam_nome= '';

%qui abbiamo inserito un controllo per l'eventuale assenza del file
for fn=fieldnames(S)'
    v =S.(fn{1});

    if isa(v,'cameraParameters')
        cam= v; 
        cam_nome= fn{1};
        break
    end

end

assert(~isempty(cam), 'Nel .mat non c''è nessun cameraParameters ');
hasCam=true;   

%Conversione in millimetri
unita ="mm";
if isprop(cam,'WorldUnits') &&~isempty(cam.WorldUnits)

    unita =string(cam.WorldUnits);
end

switch lower(unita)
    
    case {"mm","millimetro","millimetri"} 
        scaleToMM= 1;
    case {"cm","centimetro","centimetri"}
        scaleToMM= 10;
    case {"m","metro","metri"}
        scaleToMM= 1000;

    otherwise

        scaleToMM=1;
end

%Qui abbiamo richiamato le dimensioni della scacchiera passati nella
%calibrazione
WP=cam.WorldPoints;                          
if size(WP,2)== 2
    WP3 =[WP, zeros(size(WP,1),1)];         
else
    WP3 =WP;                                 
end
assert(~isempty(WP), 'WorldPoints mancanti nel cameraParameters ');
D=squareform(pdist(WP));          
DimensioneQuadraBordi= min(D(D>0));            
DimensioniQuadraConFattore= DimensioneQuadraBordi *scaleToMM;;
 [ii,jj] =find( D>0 & abs(D -DimensioneQuadraBordi)<= 2e-3*DimensioneQuadraBordi & triu(true(size(D)),1) );%abbiamo imposto una tolleranza di riferimento
assert(~isempty(ii), ' Non trovo coppie adiacenti nei WorldPoints ');

%Si stima la distanza in Pixel: abbiamo utilizzato lo Statistic toolbox di matlab 
pix_Distanza=[];
Pose=size(cam.RotationVectors,1);

for k =1:Pose
    r_vec= cam.RotationVectors(k,:);
    t_vec= cam.TranslationVectors(k,:);

    try
        Primaprova=worldToImage(cam,r_vec,t_vec,WP3);                 

    catch

        Rotazione=rotationVectorToMatrix(r_vec);
        Primaprova=worldToImage(cam, Rotazione,t_vec, WP3);                    

    end

    if any(~isfinite(Primaprova(:)))
        continue; 
    end

    d_indice=hypot(Primaprova(ii,1)-Primaprova(jj,1),Primaprova(ii,2)-Primaprova(jj,2));
    pix_Distanza =[pix_Distanza; d_indice]; 
end
assert(~isempty(pix_Distanza), 'Nessuna distanza in pixel è valida ');

%Metto a Scala: calcolo la mediana
meanPixPerEdge=median(pix_Distanza);   
pxPerMM= meanPixPerEdge/DimensioniQuadraConFattore;
mmPerPx =1/pxPerMM;

fprintf('Scala: %.6f px/mm (%.6f mm/px) | lato≈%.3f mm',pxPerMM, mmPerPx);

% Caricamento immagine con parametri calibrati
I=imread(imagine_percorso);
colori_flag=false;

% Scala di grigi
if size(I,3)== 3
    Immgrigia= rgb2gray(I);
    colori_flag= true;
else
    Immgrigia =I;
end

% Adattamento a scala della foto passata in input
Dimensioni1 =cam.ImageSize;   
DimensioniAdattate =size(Immgrigia);               
pxPerMM_cal= pxPerMM;

if all(DimensioniAdattate(1:2)==Dimensioni1(1:2))

    NoDist=undistortImage(Immgrigia, cam);
    pxPerMM=pxPerMM_cal;            
    mmPerPx= 1/pxPerMM;

else
    
    
    sx=DimensioniAdattate(2)/Dimensioni1(2);
    sy=DimensioniAdattate(1)/Dimensioni1(1);
    media=mean([sx sy]);

    %No distorsioni 
    NoDist=Immgrigia;
    pxPerMM =pxPerMM_cal*media;;
    mmPerPx= 1/pxPerMM;

end

immagineGrigia=im2double(NoDist);


%% RILEVAMENTO GRAFFI

%Parametri 
LARGH_MM=[0.5 1.5];% spessore graffio 

%Pixel per la calibrazione
min_px= max(1,round(LARGH_MM(1) /mmPerPx));
max_px=max(min_px+1,round(LARGH_MM(2)/mmPerPx));
pxPerMM= 1/mmPerPx;

fprintf('mmPerPx=%.4f | pxPerMM=%.2f', mmPerPx, pxPerMM);

% Conversione dell'immagine
immagineGrigia=im2double(immagineGrigia);

if colori_flag==true

    % Appiattimento e contrasto con filtri (pipeline di filtri)
    fprintf('\n Stiamo analizzando una immagine colorata')

    %questo filtro appiattisce l'illuminazione, lasciando le strutture ad
    %alta frequenza, eliminando le non uniformità di luce ed eventuali
    %vignettature
    addatt =max(15, round(6*max_px));%scelta del raggio maggiore della dimensione massima dell'imperfezione cosi da poter non assorbire il difetto nel processo di appiattimento
    Filtro_appiattimento=imflatfield(immagineGrigia, addatt);

    %Filtro per la equalizzazione attiva, attraverso una griglia 8X8 ogni
    %tile viene interpolato aumentando il contrasto locali dei dettagli
    %utili senza bruciare i bianchi (forma di equalizzazione locale) e
    %permettendo di applicare un processo di thresholding successivo
    Filtro_Contrasto=adapthisteq(Filtro_appiattimento,'NumTiles',[8 8],'ClipLimit',0.01);

    % Risposta multi-scala + multi-orientazione
    risp=zeros(size(Filtro_Contrasto));

    % Black-hat con dischi: tipo di filtro morfologico mette in evidenza le
    % strutture più scure come graffi o tagli a diversi raggi 'r', filtro
    % morfologico passa alto, sensibile ai dettagli scuri rispetto al fondo
    Set =unique(max(3, round([min_px, (min_px+max_px)/2,max_px]/2)));
    for r= Set

        risp= max(risp, imbothat(Filtro_Contrasto,strel('disk', r)));

    end

    % Black-hat con linee orientate, stesso procedura precedente ma con
    % linne orientate rispetto al riferimento, filtro DIREZIONALE
    % mofologico
    angoloSet= 0:15:165;
    lineaLun= max(5, round(2.5*max_px));

    for t =angoloSet
        risp=max(risp, imbothat(Filtro_Contrasto,strel('line',lineaLun,t)));
    end

    %Filtro del tool del pacchetto image preocessing matlab, esalta strutture filamentose scure
    %complementando i filtri morfologici rendendo più robusta la risposta a
    %graffi fini
    if exist('fibermetric','file')==2
        fm=fibermetric(Filtro_Contrasto, 'ObjectPolarity','dark','StructureSensitivity', max(1,round(0.7*max_px)));
        risp =max(risp,fm);
    end

    % Normalizzazione e compressione dinamica: i graffi scuri diventano
    % bianchi brillanti su fondo scuro e quindi è più facile sogliarli
    score_norm= mat2gray(risp).^0.5;

    % Mostra immagine originale
    figure('Name','Immagine originale');
    imshow(immagineGrigia); 
    title('Immagine originale');

    % Mostra Immagine Normalizzata 
    figure('Name','Score (graffi bianchi)'); 
    imshow(score_norm); 
    title('Filtro score');

    if colori_flag == true

        % Estrazione con soglia fissa rispetto alla binarizzazione
        soglia =0.79;
        bwGraffi= score_norm>soglia;

    end

    % Rimozione chiazze attraverso l'image processing tool di Matlab,
    % morfologia di pulizia, 
    bwGraffi =bwareaopen(bwGraffi,105);

    % Visualizzazione graffi
    figure('Name','Graffi estratti'); 
    imshow(bwGraffi);
    title(sprintf('Graffi estratti (score >%.2f)', soglia));


    % Sovrapposizione dei graffi su immagine reale
    figure('Name','Overlay su originale');
    imshow(immagineGrigia);
    hold on;;
    visboundaries(bwGraffi, 'Color','r', 'LineWidth',1.5);
    title(sprintf('Graffi rilevati (threshold score> %.2f)',soglia));
    hold off;

    %PIPELINE COMPLETA
    %Normalizzazione + potenziare  i  dettagli deboli (mat2gray.^0.5),
    %Visualizzazione,
    %Applicazione della soglia fissa (graffi diventano bianchi),
    %Fare pulizia morfologica (bwareaopen),
    %Immagine binaria finale dei graffi.


else

    fprintf('Stiamo analizzando una immagine bianco e nero')
    %il filtro seguente è: filtro di equalizzazione adattiva, si divide
    %l'immagine in blocchi, si equalizza ciascun tile e si ricompone
    %l'immagine in blocchi per evitare discontinuità, il ClipLimit a 0.01
    %evita un contrasto troppo esasperato, migliora il contraso locale
    %sopratutto in casi di illuminazione non uniforme
    Filtro_Contrasto=adapthisteq(immagineGrigia, 'NumTiles',[8 8], 'ClipLimit',0.01);

    % Applica soglia all'immagine binarizzata
    soglia =0.01;  
    soglia2= 0.20;

    %Filtro del tool pre-processing: rimozione chiazze attraverso 
    % l'image processing tool di Matlab,
    % morfologia di pulizia, 
    bwGraffi= (Filtro_Contrasto >soglia) & (Filtro_Contrasto< soglia2);   
    bwGraffi= bwareaopen(bwGraffi, 195);%
    
    % Mostra immagine originale
    figure('Name','Immagine originale');
    imshow(immagineGrigia);
    title('Immagine originale');

    figure('Name','Graffi estratti');
    imshow(bwGraffi);
    title(sprintf('Graffi estratti (score< %.2f)',soglia));

    figure('Name','Overlay su originale');
    imshow(immagineGrigia);
    hold on;;
    
    visboundaries(bwGraffi,'Color','r', 'LineWidth',1.5);
    title(sprintf('Graffi rilevati (threshold score> %.2f)', soglia));
    hold off;

end


% Analisi quantitativa dei graffi 

% Etichetta ogni graffio individuato
primaetichetta=bwconncomp(bwGraffi);
vari=regionprops(primaetichetta, 'Area', 'MajorAxisLength', 'MinorAxisLength');

if isempty(vari)
    warning('Nessun graffio rilevato.');
else

    % Estrazione misure
    lunghezze_cm =[vari.MajorAxisLength]*mmPerPx;
    larghezze_mm =[vari.MinorAxisLength]*mmPerPx;
    media_lunghezza_cm= mean(lunghezze_cm);
    media_larghe_mm= mean(larghezze_mm);

    % Stampa risultati
    fprintf('\n  RISULTATI GRAFFI  \n');
    fprintf('Lunghezza media: %.2f cm\n', media_lunghezza_cm);
    fprintf('larghezza media: %.2f mm\n',media_larghe_mm);
end

toc
end
    
