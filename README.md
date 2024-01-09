# MCFProject
In questa repository sono presenti i programmi del progetto di Metodi Computazionali per la Fisica.

Per scaricare la repository spostarsi in una cartella opportuna, non di pertinenza di altre repository git, ed eseguire da terminale il comando " git clone https://github.com/mariaritagrossud/MCFProject.git "

Nella cartella "Progetto" della repository scaricata ci sono 4 file principali:
- Il file "func.py" in cui sono definite le funzioni chiamate nei programmi 
- Il file "starclass.py" in cui si definisce la classe usata per gestire la simulazione di più stelle
- Il file  "Simulazioni.py" nel quale è eseguita la parte di simulazione Monte Carlo richiesta con le stelle assegnate
- Il file "Analisi_dati.py" in cui è svolta l'analisi dei dati forniti nel file "observed_starX.csv", anche esso presente nella cartella Progetto.

Poiché i moduli func.py, starclass.py e observed_starX.csv sono dentro la stessa cartella dei file che li usano, non è necessario specificarne il pathname. 
Se per qualche motivo essi fossere spostati in una cartella a parte, per importare correttamente il modulo func.py o starclass.py bisogna specificarne il pathname nel comando "sys.path.append(" pathname ")" presente nei due file mentre per caricare i dati di observed_starX.csv nel file "Analisi_dati" bisogna specificare il pathname, anziché il nome del file, nel comando "data = pd.read_csv("pathname", sep=",")".

Il file "Simulazioni.py" è gestito  tramite il modulo argparse per scegliere di quale stella avviare la simulazione. I comandi per eseguire la simulazione sono:
- "python3 Simulazione.py -sole" per avviare la simulazione considerando come stella il Sole
- "python3 Simulazione.py -betelgeuse" per avviare la simulazione considerando come stella Betelgeuse
- "python3 Simulazione.py -bellatrix" per avviare la simulazione considerando come stella Bellatrix
- "python3 Simulazione.py -alfacrucis" per avviare la simulazione considerando come stella Alfa Crucis
- "python3 Simulazione.py -tutte" se si vuole visualizzare la simulazione di, una dopo l'altra, tutte le stelle.
