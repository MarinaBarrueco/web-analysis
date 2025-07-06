# Instucciones para usar la web

## En la terminal

- Navega a la carpeta donde esta el main.py
- Usa el comando streamlit run main.py

```

(DS) marinabarrueco@MacBook-Air-de-Marina ~ % cd Desktop/Web\ Analysis/streamlit_site/                    
(DS) marinabarrueco@MacBook-Air-de-Marina streamlit_site % ls
Cleaned_Peptides.csv
VolcanoPlot.pdf
Clustering
data
Library
gibbscluster-2.0
Upregulated_peptides.pep
main.py
Upregulated_peptides_Variable.pep
(DS) marinabarrueco@MacBook-Air-de-Marina streamlit_site % streamlit run main.py 

```

## En el codigo

- El main.py tiene las insturcciones para renderizar la web
      - todo lo que sea st.algo es un objeto que se va a renderizar, el resto es python code normal.
- En Library/lib.py estan las  funciones que usa el main.py
