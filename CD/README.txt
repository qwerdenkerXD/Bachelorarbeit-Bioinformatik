Anhang der Bachelorarbeit von Franz-Eric Sill, Abgabe am 26.08.2024, Matrikelnummer 2139315, Technische Hochschule Bingen

Die Bachelorarbeit befindet sich in Bachelorarbeit.pdf

Ordnerstruktur:
    Siehe Bachelorarbeit.pdf in Abschnitt "Anhang"

Implementierung:
    Algorithmus 1:
        https://github.com/qwerdenkerXD/prot-fin-dev/blob/v0.4/experiments/recog_with_fft/methods/actions/algorithm/kidera.py
        in function "get_aa_vector".
    Algorithmus 2:
        https://github.com/qwerdenkerXD/prot-fin-dev/blob/v0.4/experiments/recog_with_fft/methods/actions/algorithm/constellation.py
        in function "create_constellation".

    In Version 0.4 fehlt noch die Implementierung, die alle Kidera Faktoren einbezieht, weshalb die folgenden Algorithmen dahingehend leicht abweichen.
    Die Implementierung wurde mit Experiment 1, dem UniRef90 Sampling hinzugefügt, auf welchem die anderen Experimente aufbauen, da sie die Ergebnisse mindestens eines Signifikanzniveaus einbeziehen, und wurde daher für die Beschreibung des Algorithmus der Arbeit mit aufgegriffen, zumal es auch die Intention von prot-fin ist, alle Kidera Faktoren zu betrachten.

    Algorithmus 3:
        https://github.com/qwerdenkerXD/prot-fin-dev/blob/v0.4/experiments/recog_with_fft/methods/actions/algorithm/hash_gen.py
        in function "create_hashes".
        In Experiment 1: https://github.com/qwerdenkerXD/prot-fin-dev/blob/17e7482d7571dd0c2db433fbd644f81197dbb5a5/experiments/recog_with_fft/methods/actions/algorithm/hash_gen.py
    Algorithmus 4:
        https://github.com/qwerdenkerXD/prot-fin-dev/blob/v0.4/experiments/recog_with_fft/methods/actions/create_db.py
        in function "create_db".
        In Experiment 1 unterscheidet sich dieser Algorithmus lediglich in der function "hashes_from_seq", hier: https://github.com/qwerdenkerXD/prot-fin-dev/blob/17e7482d7571dd0c2db433fbd644f81197dbb5a5/experiments/recog_with_fft/methods/actions/algorithm/__init__.py

Experiment-Code:
    Experiment 1:
        https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-uniref_sampling/experiments/recog_with_fft
    Experiment 2:
        https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-filter_hashes/experiments/recog_with_fft
    Experiment 3:
        https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-target_zone/experiments/recog_with_fft
    Experiment 4:
        https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-selection_method/experiments/recog_with_fft

Reproduzierbarkeit:
    Siehe Unterabschnitt "Reproduce" in Abschnitt "Results" im jeweiligen README:
        Dateien in Supplement/material/Filter Hashes:
            https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-filter_hashes/experiments/recog_with_fft/README.md
        Dateien in Supplement/material/previous:
            https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.3-exp-stft_params/experiments/recog_with_fft/README.md
        Dateien in Supplement/material/Selection-Method:
            https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-selection_method/experiments/recog_with_fft/README.md
        Dateien in Supplement/material/Target-Zone:
            https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-target_zone/experiments/recog_with_fft/README.md
        Dateien in Supplement/material/UniRef90 Sampling:
            https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-uniref_sampling/experiments/recog_with_fft/README.md
    
    Siehe Abschnitt "Generate table of Kidera Factors" im README:
        Supplement/material/Amino_Acid_Kidera_Factors.csv:
            https://github.com/qwerdenkerXD/prot-fin-dev/blob/v0.4/README.md

    Erstellung der Plots in Linux Shell mit Zielordner Supplement/results mit der Linux Shell in dem Ordner dieses README.txt geöffnet:
        Programmversionen:
            System: Ubuntu 20.04.6 LTS
            Shell: zsh 5.8
            git: 2.25.1
            python3: 3.8.10
            scipy: 1.10.1
            numpy: 1.23.0
            pandas: 2.0.1
            tqdm: 4.66.2
            matplotlib: 3.5.2
            R: 3.6.3
            tibble: 3.2.1
            ggplot2: 3.5.0
            ggdist: 3.3.2
            dplyr: 1.1.4

        Abbildung 3:
            cd Supplement/methods
            Rscript ../methods/plot_sequence_lengths.R protein.fa ../results/seq_lens.png
        Abbildungen 1, 2, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15:
            cd Supplement/methods
            python3 ../methods/plot_sequence_lengths.R protein.fa ../results/seq_lens.png
        Abbildung 6:
            # Dies ist eine einfache LaTeX-Tabelle.
            # Um den LaTeX Code dazu zu sehen, in:
            # https://github.com/qwerdenkerXD/Bachelorarbeit-Bioinformatik/blob/master/Bachelorarbeit/Material_und_Methoden.tex
            # nach "\label{fig:hash}" suchen
        Abbildung 8a, 8b, 8c, 8d:
            cddir=$(pwd)
            git clone https://github.com/qwerdenkerXD/prot-fin-dev.git
            cd prot-fin-dev
            git checkout v0.4-exp-uniref_sampling
            cd experiments/recog_with_fft
            for i in 5 0.1 0.01 0.001; do;
                python3 evaluation.py plot-frequencies "${cddir}/Supplement/material/protein.fa" "${cddir}/Supplement/results/frequencies_${i}.png"  # eventuell mit -c oder --cpu die Anzahl CPU-Cores erhöhen, sofern möglich, Standard ist 1, bei der Erstellung der 4 Plots verwendet wurden 80 CPU-Cores
            done
        Abbildung 8e:
            cddir=$(pwd)
            git clone https://github.com/qwerdenkerXD/prot-fin-dev.git
            cd prot-fin-dev
            git checkout v0.3-exp-stft_params
            cd experiments/recog_with_fft
            python3 evaluation.py plot-frequencies "${cddir}/Supplement/material/protein.fa" "${cddir}/Supplement/results/8e.png"  # eventuell mit -c oder --cpu die Anzahl CPU-Cores erhöhen, sofern möglich, Standard ist 1, bei der Erstellung des Plots verwendet wurden 80 CPU-Cores

Dateiformate: (Bei Verständnisfragen von Spalten, siehe Bachelorarbeit.pdf)
    Comma Separated Values (CSV):
        Supplement/material/*/runtimes_createdb.csv:
            Spalten:
                Winsize -> Fenstergröße der STFT
                n_Peaks -> n_peaks, siehe Bachelorarbeit.pdf in Unterabschnitt "Durchführung"
                Overlap -> Überlappung der STFT Fenster
                Time -> Dauer der Datenbankerstellung für die genannten Parameter
                DB_Size -> Datenbankgröße für die genannten Parameter in Megabyte
            Zusätzliche Spalten:
                Supplement/Target-Zone/runtimes_createdb.csv:
                    Target_Zone -> Größe der Target-Zone der jeweiligen Datenbank
            Zusatzinformation:
                In Experiment 2 werden die Hashes im Matching-Modus vor der Iteration durch die Eingabesequenzen gefiltert, also nicht bei der Datenbankerstellung. Das spart Zeit, weil dasselbe Datenbank-Training nicht mehrfach durchgeführt werden muss. Deshalb gibt es für die Datenbank keinen Parameter für das Filterquantil in der CSV-Datei.
        Supplement/material/*/summary.csv:
            Enthalten zusammengefasste Ergebnisse vom Single-Protein-Matching.

            Mögliche Spalten:
                Window_Size -> Fenstergröße der STFT
                N_Peaks -> n_peaks, siehe Bachelorarbeit.pdf in Unterabschnitt "Durchführung"
                Overlap -> Überlappung der STFT Fenster
                Median_Hits -> mediale Anzahl gefundener Treffer
                Average_Hits -> mittlere Anzahl gefundener Treffer
                Self_Matches -> Eingabeprotein hatte besten Score (True) oder nicht (False)
                Unique_Self_Matches -> Kein anderer Treffer als das Eingabeprotein hatte den besten Score (True)
                Mean_F1_Score -> siehe Gleichung 2 in Bachelorarbeit.pdf, nur dass FP auch einen schlechteren Score als das letzte Familienmitglied haben dürfen
                Mean_Precision -> siehe Gleichung 2 in Bachelorarbeit.pdf, nur dass FP auch einen schlechteren Score als das letzte Familienmitglied haben dürfen
                Mean_Liberal_F1_Score -> siehe Gleichung 2 in Bachelorarbeit.pdf
                Mean_Liberal_Precision -> siehe Gleichung 2 in Bachelorarbeit.pdf
                Mean_Sharpness -> siehe Gleichung 2 in Bachelorarbeit.pdf
                DB_Size_MB -> Datenbankgröße in Megabyte
                Runtime_find-matches -> Laufzeit des Single-Protein-Matchings in (Stunden:)Minuten:Sekunden
                Hash_Use_First_Appearance -> wird die erste Position eines Hashes gespeichert (True) oder die letzte (False)
                Quantile -> Filterquantil, welche Hashes behalten werden
                Skip_First_K_Freqs -> Parameter k aus Experiment 4
                Selection_Method -> Selektionsansatz aus Experiment 4:
                    1. none
                    2. absolute
                    3. deviation
                Significance -> Signifikanzniveau der UniRef90 Sampling Amplituden
                Median_Hits_Top_Score -> mediale Anzahl Treffer mit bestem Score
                Runtime_create-db -> Laufzeit der Datenbankerstellung in (Stunden:)Minuten:Sekunden
                Target_Zone -> Größe der Target-Zone

        Supplement/material/*/summary_match_family.csv:
            Enthalten zusammengefasste Ergebnisse vom Family-Matching.

            Mögliche Spalten:
                Window_Size -> Fenstergröße der STFT
                N_Peaks -> n_peaks, siehe Bachelorarbeit.pdf in Unterabschnitt "Durchführung"
                Overlap -> Überlappung der STFT Fenster
                Hash_Use_First_Appearance -> wird die erste Position eines Hashes gespeichert (True) oder die letzte (False)
                Quantile -> Filterquantil, welche Hashes behalten werden
                Average_F_Score -> siehe Gleichung 2 in Bachelorarbeit.pdf
                Average_Precision -> siehe Gleichung 2 in Bachelorarbeit.pdf
                Average_Sharpness -> siehe Gleichung 2 in Bachelorarbeit.pdf
                avg(Family_Members/Match_Count)
                Average_Hash_Intersection_Size -> Mittlere Anzahl Hashes, die sich eine Familie teilt
                Runtime_match-family -> Laufzeit des Family-Matchings in (Stunden:)Minuten:Sekunden
                Skip_First_K_Freqs -> Parameter k aus Experiment 4
                Selection_Method -> Selektionsansatz aus Experiment 4:
                    1. none
                    2. absolute
                    3. deviation
                Significance -> Signifikanzniveau der UniRef90 Sampling Amplituden

        Supplement/material/UniRef90 Sampling/sample_WINSIZE_*.csv:
            Enthält die Sampling-Ergebnisse von Experiment 1. WINSIZE steht für die Fenstergröße, mit der Samples gewählt wurden.

            Spalten:
                Window_Size -> Fenstergröße
                Sample_Count_Mio -> Anzahl Samples in abgerundeten Millionen
                Frequency -> STFT Frequenz-Nummer, siehe Caption Abbildung 6 in Bachelorarbeit.pdf
                First_Lower_Outlier -> Erster unterer Outlier der Amplituden, wenn vorhanden
                First_Upper_Outlier -> Erster oberer Outlier der Amplituden, wenn vorhanden
                Std_Deviation -> Standardabweichung der Amplituden
                Mean -> Mittelwert der Amplituden
                x%%-Quantile -> Das x%-Quantil der Amplituden

        Supplement/material/Amino_Acid_Kidera_Factors.csv:
            Spalten:
                Kidera_Factor -> Bezeichner des Kidera Faktors
                Kidera_Factor_Description -> Beschreibung des Kidera Faktors (wie in Tabelle 1 in Bachelorarbeit.pdf)
                A-Z ohne B,J,O,U,Z -> Wert des Kidera Faktors für die entsprechende Aminosäure

    Tab Separated Values (TSV):
        Supplement/material/mapman.tsv:
            Spalten:
                BINCODE -> MapMan-Bin
                NAME -> MapMan-Bin mit menschenlesbaren Bezeichnern statt Nummern
                IDENTIFIER -> Protein ID, wenn Blattknoten
                DESCRIPTION -> Beschreibung der Protein-Funktion, wenn Blattknoten
                TYPE -> T, wenn Blattknoten

    Waveform Audio (WAV):
        Supplement/A0A1D8EJF9.wav:
            Aminosäuresequenz mit Protein-ID "A0A1D8EJF9" aus Supplement/material/protein.fa als Audio.
            Siehe function "seq_to_audio" und "get_aa_chords" in:
                https://github.com/qwerdenkerXD/prot-fin-dev/blob/134663504c1e21fe579b798cd41c016ad81e1284/experiments/aa_as_chords/methods/protfin.py

    FASTA:
        Supplement/material/protein.fa