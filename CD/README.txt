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
            Siehe https://github.com/qwerdenkerXD/prot-fin-dev/tree/v0.4-exp-filter_hashes/experiments/recog_with_fft/README.md
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

    Erstellung der Plots in Linux Shell mit Zielordner Supplement/results:
        Programmversionen:
            System: Ubuntu 20.04.6 LTS
            Shell: zsh 5.8
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
            # Dazu in:
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
