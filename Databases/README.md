# Databases

Two databases were used in the analysis. Note that the databases represent the content of RefSeq or the NCBI nt sequence set at the point in time when they were downloaded.


std database built on Oct 27, 2017. Command line:
```
make std
```

Content:
```
Name	Rank	Taxonomy ID	# of sequences	# of k-mers	k-mer duplication
Bacteria	superkingdom	2	15382	1.35E+10	2.27
Homo sapiens	species	9606	232355	2.67E+09	2.49
Archaea	superkingdom	2157	369	4.75E+08	1.26
Viruses	superkingdom	10239	139477	2.94E+08	4.47
synthetic construct	species	32630	9740	2.39E+06	3.91
```

nt database built February 23, 2018. Command line:
```
make nt
```

Content:
```
Name	Rank	Taxonomy ID	# of sequences	# of k-mers	k-mer duplication
Bacteria	superkingdom	2	6549136	1.72E+10	2.67
Fungi	kingdom	4751	3836040	4.85E+09	1.76
Alveolata	no rank	33630	472264	8.47E+08	1.68
Apicomplexa	phylum	5794	299510	6.97E+08	1.77
Archaea	superkingdom	2157	326186	6.64E+08	1.38
Viruses	superkingdom	10239	2050013	5.42E+08	6.49
Euglenozoa	no rank	33682	150499	3.19E+08	2.09
Kinetoplastida	order	5653	146387	3.14E+08	2.1
Stramenopiles	no rank	33634	330637	3.12E+08	1.47
Chlorophyta	phylum	3041	159680	2.35E+08	1.2
Amoebozoa	no rank	554915	133625	1.45E+08	1.23
Rhodophyta	no rank	2763	60582	6.93E+07	1.52
Choanoflagellida	order	28009	21662	4.46E+07	1.08
Parabasalia	no rank	5719	57906	3.74E+07	1.47
Entamoeba	genus	5758	41528	3.45E+07	1.44
Cryptophyta	class	3027	29269	3.40E+07	1.09
Haptophyceae	no rank	2830	41463	3.08E+07	1.44
Heterolobosea	class	5752	16885	2.40E+07	1.07
Apusozoa	no rank	554296	10877	2.07E+07	1.06
Fornicata	no rank	207245	14657	9.35E+06	1.45
Rhizaria	no rank	543769	19683	8.83E+06	2.07
Jakobida	no rank	556282	459	1.39E+06	1.09
Glaucocystophyceae	class	38254	352	6.69E+05	1.16
Syndiniales	order	88547	2317	3.39E+05	2.12
Oxymonadida	order	66288	472	2.35E+05	2.24
Malawimonadidae	family	136087	43	1.58E+05	1.04
Centroheliozoa	no rank	193537	295	1.34E+05	1.59
Telonemida	order	589438	177	8.36E+04	1.21
Palpitomonas	genus	759891	14	8.20E+04	1.06
Collodictyonidae	family	190322	22	7.55E+04	1.52
Picozoa	phylum	419944	104	6.48E+04	1.21
Tsukubamonadidae	family	1084709	9	5.54E+04	1
Katablepharidophyta	class	339960	248	4.57E+04	2.36
Breviatea	no rank	1401294	21	1.92E+04	1.03
Trimastix	genus	137418	13	1.29E+04	1.03
```