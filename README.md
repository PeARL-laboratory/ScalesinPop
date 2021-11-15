# ScalesinPop

This Github repository includes the code to produce LDA topic models for the triadic corpora reported in: 

Sears & Glosson (2021), "Topic models reveal scale systems in popular music." Paper presented as a late-breaking demo at the International Society for Music Information Retrieval Conference (2021). https://archives.ismir.net/ismir2021/latebreaking/000048.pdf

## Code

After cloning/pulling the repository, run LDA_cp.m to produce the LDA models for the common-practice data sets.

Dependencies: 
RNs_2_SDs.m (included)
MidiToolbox https://github.com/miditoolbox/1.1 (download and install)

Users may then run LDA_pop to produce the LDA models for the popular music data sets.

Note that these models are non-deterministic, so the data will not exactly reproduce the results reported in Sears & Glosson (2021).

## Data Sets
The paper includes the following data sets:

Annotated Beethoven Corpus (ABC): https://github.com/DCMLab/ABC

Bach Chorales Melody-Harmony Corpus (BCMH): https://github.com/PeARL-laboratory/BCMH

Theme and Variation Encodings with Roman Numerals (TAVERN): https://github.com/jcdevaney/TAVERN

McGill Billboard: https://ddmal.music.mcgill.ca/research/The_McGill_Billboard_Project_(Chord_Analysis_Dataset)/

RollingStone-200: http://rockcorpus.midside.com/

Further details about the relative root encoding scheme that was used to parse these data sets appears in Sears & Forrest (2021): https://osf.io/kdzm3/

