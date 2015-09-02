## predictFromSequence

Some info on the project here.


Predict transcription (or any positional mark) from genome sequence by modelling the factors that can bind to DNA and add or remove states. 

The raw DNA sequence cannot be altered but a series of layers (binary vectors of same length of sequence) can have regions of any length switched (0 <-> 1). 

Some factors may also recognise patterns on the layers (e.g. regions in state 1) in addition to or instead of the underlying DNA sequence. 

The system models a series of very many sequential binding and modifying events. 

The ability of factors to bind changes through this series so that the number and positions of possible binding events for each factor change during the sequence.


### Notes

2015-09-02:  [predictFromSequence.devFunc.R] Wrote alternative to runLayerBinding() called runLayerBinding.fast() that does not apply mods in random order and re-calculate, it just does them all for a single factor at a time. This breaks the idea of sequential mods creating pattern matches for other factors. The point was to increase speed enough to see if the system could find any TSS reliably. I also widened the goalposts by extending the TSS to a promoter with 200bp upstream and 200bp downstream. I created a random set of factors and a 1-layer LayerSet (layerList) and ran optimiseFactorSet() [n.iter=1000, mut.rate=0.1, modsPerCycle=10000].  This quickly scored above 1% and then to 3% (c.f. ~ 0.1%, previously), which is not too surprising with the widening of the targets.  The best score achieved was 12% and the plot of scores (by iteration index) showed a step-change around iteration 700.  Due to the stochastic nature of runLayerBinding.fast(), the set of factors does not produce the same result each time (range 3% - 12%). I therefore re-ran runLayerBinding.fast() 100 times without modifying the factor set, counting for each base the number of times it was marked. This gave a very nice set of peaks and troughs, some of the peaks being directly below TSS/promoters. The cor score for the summed vector and the tss vector was 12%. Using a threshold to divide the summed scores into bound/unbound did not help the cor score, possibly because the zeros are informative - in the test vector (HOXA cluster), one obvious pattern is the absense of prediction in the large gene desert to the left (chromosomal view) to the left of the HOXA cluster.

LayerSet became encapsulated within LayerList (to store extra annotations) and FactorSet needs to be encapsulated to carry the optimScores table (and in future other data) returned by optimiseFactorSet() TODO. 

Independently, I was looking at known motifs (loadMotifLibrary.XK2005.R). It would be better to see notes in context but here they are anyway:-
- the random patterns are longer (mean 12.24 vs 9.9) but often have 1000s of hits in the test_sequence
- all the real-life (but 'discovered') motifs have <= 130 hits and most have fewer than 20
- Information content? The real motifs are mostly more precise, with few Ns, and redundant codes.
- Could specify the letter freqeencies to match? 
- So, this selection of 'real' motifs are far more specific than my letter-frequency and length matches 'random' sample.

### Challenges

Speed - the number of modifications 

Random vs. Real - Using 'random' binding factors and evolving the system will be (IS!) extremely slow. The searchable space is just too huge. 

### TODO list.
Specify factorSet to include optimscores (e.g. make a list like layerList)

Spilt off calculation of binding positions from binding mods to speed up a series of mods.  Many TFs bind at independent locations.  Only a problem for binding events that depend on previous events but genome is so large and number of event so large, that many could be allowed to happen without re-calculation. e.g. if running 100k mods, could recalculate after every 1000. Practically, could split off binding calcs into separate thread. ASSUMPTION: # of events far greater than number of types of factors. [speed-up|design]

Additionally, create alternative algorithm ('FAST') that changes a high proportion of possible sites (-> 1.0) with every iteration. BUT many compositional factors will match all along the sequence or at many thousands of overlapping locations. These would then create many thousands of redundant hits to be processed, and also may swamp other signals. 

If optim starts to work. Track when factors get added to factor set. Also need way of measuring importance of each factor.