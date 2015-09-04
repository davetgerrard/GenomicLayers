## predictFromSequence


Predict transcription (or any positional mark) from genome sequence by modelling the factors that can bind to DNA and add or remove states. 

The raw DNA sequence cannot be altered but a series of layers (binary vectors of same length of sequence) can have regions of any length switched (0 <-> 1). 

Some factors may also recognise patterns on the layers (e.g. regions in state 1) in addition to or instead of the underlying DNA sequence. 

The system models a series of very many sequential binding and modifying events. 

The ability of factors to bind changes through this series so that the number and positions of possible binding events for each factor change during the sequence.


### Notes (reverse chronological)
__2015-09-04__: Downloaded HOXD sequence and GENCODE genes in region (chr2:176893697-177091402) from UCSC (table browser for latter).
I ran runLayerBinding() ten times using the original and final factorSets for the 5_layer set run.
	> cor(tss.vector, result.vector.initial/reps)
	[1] 0.03020895
	> cor(tss.vector, result.vector.final/reps)
	[1] 0.09556068
It seems there was an improvement. Currently quite slow to re-run runLayerBinding() on my machine over the 5_layer set. 

I don't really expect these optimised sets to be applicable to every sequence (esp. if strand specific) but am hoping they are recognising promoters and not just optimising for the coverage needed to hit as many promoters as possible with the 'right' amount of coverage. The large difference between the 1-layer and 5-layer tests I think shows that it is not just optimising for coverage. 

I wonder if can simulate a maximal coverage finder?

- randomise the promoter positions
- cover different amounts of the target sequence
- what is the maximum correlation that I can get? 

P.S. the final score from yesterday's run was 0.24! (N.B. it successfully self-exited after 1000 iterations without improvement)

__2015-09-03__: Ran first script on the cluster for 10,000 iterations. No better than the first 1000 run a day earlier. May have been a bit of luck involved there. Still progressed up from 4% (iter 20) to ~9.5%  (although most achieved in first 1000 iters). Mostly just a requirement test for such a script. Parsing the resulting factors is still not very easy, may need better summarisation scripts. Also noted that many of the factors were to bind anywhere and set the target layer to 0 - not a bad strategy when the target size is a relatively small proportion of the full sequence. 

Also started a 5-layer test today (10,000 iters). Despite expectation that this might slow the progression, early inications are that it is already scoring high!

	[1] "Round 187 oldScore 0.135323991505638 newScore 0.173346165894045 Better?"
	[1] "Yes!"
	
Some of the non-improved scores are still very low, which might indicate the system is easily broken but more likely is due to the fact that scoring is still stochastic and is picking rounds that largely by chance were high. TODO remove randomness in layer modificiation (but how? saturation of all possible sites is only way I can think and would be VV. slow). OR mutate at lower rate and run modification several times to get average score. This second option would make the scores more robust but would slow the system too. The second option also would avoid saturation, which is good as eventually I intend to introduce factor abundance as a parameter - i.e. I intend that not all sites are occupied.

Looking at the factor set of the five layer case. It seems it is going to be very difficult to read what is happening in these traces. It might be a good idea to make a 'movie' version of runLayerBinding to watch the sequence of changes that occur.  With 5 layers, and 30 factors, could be tricky to visualise.

If this is 'working', need to start doing proper training and testing on different sets of sequences.

- sort out stochastisity of tests
- add in multi-sequence training and testing sets, with cross-validation.


__2015-09-02__:  [predictFromSequence.devFunc.R](scripts/predictFromSequence.devFunc.R) Wrote alternative to runLayerBinding() called runLayerBinding.fast() that does not apply mods in random order and re-calculate, it just does them all for a single factor at a time. This breaks the idea of sequential mods creating pattern matches for other factors. The point was to increase speed enough to see if the system could find any TSS reliably. I also widened the goalposts by extending the TSS to a promoter with 200bp upstream and 200bp downstream. I created a random set of factors and a 1-layer LayerSet (layerList) and ran optimiseFactorSet() [n.iter=1000, mut.rate=0.1, modsPerCycle=10000].  This quickly scored above 1% and then to 3% (c.f. ~ 0.1%, previously), which is not too surprising with the widening of the targets.  The best score achieved was 12% and the plot of scores (by iteration index) showed a step-change around iteration 700.  Due to the stochastic nature of runLayerBinding.fast(), the set of factors does not produce the same result each time (range 3% - 12%). I therefore re-ran runLayerBinding.fast() 100 times without modifying the factor set, counting for each base the number of times it was marked. This gave a very nice set of peaks and troughs, some of the peaks being directly below TSS/promoters. The cor score for the summed vector and the tss vector was 12%. Using a threshold to divide the summed scores into bound/unbound did not help the cor score, possibly because the zeros are informative - in the test vector (HOXA cluster), one obvious pattern is the absense of prediction in the large gene desert to the left (chromosomal view) to the left of the HOXA cluster.

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

Evolve a base set of factors that can predict basal TSS (does this distinction really exist?). Then evolve from this state new sets of factors that each predict tissue-specific sets of genes. Then re-combine with the basal set (if that is possible) or have sequential series. Hmmm. 

### TODO list.
Specify factorSet to include optimscores (e.g. make a list like layerList)

Spilt off calculation of binding positions from binding mods to speed up a series of mods.  Many TFs bind at independent locations.  Only a problem for binding events that depend on previous events but genome is so large and number of event so large, that many could be allowed to happen without re-calculation. e.g. if running 100k mods, could recalculate after every 1000. Practically, could split off binding calcs into separate thread. ASSUMPTION: # of events far greater than number of types of factors. [speed-up|design]

__DONE__ Additionally, create alternative algorithm ('FAST') that changes a high proportion of possible sites (-> 1.0) with every iteration. BUT many compositional factors will match all along the sequence or at many thousands of overlapping locations. These would then create many thousands of redundant hits to be processed, and also may swamp other signals. 

If optim starts to work. Track when factors get added to factor set. Also need way of measuring importance of each factor.

__DONE__ In script logs, log initial score.

Multi-sequence (random) training and testing sets. Concatenate seqs? - no. How to avoid over-fitting? Need fast way to get a bunch of sequence. May be worth using BSGenome and write functions to retrieve arbitrary segments. Would also then need to pull all relevant data for each region.

Write parallel version of optimiseFactorSet to use scores from a range of sequences (take mean or median score).

Are the factorSets strand specific? Could program to make searches on both strands.

