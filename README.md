## predictFromSequence


Predict transcription (or any positional mark) from genome sequence by modelling the factors that can bind to DNA and add or remove states. 

The raw DNA sequence cannot be altered but a series of layers (binary vectors of same length of sequence) can have regions of any length switched (0 <-> 1). 

Some factors may also recognise patterns on the layers (e.g. regions in state 1) in addition to or instead of the underlying DNA sequence. 

The system models a series of very many sequential binding and modifying events. 

The ability of factors to bind changes through this series so that the number and positions of possible binding events for each factor change during the sequence.

## TODO List

DateComplete|DateAdded|Priority|Description
: ------ :|:-------------:|: -----:|:--------------------------
TODO|2015-12-22|1|GitHub repo (requires open access?)
TODO|2015-12-22|1|CRAN/Bioconductor
TODO|2015-12-22|1|Implement factor abundance on a scale from 0 to max. Then mods can be applied by the proportion of each factor. 
TODO|2015-12-22|1|FactorSet quantities are mutable
TODO|2015-12-22|1|Factors can be selected from libraries (e.g. Jolma, TRANSFAC)
TODO|2015-12-22|1|Improve summarisation of factorSets (e.g. base composition, types)
TODO|2015-12-22|1|Allow factor type to mutate during optimisation
TODO|2015-12-22|1|Either re-write the optimiisation function to export all functions in the environment to cluster nodes OR package the code properly. OR add a customisable watch_function to the optimisation function.
TODO|2015-12-22|1|move on to testing against specific TSS being activated (quantitative?) or de-activated.
TODO|2015-12-22|1|no long range in-cis effect (could be added in by enabling offset parameter in bindingFactor)
TODO|2015-12-22|1|FactorSet needs to be encapsulated to carry the optimScores table (and in future other data) returned by optimiseFactorSet()
TODO|2015-12-22|1|Factors can be converted from common motifs formats (e.g. meme). Enables use of libraries. How to handle layer specificities?
TODO|2015-12-22|1|Allow duplication/incorporation of the same factor within a factorSet so that it can be applied more than once. WARNING: this will require linked abundance measure to be stored separately. 
TODO|2015-12-22|1|Test if system can rediscover some promoter motifs (e.g. HNF4A in liver specific promoters). Might be better to use existing data from similar papers.
TODO|2015-12-22|1| Run with many more factors to find more TSS? Or allow number of factors to mutate (start small, allow to grow within limits).
TODO|2015-12-22|1| Separate factor names from indexes within a factorSet, allow multiple factors to share names (and, hence abundance).  Hmmmm, design issue here. 
TODO|2015-12-22|1| Overhaul factorSet objects
TODO|2015-12-22|1| Test if even proportions of mutation types (duplicate, insert, delete) promotes greater numbers of factors. (Can influence be separated from optimisation? - only if number goes down). Could set deletion rate higher than combined duplicate+insert.
TODO|2015-12-22|1| Alter optimise(), createBindingFactor and mutateBindingFactor to use a prefix to name factors. For optimise, this could be the run and iteration name.
TODO|201X-XX-XX|1|
TODO|201X-XX-XX|1|



### Notes (reverse chronological)

__2015-12-22__ 

Feeling a bit negative about this idea. Part overwhelmed, part skeptical. There seem too many factors that cannot be accounted for in this framework. 

- post-transcriptional modification
- extra-cellular regualation
- long-range enhancer effects
- cooperative binding
- stable vs transient binding
- local tethering of factors (e.g. factors may stick around to bind again nearby)
- changing factor abundances with feedback
- regulatory loops (positive and negative feedback)

Some of the above may be approximated by latent variables or careful abstraction.


On the bright side, I have made good progress in implementing quite a few features. 

- Optimisation does improve the factor set to be better than random. 
- Using IRanges, I can run layerBinding over a whole chromosome in a few minutes or less (30 factors)
- Parallelised to optimise faster (e.g. testing 16 mutations on 16 cores)
- Can be run on much smaller memory system (with loss of speed).
- Several mutation types and more to follow.
- Can pass custom scoring functions  (need to tidy this up though).


Thinking about allowing number of factors to mutate. If allow insertion, deletion, duplication with equal probabilities, then there may be a mutational pressure to increase overall numbers (insert + duplicate vs delete). However, my current selection regimen (winner takes all) does not really promote neutral mutations. Could easily test for this (add to TODO table).


__2015-12-21__  Progress on runs

All runs finished over the weekend

- L5_chr19.par (1kb)	original scoring method on 1kb targets centred on TSS.	Exactly 7 days = time limited on cluster
- L5_chr19.400	original scoring method on 400bp targets centred on TSS.	5 days and 22 hours
- L5_c22.acc	Accuracy	1 day and 1 hour (limited to 100 unchanging iters)
- L5_c22.ppv	Positive predictive value	5 hours (limited to 100 unchanging iters)
- L5_c22.tpr	True Positive Rate	12 hours (limited to 100 unchanging iters)

__L5_chr19.par (1kb):__ Has covered 9Mb (~15%) of the chromosome using 181k regions to hit 11273/11451 (98%) of targets. 

	[1] "x16"
	[1] "Round 2996 . Marks on target layer: 181047 , Coverage: 9248198 , Regions with a hit: 30369 , Targets Hit: 11273 , Chrom size: 59128983 , Target count: 11451 , Target coverage: 11451000"
	[1] "Round 2996 . OldScore 0.132442166384536 NewScore 0.131884788501455 Better?"
	[1] "No!"


 
__ L5_chr19.400:__ Has covered 4.8Mb (~8%) of the chromosome using 81k regions to hit 9277/11451 (81%) of targets. 

	[1] "Round 2735 . Marks on target layer: 81051 , Coverage: 4792154 , Regions with a hit: 9770 , Targets Hit: 9277 , Chrom size: 59128983 , Target count: 11451 , Target coverage: 4580400"
	[1] "Round 2735 . OldScore 0.0865191892290033 NewScore 0.0861258953846387 Better?"
	[1] "No!"
	[1] "No improvement in 1001 iterations, exiting!"
	     user    system   elapsed
	105276.58  73345.51 510483.86

 
__ L5_c22.acc:__ Has covered 10kb (~2%) of the chromosome using 301 regions to hit 177/4070 (4.3%) of targets. The high score mainly represents the low false negative rate as most of the chromosome is (correctly) unmarked. > 0.97 was achieved already by iteration 3.

	[mqbssdgb@login2(hydra) seqPredict]$ tail ../q_out/L5_c22.acc.400.o93111
	[1] "Round 489 . Marks on target layer: 301 , Coverage: 10261 , Regions with a hit: 86 , Targets Hit: 177 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 489 . OldScore 0.970765734709443 NewScore 0.970742726486378 Better?"
	[1] "No!"
	[1] "No improvement in 101 iterations, exiting!"
	    user   system  elapsed
	25219.92 10893.50 89939.56

 
__ L5_c22.ppv:__ Has covered <2kb (0.0037%) of the chromosome with 95 regions to hit 106/4070 2.6% of the targets. This seems very precise (and misses almost everything). I suspect it may be hitting one or two clusters of many TSS (but I should check). Subjectively, I would say this is not as good as the accuracy score.

	[mqbssdgb@login2(hydra) seqPredict]$ tail results/pfs_layer5_chr22_400bp_ppv/L5_c22.ppv.400.o93132
	[1] "Round 126 . Marks on target layer: 95 , Coverage: 1903 , Regions with a hit: 54 , Targets Hit: 106 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 126 . OldScore 0.657402728650834 NewScore 0.601681555438781 Better?"
	[1] "No!"
	[1] "No improvement in 101 iterations, exiting!"
	     user    system   elapsed
	 6073.188  2830.469 20344.988


 
__ L5_c22.tpr:__  Has covered 44Mb (86%) of the chromosome using 68k regions to hit 3269/4070 (80%) of targets. Optimised for finding targets doesn't care about false positives.
	
	[mqbssdgb@login2(hydra) seqPredict]$ tail results/pfs_layer5_chr22_400bp_tpr/L5_c22.tpr.400.o93180
	[1] "Round 220 . Marks on target layer: 68459 , Coverage: 44350954 , Regions with a hit: 3269 , Targets Hit: 4070 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 220 . OldScore 0.894632868425867 NewScore 0.88500254965998 Better?"
	[1] "No!"
	[1] "No improvement in 101 iterations, exiting!"
	    user   system  elapsed
	13441.01  5836.35 45145.92


The performance of the different scoring metrics is dependent on the relative numbers of TP and TN. When searching for small regions in a large chromosome, TP/TN is very low. 

It is also not obvious whether the 'standard' scoring metrics are any better than the custom method I came up with. What is more important, starting transcription around about the right place, or preventing transcription in the wrong place?


__TODO__ Functions to inspect output factorSets from above runs. Are they targeting randomd DNA?  Would a larger number of factors get better scores? More layers?
__TODO__ Implement factor abundance on a scale from 0 to max. Then mods can be applied by the proportion of each factor. Factor abundances can mutate.

__2015-12-18__ Using accuracy (acc) as scoring method

This was really interesting. I fudged the code a bit but managed to get optimiseFactorSet() to work with acccuracy (TP + TN) / (TP +FP +FN +TN) as a score to optimise. Accuracy is basically the proportion of all sequence in the correct class.

It began with what looked like good scores for the targets. 

	[1] "Round 1 . Marks on target layer: 77537 , Coverage: 8026071 , Regions with a hit: 2068 , Targets Hit: 2553 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 1 . OldScore 0.554783338607464 NewScore 0.768401309130101 Better?"
	[1] "Yes!"
	[1] "x16"
	[1] "Round 2 . Marks on target layer: 57081 , Coverage: 3569626 , Regions with a hit: 1444 , Targets Hit: 1762 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 2 . OldScore 0.768401309130101 NewScore 0.876863460548262 Better?"
	[1] "Yes!"

but very quickly 'decided' to ignore the promoters and optimised to mimismise the False positives, which is easy if it just doesn't mark very much.

	[1] "Round 3 . Marks on target layer: 186 , Coverage: 12566 , Regions with a hit: 4 , Targets Hit: 5 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 3 . OldScore 0.876863460548262 NewScore 0.970420397660052 Better?"
	[1] "Yes!"
	[1] "x16"
	[1] "Round 4 . Marks on target layer: 180 , Coverage: 11828 , Regions with a hit: 4 , Targets Hit: 5 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
	[1] "Round 4 . OldScore 0.970420397660052 NewScore 0.970437695540015 Better?"
	[1] "Yes!"
	[1] "x16"

Note the sudden drop from 57081 regions (3.5Mb) to 186 regions (12.5kb). It reached accuracy of 0.97 in 4 rounds (x16 cores) and I stopped the job.

So, 'accuracy' is very good at not marking bits of the genome containing promoters when the promoters are a small fraction of the genome. However, for defining where such promoters are 'accurately', it is not much use.

__ADDENDUM__ I subsequently noted that the number of TSS hit had begun to increase again by the 8th cycle. So I renamed the run (pfs_layer5_chr22_400bp_acc) and restarted it with logging every iteration.

__DONE__  Start directing the STDOUT and STDERR error logs to the output directory for each run. They are too informative to lose.
__TODO__  I think I've considered this before, but will need some way to (visually?) summarise the resulting factorsets. Do they have an excess of a type of factor (currently limited by me)? Do they feature certain types of binding/modifying early or late in the binding order? Are they simply finding GC rich (or AT poor) regions?  Do some have more influence on the output than others?

Also set of a similar job (pfs_layer5_chr22_400bp_ppv) to use positive predictive value (precision, TP/(TP+FP)) as optimisation score. 

Both ppv and acc runs seem to be hovering around a small number of true-positive by only marking a small fraction of the genome. I'm slightly worried that they are very precisely matching a small set of (possibly overlapping TSS). Over-fitting?

I am logging every round though, so in theory could visualise over-fitting by re-run each cycle's factorSet against a different chromosome. Might expect to see the performance increase as general (true) characteristics are picked up and then drop as the factorSet over-fits to the training chromosome.


__Future__ Thinking about moving to finding selected TSS and perhaps not finding unselected TSS (1,0,-1 ? ). I was going to attempt to train towards either broadly expressed promoters (house-keepers) or the transcriptional state of a pluripotent stem cell (or ESC). However, I'm not sure how many TSS can be found reliably (current methods seem to be finding a small proportion accurately). 

An alternative method would be to work toward tissue-specific TSS sets. I could specifiy partially overlapping sets that are enriched in either(both) of two tissue types. Then optimise for their marking separately and finally combine or contrast the optimsied sets. Liver and Brain might be two distinct tissue types. OR I could find a dataset of known targets of a certain TF and try to learn the TFs binding factor. The problem with both of these is tissue specificity and I don't expect pfs to cope yet with regulation from enhancers outside the immediate promoter. 

I could use a dataset that relied on promoter-proximal binding: TATA, REST?

[Paper featuring TFs that bind close to the TSS](http://nar.oxfordjournals.org/content/36/21/6795.short). It is possible that different subsets of these promoters have proximal and distal binding sites that may operate at different times/conditions. HNF4A and liver are candidates.

the above was cited by another paper using primary [DNA sequence to predict promoter activity!](http://nar.oxfordjournals.org/content/39/11/e75.full)

[Same team 2011](https://www.jstage.jst.go.jp/article/gi/25/1/25_1_53/_article)
[Same team 2011b](http://genome.cshlp.org/content/21/5/775.full.pdf+html)

Do they model changing binding specificity of factors through time?  




__2015-12-17__ 
Killed the chr7 parallel job started on 10/12/2015. Only managed 500 rounds in a week. Are large chromosome going to be too large to handle. Could they mutate into slow runs? Maybe build in an execution time limit for each modification set.


__2015-12-16__ 

Yesterday I implemented some standard scoring functions (accuracry, tpr, fpr etc). I wanted to use these as optimisation scores. I added them to a script (pfs.scoreMeasures.R) and included their use in a run script (pfs.chr22.400bp.par.R). It failed though because the functions were not visible on the parallel nodes. The optimisation function does the export. 

__TODO__ Either re-write the optimiisation function to export all functions in the environment to cluster nodes OR package the code properly. OR add a customisable watch_function to the optimisation function.

__2015-12-14__ Progress on latest run (chr19, more detailed output)

	 "Round 668 . Marks on target layer: 76205 , Coverage: 7133959 , Regions with a hit: 14718 , Targets Hit: 10565 , Chrom size: 59128983 , Target count: 11451 , Target coverage: 11451000"
	[1] "Round 668 . OldScore 0.125974998293263 NewScore 0.125704617482101 Better?"

So, the current scoring is somewhat misleading. I think the max score for this chrom may be < 13%. 
Layer 5 is marked over 7Mb (of 51Mb chrom (13%) and 11.4 Mb of target (expanded TSSs)). 
There are 76k regions and ~15k of them overlap a TSS region.
Of 11451 TSS regions, 10565 (92%) are hit by the regions.

This is a multi-optimisation problem: the # of TSS hit  versus the proportion of non-TSS chromosome marked.

- __TODO__ alter scoring and/optimisation
- __TODO__ produce graphics of chr19 layer5 vs TSS.IR to show fit.
- __TODO__ run with reduced region as target.
- __TODO__ test if chr19 optimised factor set is any good at finding TFs on other chroms.
 
Explored the currentFactorSet of the unfinished chr19 run (see load.pfs.chr19.1kb.par.R). Ran runLayerBinding with a custom 'watch.function()' to print out the number of TSS hit as each factor is applied.

	> modLayerSet2 <- runLayerBinding(layerList=layerList.1, currentFactorSet, iterations=modsPerCycle,  collect.stats=TRUE, target.layer=target.layer, verbose=TRUE, 
	+                                watch.function=function(x, target.vec=tss.IR , ...) {print(paste("Watch function", length(x[[target.layer]]), sum(overlapsAny(target.vec,x[[target.layer]]))))})
	[1] "2015-12-14 17:37:34 runLayerBinding pos 1"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.25 673756 33334[1] "Watch function 0 0"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.7 19434 19434[1] "Watch function 0 0"
	[1] "Hits match whole chromosome"
	bf.16 1 33334[1] "Watch function 0 0"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.4 43332 33334[1] "Watch function 33334 8445"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.12 396 396[1] "Watch function 33568 8478"
	[1] "Hits match whole chromosome"
	bf.27 1 33334[1] "Watch function 33568 8478"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.17 4715 4715[1] "Watch function 33568 8478"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.3 42890 33334[1] "Watch function 64087 10225"
	[1] "Hits match whole chromosome"
	bf.24 1 33334[1] "Watch function 64087 10225"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.11 23544 23544[1] "Watch function 68740 10375"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.6 28080 28080[1] "Watch function 69879 10428"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.14 55643 33334[1] "Watch function 69879 10428"
	[1] "Hits match whole chromosome"
	bf.5 1 33334[1] "Watch function 69879 10428"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.1 2727 2727[1] "Watch function 69879 10428"
	[1] "Hits match whole chromosome"
	bf.22 1 33334[1] "Watch function 69879 10428"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.20 171 171[1] "Watch function 69879 10428"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.29 721 721[1] "Watch function 69879 10428"
	[1] "Hits match whole chromosome"
	bf.26 1 33334[1] "Watch function 69879 10428"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.19 25472 25472[1] "Watch function 71224 10428"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.8 19523 19523[1] "Watch function 70708 10436"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.9 21443 21443[1] "Watch function 73308 10514"
	[1] "Sequence of length  59128983 , using  6 windows of length 1e+07"
	bf.23 33685 33334[1] "Watch function 75920 10568"
	[1] "2015-12-14 17:42:00 runLayerBinding.fast pos 2"
 
BF.4 is a poly-B run (not A) and BF.3 is a poly-V run (not T). other jumps include other similar factors. Commonly GC-rich regions, not sure if not-A or not-T are more discerning than all-GC.

__TODO__ move on to testing against specific TSS being activated (quantitative?) or de-activated.
__TODO__ add ... to optimise function to allow watch function on layer binding (probably won't work with parallel jobs).

 
__2015-12-11__  Progress on latest run

L5_chr19.p mqbssdgb     r     12/12/2015 03:23:05	17 hrs
 "Round 213 . Marks on target layer: 120589 , Coverage: 8123230"
[1] "Round 213 . OldScore 0.127049312140155 NewScore 0.126935690792474 Better?"
[1] "No!"

# Started with a score of 6%, quickly attained 12.7%.  Need more output on number of TSS with ANY overlap and # of regions with a TSS target overlap.


__2015-12-11__: Progress on long runs 08:50

	L5_chr7.par.o91126	Round 86	0.0645244973393824	12/10/2015 11:17:31	~22hrs
	L5_chr22.par	Round 616	0.0885717943907263	~47hrs

What is being optimised?
On chr22, there are ~4500 transcripts and the chrom is ~51Mb long. When extending the targets to 1kb, the total target area is 4.5Mb. 4.5/51 = 0.08823529
Certainly need to collect more stats on coverage of target.layer during optimisation.

Current test function:-

	test_function <- function(layerList, targetLayer=target.layer, target.vec)  {
	  inter.size <- sum(width(intersect(layerList$layerSet[[targetLayer]], target.vec)))
	  union.size <- sum(width(union(layerList$layerSet[[targetLayer]], target.vec)))
	  return(inter.size/ union.size)
	}


I re-load the current best factorSet from results/pfs_layer5_chr22_1kb_par/currentFactorSet.600.Rdata (manually but commands stored in scripts/pfs.reloadResultOnHydra.R).
I ran runLayerBinding a couple of times. Each time, Layer 5 (target) had only (and exactly) 716 bp coverage in 10 regions (c.f. 4070000 for target). 

I don't think it is marking as much of the chromosome as I expected and there may be a clue in output from runLayerBinding():-

	> modLayerSet2 <- runLayerBinding(layerList=layerList.1, factorSet = currentFactorSet, verbose=TRUE)
	[1] "2015-12-11 10:32:43 runLayerBinding pos 1"
	[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
	bf.1 25885 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
	bf.2 6878 1bf.3 1 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
	bf.4 17615 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
	bf.7 37536 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
	bf.8 287 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"

I don't think it is marking as much of the chromosome as I expected and there may be a clue in output from runLayerBinding():-

The tss may not be unique 
	
	> sum(width(intersect(tss.IR, tss.IR)))
	[1] 2287192
	> sum(width(tss.IR))
	[1] 4070000

NO, checked that and they are made unique. The problem is that using 1kb, they overlap each other

	> sum(width(tss.IR))
	[1] 4070000
	> sum(width(reduce(tss.IR)))
	[1] 2287192


The layer_region types are only marking a single central hit. Because:-

    stateWidth <- bindingFactor$mods[[thisLayer]]$stateWidth
    hits <- resize(hits, width=stateWidth, fix="center")    # adjust the width

For a layer_region matching 'N' and no other layers, matchBindingFactor returns an IRanges spanning whole chrom. Resize converts this to a small span (stateWidth) in the centre of the chromosome.

Should such bindingfactors be allowed?  If so, need to preserve number of hits (e.g. 1:N for length(chrom)). May add to memory/time requirements (probably why I cut it down to a single region in the first place).

I added an edge-case to runLayerBinding such that if the hits object is an IRanges object spanning the whole chromosome, then it will generate max.hits hits randomly across the chromosome. 

__2015-12-10__: Progress on long runs. 11am

	Lay5_chr22.1kb	Round 2122	0.0885956837508925	12/07/2015 22:33:35	~ 61 hrs, looks like no improvement for 900 iters, expect to terminate today.
	L5_chr22.par	Round 370	0.0874762103883457	12/09/2015 10:03:18	~ 25 hrs

I wonder if this is limited simply by the proportion of the target that can be covered. If so, increasing the number of marks that occur duing layerBinding may improve the fit.

__2015-12-09__:  I left the chrM job running with 10,000 iterations. The computer rebooted itself to install updates but it looks like the script properly exited after 1000 iterations without improvement in about 2hours (it reached 33% (~12 of 36 TSS) by 2000 iterations ).

I stopped the original pfs run on chr22 after 34hrs/1600 iters. It was stalling on 4% for 500iters and I wanted to get a parallel job started on the server. chr22:  47 iters/hr. Be careful when comparing this with parallel jobs because the scripts currently generate and test 16 (n.cores) mutated sets PER ITERATION. The parallel chr22 completed 6 iters in 30mins but this represents testing 6*16=96 factorSets.

__DONE__ remove verbosity from factor mutation, it is all I can see in the logs.
__DONE__ with parallel jobs, log scores more often (because many more factors are tested).
__TODO__ better reporting on optimisation e.g. score, size of seq, width of target regions, width of target seq modified, number of regions.
__TODO__ more liberal mutations: duplications (limited?), switching, quantity/proportion of each factor

__2015-12-08__: Last night began the first long optimisations (10,000 iters) using pfs code on whole of chr22. The first set-off @ >500 iters by morning (approx 40/hr). Will take 20 days (so will probably be killed).

[DONE 2015-12-08 ] more speed up (parallelise optimisation: mutate each accepted set 10 times and run each independently, keep the best.

The two runs were for 200bp centred on TSS and 1kb centred on TSS. The 200bp is slow to improve (<4% after 500i), the 1kb started slow, scoring 0.0 for most of the early runs but jumped to 7% within the first 100i and then is still <9% after >300i. One positive is that the scores seem quite robust, most failed runs are very close to the current optimum. Similar code run on chrM ran 100 iterations in 3 minutes and achieved 18% (of only 36 TSS).

__N.B.__ I added the strict ordering to make computation tractable. I wonder if ordered regulation in biology is adaptive because of robustness (?). 

Is there a theory of gene regulation that suggests how wide a promoter would be good to aim for. 

The current system has many flaws:-
- no long range in-cis effect (__TODO__ could be added in by enabling offset parameter in bindingFactor)
- improper control over the numbers of each bindingFactor (currently an equal number of each (or fewer).
- It is very slow to get started because factors binding higher layers won't in the initial state. (unless they bind to empty state).



__2015-12-07__: The new windowing approach was successful on chr22 (DNA-motif in 25 seconds with 5Mb window, 17 secs with 10Mb window) but revealed 'matches' across the first 16Mb of the chromosome (the telomere is NNNs).  DESIGN DECISION: whether to allow matches to Ns?  I think not, because the telomeres and centromeres are a relatively small proportion of the genome but could suck up a large number of hits, when there is no evidence there should be hits there.

[DONE 2015-12-07] filter out Ns from telomeres etc. (or don't match to them).   Solved using fixed='subject'

[DONE 2015-12-07] Window size against a single layer layerSet seemed ok up to 10Mb. Need to test on 5 layer system.

Ran 1 layer chr22 with 30 binding factors in 3m 20secs. No strain on memory on my PC.

Also attempted short optimisation, noticed that many factors giving exactly one hit:-

	[1] "2015-12-07 17:02:48 runLayerBinding.fast thisBF = bf.10 N"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast n.hits = 1"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast n.hits.used = 1"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast thisBF = bf.11 N"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast n.hits = 1"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast n.hits.used = 1"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast thisBF = bf.12 N"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast n.hits = 1"
	[1] "2015-12-07 17:02:48 runLayerBinding.fast n.hits.used = 1"


__2015-12-04__: I spent quite a few hours re-writing functions to use IRanges() objects instead of vectors of 0/1 to mark the layers (see pfs.functions.R). I developed using chrM but as soon as I tried an optimisation with chr22, I hit the same memory problems as before. One thing to appease it was to return a single IRanges spanning the whole sequence for "layer_region" and "layer_island" matches. But the main problem is that DNA matches take a lot of memory. Maybe, split the sequence into segments of a given size (e.g. 100kb), create the IRanges for each and then concatenate them (parallel jobs?). Make sure the splits overlap by more than patternLength. Then re-combine the IRanges with reduce().

I did that, I added max.window parameter to matchBindingFactor. Left it running on chr22. Need to use as large a window as possible (e.g. 1-10Mb).



__2015-12-02__: The chr7 runs were failing on a Hydra full node (512Gb) without an error message. I scaled down to chr22 and got this error message:-

	[1] "2015-12-01 21:54:21\trunLayerBinding.fast thisBF =\tbf.9"
	Timing stopped at: 429.981 33.544 463.486
	[1] "Error in .Call2(fun, object@ptr, ..., PACKAGE = \"IRanges\") : \n  cannot allocate memory block of size 67108864 Tb\n"
	attr(,"class")
	[1] "try-error"

Hmmm. Something is very wrong if trying to use 67 million terabytes of memory!

Speed up/ memory down options:
	matchBindingFactor() I think could be faster. It might also be possible to return a list of potential modifications rather than a list of hits. This would then remove the need work out what the mods are. It depends on whether it's faster to (a) work out all the hits but only calculate a subset of mods, or whether calculate all the hits and their mods in one go and apply a sample (probably still the former?).

I added in Rprof() to try and profile the code, but it didn't give any output (even in logging file).

With extra verbosity, I was able to see that the script failed when trying to matchBindingFactor() for factors with long strings of degenerate bases (e.g. TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT). This may be a case where BSgenome works better than loading a chromosome as a DNAstring set? I did a test and found that searching whole chromosomes was extremely slow for this, mostly because of the emptyLayer matching. Decided to re-code into just using HITS objects and maybe remove island types? (though IRanges are rle objects..). One problem is that I wanted to allow mismatches. 

__2015-12-01__: Just trying to do something with this, not sure what. Set off a training/optimsiation run against the whole of chromosome7 ([script](../scripts/predictFromSequence.chr7.R)). This is approx. x1000 the size of previous regions so I upped the number of mods from 1,000 to 1,000,000. The first attempt on Hydra failed after about 4 hours with a core dump. That was using 8 cores (~256Gb RAM). It also only managed about 3 iterations in that time. This morning I re-started with 16 cores. May need to re-write the modifyLayerBinding again to make it __MUCH__ faster

Should probabably accept mutations that are no worse than a given % (rather than only accepting muts that are better).

__2015-09-14__: Need a set of targets to aim for. Can demonstrate learning of TSS. Might be good to have set of broadly expressed TSS and optimise first for these. Then bring in tissue-specific or developmental TSS and add more factors to explain tissue-specific expression. Alternatively, could train against all of them, expecting it to perform better for broad TSS. Show that it trains to TSS but also that it is better at predicting the broad promoters. Perhaps get lists from other papers on TSS classes. OR just use expression data across a set of tissues. CAGE data may be best as based on speicifc TSS. 

__2015-09-11__: Reading Irie2011Predicting (promoter prediction using all TFBS).  There is lots of data out there on TSS activity profiles (transcriptomes). They claim high predictive power (>80% just using TFBS sites).  For this project, would need to show that cell-lineage memory (through layer marking) can improve on simple knowledge of TFBS.

__2015-09-08__: Reading a review on promoters (Lenhard et al., NRG, 2012) detailing different clases of mammalian promoters (or perception thereof) e.g. tissue-specific, ubiquitous, developmentally regulated. It might be 'novel' to skip these distinctions and just generate transcriptomes from sequence. First for ES cell (or progenitor gut cell?) then with that as baseline (marking and/or factors) develop coding for distinct tissues. Hope to see silencing of many tissue-specific genes in non-specific tissues.

What about ribosomal promoters and other 'high-performance' promoters (e.g. translation initiation). The above review suggests they may have distinct pattern too. If I use RNA-seq data to create model of expression, then these genes may have been depleted prior to sequencing and have artificially deflated counts - impacts on model testing. Same goes for mitochondria?

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

Should probably also look at prior art

- [list of softwares](http://molbiol-tools.ca/Promoters.htm)
- [a larger list of softwares](http://www.genetools.us/genomics/Promoter%20databases%20and%20prediction%20tools.htm)
- [An even larger list at geneprediction.org](http://www.geneprediction.org/software.html)
- [Fruitfly](http://www.fruitfly.org/seq_tools/promoter.html)

>Training set:
Our training and test sets of human and Drosophila melanogaster promoter sequences are available to the community for testing transcription start site predictors. These sites also contain our representative, standardized data sets of human and Drosophila melanogaster genes.

- [Ohler lab](https://ohlerlab.mdc-berlin.de/software/)


Looks like I'm 15 years too late!  However, most of these focus on small patterns over a few 100 bps. They (mostly) don't model a populatin of large DNA binding factors over the 100kb scale. I suspect they don't model availability of binding sites (genome- or region- wide) to influence lieklihood of bindng.

Things that may still be novel:-

- any part of the genome can be a tss! Need to predict not only where active tss are, but also where they are not.
- model regional and/or genome-wide competition for binding sites
- model 'memory' of system provided by marks placed on the sequence (histone marks, dna-methylation, nucleosome presence, compaction, competitive binding, etc).
- predict __when__ promoters are active.


__2015-09-03__: Ran first script on the cluster for 10,000 iterations. No better than the first 1000 run a day earlier. May have been a bit of luck involved there. Still progressed up from 4% (iter 20) to ~9.5%  (although most achieved in first 1000 iters). Mostly just a requirement test for such a script. Parsing the resulting factors is still not very easy, may need better summarisation scripts. Also noted that many of the factors were to bind anywhere and set the target layer to 0 - not a bad strategy when the target size is a relatively small proportion of the full sequence. 

Also started a 5-layer test today (10,000 iters). Despite expectation that this might slow the progression, early inications are that it is already scoring high!

	[1] "Round 187 oldScore 0.135323991505638 newScore 0.173346165894045 Better?"
	[1] "Yes!"
	
Some of the non-improved scores are still very low, which might indicate the system is easily broken but more likely is due to the fact that scoring is still stochastic and is picking rounds that largely by chance were high. __TODO__ remove randomness in layer modificiation (but how? saturation of all possible sites is only way I can think and would be VV. slow). OR mutate at lower rate and run modification several times to get average score. This second option would make the scores more robust but would slow the system too. The second option also would avoid saturation, which is good as eventually I intend to introduce factor abundance as a parameter - i.e. I intend that not all sites are occupied.

Looking at the factor set of the five layer case. It seems it is going to be very difficult to read what is happening in these traces. It might be a good idea to make a 'movie' version of runLayerBinding to watch the sequence of changes that occur.  With 5 layers, and 30 factors, could be tricky to visualise.

If this is 'working', need to start doing proper training and testing on different sets of sequences.

- sort out stochastisity of tests
- add in multi-sequence training and testing sets, with cross-validation.



__2015-09-02__:  [predictFromSequence.devFunc.R](scripts/predictFromSequence.devFunc.R) Wrote alternative to runLayerBinding() called runLayerBinding.fast() that does not apply mods in random order and re-calculate, it just does them all for a single factor at a time. This breaks the idea of sequential mods creating pattern matches for other factors. The point was to increase speed enough to see if the system could find any TSS reliably. I also widened the goalposts by extending the TSS to a promoter with 200bp upstream and 200bp downstream. I created a random set of factors and a 1-layer LayerSet (layerList) and ran optimiseFactorSet() [n.iter=1000, mut.rate=0.1, modsPerCycle=10000].  This quickly scored above 1% and then to 3% (c.f. ~ 0.1%, previously), which is not too surprising with the widening of the targets.  The best score achieved was 12% and the plot of scores (by iteration index) showed a step-change around iteration 700.  Due to the stochastic nature of runLayerBinding.fast(), the set of factors does not produce the same result each time (range 3% - 12%). I therefore re-ran runLayerBinding.fast() 100 times without modifying the factor set, counting for each base the number of times it was marked. This gave a very nice set of peaks and troughs, some of the peaks being directly below TSS/promoters. The cor score for the summed vector and the tss vector was 12%. Using a threshold to divide the summed scores into bound/unbound did not help the cor score, possibly because the zeros are informative - in the test vector (HOXA cluster), one obvious pattern is the absense of prediction in the large gene desert to the left (chromosomal view) to the left of the HOXA cluster.

LayerSet became encapsulated within LayerList (to store extra annotations) and FactorSet needs to be encapsulated to carry the optimScores table (and in future other data) returned by optimiseFactorSet() __TODO__. 

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

### __TODO__ list.
Specify factorSet to include optimscores (e.g. make a list like layerList)

Spilt off calculation of binding positions from binding mods to speed up a series of mods.  Many TFs bind at independent locations.  Only a problem for binding events that depend on previous events but genome is so large and number of event so large, that many could be allowed to happen without re-calculation. e.g. if running 100k mods, could recalculate after every 1000. Practically, could split off binding calcs into separate thread. ASSUMPTION: # of events far greater than number of types of factors. [speed-up|design]

__DONE__ Additionally, create alternative algorithm ('FAST') that changes a high proportion of possible sites (-> 1.0) with every iteration. BUT many compositional factors will match all along the sequence or at many thousands of overlapping locations. These would then create many thousands of redundant hits to be processed, and also may swamp other signals. 

If optim starts to work. Track when factors get added to factor set. Also need way of measuring importance of each factor.

__DONE__ In script logs, log initial score.

Multi-sequence (random) training and testing sets. Concatenate seqs? - no. How to avoid over-fitting? Need fast way to get a bunch of sequence. May be worth using BSGenome and write functions to retrieve arbitrary segments. Would also then need to pull all relevant data for each region.

Write parallel version of optimiseFactorSet to use scores from a range of sequences (take mean or median score).

Are the factorSets strand specific? Could program to make searches on both strands.

Memory. The multi-layer objects take up lots of memory and can be a bit slow to work with.  Is there a way to use 'views'. Also, instead of retaining very long binary vectors, is it possible to use GRanges and set operations to flip ranges between 1 and 0? I think that might be a lot quicker than altering ranges of a BString. 

How to 'train' the system?  
	Divide into bins? Training Set, Testing Set and Validation set?
What about tissue/stage specific expression? Also want to 




