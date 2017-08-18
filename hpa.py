#!/usr/bin/env python

"""
Standalone version of HPA algorithm.

Some code from hmatrix3.py was duplicated here to make this self contained and more accessible,
and also modified (for example this code uses 0..ntaxa-1 internal taxa indices instead of 1..ntaxa.

"""

__project__ = "HPA"
__version__ = "1.0"
__author__ = "Guy Hoelzer and Rich Drewes"
__copyright__ = "Copyright 2017 Guy Hoelzer and Rich Drewes."
#__citation__ = ""

import time
import random
import sys
from collections import Counter, defaultdict
import dendropy
from itertools import combinations
import operator
import numpy
from numpy import where
from math import sqrt

def FastCompatibleJustSynsets(c1, c2, verboselevel=0):
    """ Check if lists of partitions (or single partition) in c1 and c2 are compatible with each other (each partition in c1 doesn't overlap another partition in c2, or is a subset of another partition in c2, and vice-versa).
    """
    # c1 and c2 have been called "synsets" in the past
    # note: the c1 and c2 are each assumed to be self compatible and there are no duplicate taxa listed within c1 or within c2
    if verboselevel > 1:  print "looking for compatibility between:", c1, "AND:", c2
    compatible=True
    if type(c1) is not type([]): c1=[c1]
    if type(c2) is not type([]): c2=[c2]
    # compare each partition in the first list of partitions (c1) to each partition in the second list of partitions (c2)
    for p1 in c1:
        for p2 in c2:
            if verboselevel > 1: print "FC:", p1, p2
            if len(set(p1) & set(p2)):  # share an element, and . . .
                if verboselevel > 1: print "the partitions share an element, so one must be a subset of the other for campatibility"
                if len([p1t for p1t in p1 if p1t not in p2]) == 0:
                    if verboselevel > 1: print "and the first is a subset of the second, so OK so far"
                    continue
                if len([p2t for p2t in p2 if p2t not in p1]) == 0:
                    if verboselevel > 1: print "and the second is a subset of the first, so OK so far"
                    continue
                if verboselevel > 1: print "neither is a subset of the other, so these partitionings are incompatible"
                compatible=False
                break
        if not compatible: break
    if verboselevel > 1: print "compatible:", compatible
    return compatible

def CoOccur(lt, ntaxa, allfreqsd, truetree, restrictsize=None, verboselevel=0, taxanames={}, includesingle=False):
    """
        See how many times taxon OR partition lt cooccurs with each other taxon, among the partitions in the global partition frequency occurrence table allfreqsd.

        allfreqsd is dict of frequency of occurrence of partitions given by keys.
        As always, the partitions are represented as tuples of increasing taxa numbers e.g. (1, 4, 9) but not (1, 9, 4).

        Returns list of (tn, n), sorted by decreasing cooccurrence count n, where tn is set of target taxon (those taxa not in lt) listed in increasing order.

        if restrictsize is not None, only look at partitions of size restrictsize

        lt can be a tuple or int taxon number.

        Looking at simple cooccurrence should isolate the effects of long branches which can otherwise cloud
        partition frequency counts by e.g.: long branch 11 creating partition (11, 18, 19, 20) with some
        frequency f1 in addition to a real partition (18, 19, 20) with frequency f2, when 18, 19, 20 really
        cooccur f1+f2 times

        Note: internal taxa indices are 0..ntaxa-1
    """
    intt=" "
    if type(lt)==type((1,)):
        if len(lt)==1 or lt in truetree:
            intt="*"
    if restrictsize==None:
        if verboselevel > 1:
            print "how many times does taxon/partition", intt+`lt`,"cooccur with other taxon (in all partitions)?"
    else:
        if verboselevel > 1:
            print "how many times does taxon/partition", intt+`lt`,"cooccur with other taxon (in all partitions of size %d)?" %restrictsize

    if type(lt)==type((1,)):
        # lt is a partition stored as a sorted tuple, not just a single taxon
        ot=[tn for tn in range(0, ntaxa) if tn not in lt]
        cocc={}
        for tn in ot:
            cocc[tn]=0
        maxi=0
        maxn=0
        # count cooccurrences of the partition lt in partitions with each other taxon tn
        for p in allfreqsd.keys():
            if restrictsize is not None and len(p) != restrictsize: continue
            for tn in ot:
                if set(lt) < set(p) and tn in p:
                    #print "cocc", p, tn, allfreqsd[p]
                    cocc[tn]+=allfreqsd[p]
        otl=[(tn, cocc[tn]) for tn in ot]
        if len(otl)==0: return None
        otl.sort(key=operator.itemgetter(1), reverse=True)
        tn, lastn=otl[0]
        dropthresh=2            # this percent drop is deemed important. USED ONLY FOR PRINTING!
        if verboselevel > 1:
            for tn, n in otl:
                #if tn==maxi: ismax="*"
                #else: ismax=" "
                marker=None
                if n>0:
                    #droppct=float(n)*100.0/float(lastn)
                    droppct=(float(lastn)-float(n))*100.0/float(lastn)
                    if droppct>dropthresh: marker="-------"
                    droppct='%2.1f' % droppct
                else:
                    droppct=" "
                if marker is not None: print marker
                print lt, '\t', tn, '\t', n, '\t', droppct
                lastn=n
        return otl

    # lt is just a single taxon, not a partition. this is the normal situation.
    ot=[tn for tn in range(0, ntaxa) if tn!=lt]     # ot means 'other taxa' besides the one we are comparing against (lt)
    cocc={}
    maxi=0
    maxn=0

    # for good sized data set (40 taxa, 200K or 500K bases):
    # 11m02s python on original version
    #  6m51s python on second version
    #  2m25s pypy on second version
    #  1m46s pypy on second version with no long branch check, so CoOccur called once instead of twice
    # (effect on memory usage unknown, but probably not worse and maybe even better)
    # (may actually slow things down to use pypy on small ntaxa or small seqlen, but makes up for it on large data)

    # original version, before optimization attempt:
    if 0:
        for tn in ot:
            cocc[tn]=0
        for p in allfreqsd.keys():
            if restrictsize is not None and len(p) != restrictsize: continue
            for tn in ot:
                if lt in p and tn in p:
                    #print "cocc", p, tn, allfreqsd[p]
                    cocc[tn]+=allfreqsd[p]

    # optimization attempt 1, assuming dictionary lookups are expensive, only do that once:
    if 1:
        if restrictsize is not None:
            print "CoOccur: restrictsize no longer supported in optimized version of CoOccur! exiting"
            exit()
        for tn in ot:
            cocc[tn]=0
        #lll=allfreqsd.items()
        lll=[(p, allfreqsd[p]) for p in allfreqsd.keys() if lt in p]
        for tn in ot:
            cocc[tn]+=sum([v for p, v in lll if tn in p])
        # NOTE! the "cooccurrences of taxon with itself" will always equal the total number of variable positions,
        # since any single taxon will always occur in at least one partition in each VP!
        if includesingle:
            cocc[lt]=sum([v for p, v in lll])
            #print "cooccurrences of taxa with itself:", cocc[lt]
            #print [v for p, v in lll]

    otl=[(tn, cocc[tn]) for tn in ot]
    if includesingle: otl.append((lt, cocc[lt]))
    if len(otl)==0: return None
    otl.sort(key=operator.itemgetter(1), reverse=True)
    tn, lastn=otl[0]
    dropthresh=2            # this percent drop is deemed important. THRESHOLD USED ONLY FOR PRINTING!
    if verboselevel > 1:
        for tn, n in otl:
            #if tn==maxi: ismax="*"
            #else: ismax=" "
            marker=None
            if n>0:
                #droppct=float(n)*100.0/float(lastn)
                droppct=(float(lastn)-float(n))*100.0/float(lastn)
                if droppct>dropthresh: marker="-------"
                droppct='%2.1f' % droppct
            else:
                droppct=" "
            if marker is not None: print marker
            if lt in taxanames and tn in taxanames:
                print "%d:%d" % (lt, tn), taxanames[lt], '\t', taxanames[tn], '\t', n, '\t', droppct
            else:
                print lt, '\t', tn, '\t', n, '\t', droppct
            lastn=n
    return otl

def BuildCorrelationDropList(ntaxa, freqsd, truetree, verboselevel, dropthresh=0.0, taxanames={}):
    """
        Build a correlation drop list for each taxon considered as a source, to each other taxon, among the set of {all partitions} given in freqsd or {all 'significant' partitions} selected from freqsd.

        With a 'dropthresh' of 0.0, we are doing an ntaxa^2 number of correlation computations.
        That could become time consuming for lots of taxa, obviously.

        Internal taxa indices from 0..ntaxa-1 as always.

        Returns tuple (overallbigdropsavgl, overallbigdrops, allcdls):

            overallbigdropsavgl, a sorted list of all big drops averaged over all source partitions
                looks like:
                    [((4, 10, 17), 46.20803989751318, [45.90303373972214, 46.12985540119082, 46.591230551626595]),
                     ((1, 11, 19), 41.76813293957761, [41.68967421314191, 41.765845557708275, 41.84887904788265]), ...
                    which says:
                        partition (4, 10, 17) has average drop of 46.208% to next most correlated taxon,
                            made up of 45.9, 46.12, 46.59 (one numeric value for each source taxon within the partition)
                        partition (1, 11, 19) has average drop of 41.768 made up of 41.76, 41.76, 41.84
                        ...
            overallbigdrops, a sorted list of all big drops with each source taxon listed separately
                looks like:
                    [((3, 14, 16, 18), 55.62438865340724, 13.90609716335181),
                     ((3, 14, 16, 18), 55.40103492884864, 13.85025873221216), ...
                    the second numeric value is an experimental partition length adjustment to the first value
                    just the first value divided by the partition length

            allcdls, a dict with keys for each source taxon giving cooccurrences of all taxon in full set of taxa with 
                each other taxon, among the partitions in the partition count dict freqsd.

        NEW: now we include the self-correlation of partition to itself as an entry in big drops list.
    """
    if dropthresh > 0.0:
        print "Only dropthresh=0 is now supported in BuildCorrelationDropList"
        exit(0)

    allbigdrops={}              # drops separated by taxon, allbigdrops[i] holds the big drops for taxon i
    overallbigdropsl=[]         # all big drops in one big combined list
    overallbigdropsavg={}       # with data of all source taxa for a given partition averaged together for that partition

    allcdls={}

    if verboselevel > 1: print "building correlation drop percentage list"
    includesingle=True
    for i in range(0, ntaxa):
        if verboselevel > 1:
            if i in taxanames: realname=" (taxon name %s)" % taxanames[i]
            else: realname=""
            print "on taxon index %i %s" % (i, realname)
        allbigdrops[i]=[]
        #print "looking at taxon cooccurrences to build a table of big drops", i
        # NOTE: with "includesingle=True" we now include "self correlation" in the drops list, which is 
        otl=CoOccur(i, ntaxa, freqsd, truetree, verboselevel=verboselevel, taxanames=taxanames, includesingle=includesingle)
        allcdls[i]=otl
        # returns sorted list of (taxon, number of correlations to that taxon)
        tn, lastn=otl[0]
        ti=1
        while ti < len(otl):
            tn, n=otl[ti]
            if n>0:
                droppct=(float(lastn)-float(n))*100.0/float(lastn)
                if droppct > dropthresh:
                    # NOTE: we are including the key taxon in the tuple here. elsewhere we don't
                    targetp=[otn for (otn, tn) in otl[0:ti]]
                    if not includesingle:
                        targetp.append(i)      # no longer do this since we called CoOccur to build drop list with source taxon included!
                    # source taxon will always appear at top of sorted drop list because it is not possible for taxon to occur more frequently
                    # with any other single taxon that it occurs with *all* other taxa and *alone* (which is how we define self cooccurrence)
                    targetp=tuple(sorted(targetp))
                    #print "drop of > %d percent (%.1f) after taxon %d combined with taxa %s" % (dropthresh, droppct, i, `targetp`)
                    allbigdrops[i].append((targetp, droppct))
                    # trying to add an adjustment to favor smaller partitions which are more likely to be correct
                    overallbigdropsl.append((targetp, droppct, droppct/len(targetp)))
                    if targetp not in overallbigdropsavg:
                        overallbigdropsavg[targetp]=[]
                    overallbigdropsavg[targetp].append(droppct)
            lastn=n
            ti+=1

    overallbigdropsl.sort(key=operator.itemgetter(1), reverse=True)

    # Compute average drops for a partition.
    # A given partition will appear once for each source partition within it.
    for p in overallbigdropsavg:
        v=overallbigdropsavg[p]
        a=sum(v)/float(len(v))
        overallbigdropsavg[p]=(a, v)        # replace the list of drops with (average, list of drops)
    overallbigdropsavgl=[(p, overallbigdropsavg[p][0], overallbigdropsavg[p][1]) for p in overallbigdropsavg]
    overallbigdropsavgl.sort(key=operator.itemgetter(1), reverse=True)     # sort on avg

    #print overallbigdropsavgl
    return overallbigdropsavgl, overallbigdropsl, allcdls

def DumpAndPlotBigDropsList(bdl):
    """ Diagnostic tool to graphically show the big drops list to help see where the breakpoint in correlations might naturally fall.
    """
    if 0:
        for p, d, dn in bdl:
            print p, d
    ys=[d for p, d, dn in bdl]
    xs=range(0, len(ys))
    f=pl.figure
    ax=pl.axes()
    ax.plot(xs, ys, 'r')
    # identify where TT partitions occur in this plot
    xv=0
    txs=[]; tys=[]
    nxs=[]; nys=[]
    for p, d, dn in bdl:
        if p in truetree:
            tys.append(d)
            txs.append(xv)
        else:
            nys.append(d)
            nxs.append(xv)
        xv+=1
    ax.plot(txs, tys, 'b.', label='True tree partitions')
    ax.plot(nxs, nys, 'g.', label='Non-true tree partitions')
    #pl.legend([patches[0][0], patches[1][0]], ["Extinct", "Still living at end"])
    pl.xlabel('Partition');
    pl.ylabel('Percent correlation drop after partition');

    pl.xlim((0, 50))
    # mark the previously identified 'markbkpt'th point as the breakpoint
    (p, ym, dn)=bdl[markbkpt]
    xm=markbkpt
    #arrowstyle='->'
    arrowstyle='fancy'
    #ax.annotate("Breakpoint", xy=(xm, ym), xycoords='data', xytext=(xm+20, ym+20), textcoords='offset points')
    ax.annotate("Estimated breakpoint (ntaxa-2 method)", xy=(xm, ym), xycoords='data', xytext=(xm+10, ym+100), textcoords='offset points',
        arrowprops=dict(arrowstyle=arrowstyle, connectionstyle="arc3", color='cyan'), color='cyan')
    # plot incompatibilities with simple puzzling
    ntc=[]
    xm=0
    for p, d, dn in bdl:
        if not FastCompatibleJustSynsets(p, ntc):
            ax.annotate("First incompatibility", xy=(xm, d), xycoords='data', xytext=(xm+30, ym+50), textcoords='offset points',
                arrowprops=dict(arrowstyle=arrowstyle, connectionstyle="arc3", color='red'), color='red')
            break
        else:
            ntc.append(p)
        xm+=1

    handles, labels = ax.get_legend_handles_labels()
    print handles
    print labels
    #ax.legend([handle for i,handle in enumerate(handles) if i in display], [label for i,label in enumerate(labels) if i in display])
    ax.legend()

    #ttpartsys=[d for p, d, dn in overallbigdropsl if p in truetree]
    pl.title('Correlation Drop Percentages')
    fn='drops.png'
    pl.savefig(fn); print "saving drops plot to", fn
    #pl.show()
    pl.close()
    #exit()

def SignalBuilderClean(sigs, allfreqsd, ntaxa, truetree, verboselevel=0, outgrouptaxa=[], taxanames={}, dolongbranch=False, omitpartitions=[], taxareversenames={}, justreturnbigdrops=False):
    """ Use HPA algorithm to return a set of compatible partitions that represent an estimated phylogeny given the raw partition frequences in allfreqsd dict.

        HPA method involves building a 'drop list' where each taxon is considered as source and we look at the next most highly correlated taxon,
        and then the next most highly correlated taxon, and so on, and build partitions from those progressive sets, and then combine the most compelling
        compatible set of those partitions together.

        Q: do we pay attention to our significance metric, that adjusts frequency by size, at all any more?
        A: nope. sigs dict is not even being referenced any more. everything based on correlation tables.
        Q: do we ever look at a reduced set of the most significant partitions, for example ntaxa*5, in cooccurrences or anywhere else?
        A: not anymore! now we look at all

        FUTURE: allow option to use significant vs. all partitions in cooccurrence step
                (do check on similarity of output)

        If justreturnbigdrops is True, we just return the sorted bigdropsavgl list (with internal taxa indices). This
        is useful when caller wants to look at intermediate results for some reason.
    """
    freqsd=allfreqsd

    lbpass=1
    lastLBset=set([])
    while True:     # keep constructing overallbigdropsl until no further LB taxa are identified
        dropthresh=0.0
        if verboselevel > 0: print "phase 1, pass %d, drop threshold %f" % (lbpass, dropthresh)
        # 'dropthresh' percent drop or greater is deemed important.
        # maybe for large datasets a threshold of 0.0 won't be workable
        # because we will have too many entries in our drop table?

        overallbigdropsavgl, overallbigdropsl, allcdls=BuildCorrelationDropList(ntaxa, freqsd, truetree, verboselevel, dropthresh=dropthresh, taxanames=taxanames)

        """
            a curve drawn of the percentage drops in the correlation table empirically shows a sharp bend around the point where
            we start to get non-TT partitions. we have tried to estimate this point in various ways, with an arbitrary 5% threshold,
            with LOF (local outlier factor), but it seems to work well to just put the breakpoint for LB taxa after the first ntaxa-2
            partitions. Any taxa that don't occur in any of the first ntaxa-2 partitions are designated LB taxa.
        """
        if 1:
            markbkpt=ntaxa-2
            if verboselevel > 1: print "using the ntaxa-2 criterion, the breakpoint for long branch check is", markbkpt

        if verboselevel > 2: # informational printout: {{{ showing variance and other info of different drops for each source taxon
            print "Average drop across all representations of a partition (for each source taxon):"
            print "Note: the x of y found means we found x source-taxa representations of that partition"
            print "Note: because of order of inclusion of taxa in partitions based on correlation levels, sometimes we don't find all partitions representated"
            for p, av, v in overallbigdropsavgl:
                if p in truetree: intt="*"
                else: intt=" "
                warn=float(len(v))/float(len(p))
                print intt, len(v), "of", len(p), "found, frac=", warn, "var=", numpy.var(v), len(v), "of", len(p), p, av         #, v   # v is the list of drops we are averaging
        #}}}

        if verboselevel > 2: # informational printout: {{{ showing 20 partitions with most significant sum percent drops
            print "the 20 partitions with the most significant sum percentage drops (after them) are:"
            for p, n, na in overallbigdropsl[:20]:
                if p in truetree: intt="*"
                else: intt=" "
                print intt, p, n      #, allbigdrops[p][1]
        #}}}

        if justreturnbigdrops:
            return None, overallbigdropsavgl

        if len(omitpartitions) > 0:
            omitpartitionsinternal=[]
            print "omitpartitions:"
            for p in omitpartitions:
                ip=tuple(sorted([taxareversenames[ti] for ti in p]))
                omitpartitionsinternal.append(ip)
                print p, ip
            print "the 20 partitions with the most significant sum percentage drops (OMITTING previously selected partitions) are:"
            # convert omitpartitions to internal index numbers for ease of comparison
            print "omitpartitions with taxa names converted to indices"
            print omitpartitionsinternal
            print "the most significant remaining:"
            #for p, n, na in overallbigdropsl:
            for p, av, v in overallbigdropsavgl:
                if p in omitpartitionsinternal:
                    print "SKIP", sorted([taxanames[ti] for ti in p]), p
                else:
                    print sorted([taxanames[ti] for ti in p]), p
            print "FIXME--want to return this data not exit here, but this is just test code"
            exit()

        """
            At this point we have the big drops for each source taxon i in allbigdrops[i], sorted by drop %
            And we have a global list of big drops in overallbidropsl, sorted by drop %
        """
        # now identify "long branch taxa":
        # find which taxa are not present in any of the drops above certain point in droplist
        # we were using a percentage drop as the threshold, but now we just look at a certain number of items in the drop list
        # breakpointthresh should probably be computed from the data, maybe using a density approach (see email)
        #breakpointthresh=5.0
        if verboselevel > 1: "finding taxa that are not present in any partition above breakpoint (the first %d partitions) in overall average drop list", markbkpt
        notfound=[]
        foundtaxa=set([])
        alltaxa=set([t for t in range(0, ntaxa)])
        for p, n, na in overallbigdropsavgl[:markbkpt]:
            if verboselevel > 1: print p, n
            #if n < breakpointthresh:
            #    print "stopping looking for breaks on source taxon", i, "because breakpointthresh reached"
            #    break
            for t in p:
                if verboselevel > 1: print "found", t, "in", p
                foundtaxa.add(t)

        LBset=alltaxa-foundtaxa

        if verboselevel > 0: print "LBINFO: on pass %d, the following taxa ('LB taxa') were not found in any big drops in that pass:" % lbpass, LBset

        if not dolongbranch:
            if verboselevel > 0: print "Long branch processing has been disabled! Nothing will be done with any LB taxa found."
            break

        if LBset==lastLBset:
            break

        if verboselevel > 2: # informational printout {{{
            print "big drops :20"
            for p, n, na in overallbigdropsavgl[:20]:
                print p, n, na
        #}}}

        """
            Optionally remove partitions containing long branch taxa.
            In an earlier implementation, when removing LB taxa we would just pull any partitions containing the LB taxa
            out of the overallbigdropsl and proceed. Now, we remove any partitions containing the LB taxa from the *entire*
            list of partitions and *recompute* the overallbigdropsl.
        """
        if 0:       # old way, removing from overallbigdropsl
            totalnparts=len(overallbigdropsavgl)
            overallbigdropsavgl=[(p, n, na) for p, n, na in overallbigdropsavgl if len(set(p).intersection(LBset)) == 0]
            lbfilterednparts=len(overallbigdropsavgl)
            if verboselevel > 1:
                print "LBINFO: after filtering out partitions containing long branch taxa from overallbigdropsl, we went from %d partitions to %d" % (totalnparts, lbfilterednparts)

        if 1:       # new way, removing from source partitions in freqsd and then recomputing overallbigdropsl
            ls=len(freqsd)
            freqsd={k:v for k,v in freqsd.items() if len(set(k).intersection(LBset))==0}
            es=len(freqsd)
            if verboselevel > 1: 
                print "LBINFO: freqsd went from %d to %d partitions after removal of partitions containing LB taxa" % (ls, es)
            lastLBset=LBset

        lbpass+=1

    # at this point we have all LB taxa removed and a list of average correlation percent drops

    if 0:
        if verboselevel > 1:
            print "we are looking at the bigdropslist where each source taxon is considered separately"
        bdl=overallbigdropsl
    else:
        if verboselevel > 1:
            print "we are looking at the average bigdropslist where each drop for each source taxon in partition are averaged"
        bdl=overallbigdropsavgl

    if 0: DumpAndPlotBigDropsList(bdl)   # informational printout and plot of correlation % drops

    tc=[]

    if verboselevel > 0 and len(outgrouptaxa) > 0:
        outgrouptaxa=tuple(sorted(outgrouptaxa))
        print "seeding trialclique with the following outgroup taxa:", outgrouptaxa
        outgrouptaxacomplement=tuple(sorted([n for n in range(0, ntaxa) if n not in outgrouptaxa]))
        print "seeding trialclique with the following outgroup taxa complement:", outgrouptaxacomplement
        tc.append(outgrouptaxa)
        tc.append(outgrouptaxacomplement)

    badadds=0
    stopatfirstincompatibility=False
    if verboselevel > 1: print "starting initial puzzling"
    if verboselevel > 1:
        if len(truetree) == 0: print "NOTE: ignore non TT errors because we have no true tree to compare against"
    # puzzle from biggest drop on down, until we reach some exit condition (might be first incompatibility encountered, full set)
    ntry=0
    nnonleafpartitions=0    # number of non leaf partitions added to trialclique
    for fp, n, na in bdl:       # bdl is either average big drop for all source taxa in partition, or just full, unaveraged list
        if nnonleafpartitions==ntaxa-2:
            if verboselevel > 1: print "we have reached ntaxa-2 partitions in trial clique and will now exit"
            break
        if fp in tc:
            continue
        if FastCompatibleJustSynsets(fp, tc):
            tc.append(fp)
            if len(fp)>1: nnonleafpartitions+=1
            if fp in truetree:
                if verboselevel > 1: print "correctly added compatible", fp, "to trialclique"
            else:
                if verboselevel > 1:
                    if len(truetree)>0: print "ERROR added compatible but non TT partition", fp, "to trialclique"
                badadds+=1
        else:
            if fp in truetree: intt="*"
            else: intt=""
            if stopatfirstincompatibility:
                if verboselevel > 1:
                    print "we got our first incompatibility on", `fp`+intt, "and hpa is stopping now"
                break
        ntry+=1
        if ntry==40 and verboselevel > 1: print "will no longer show correct rejections due to incompatibility with trialclique"

    if verboselevel > 1:
        print "length of trialclique is now", len(tc)
        print "trialclique is now:", tc
        print "%d bad adds so far (valid only if we know the true tree!)" % badadds
        print "done"

    return tc, overallbigdropsavgl, allcdls

def TaxaAtOrBelowNode(thisn):
    """ Return list of taxa names AT OR BELOW a given node in a dendropy tree. """
    if thisn is None:
        return []
    if thisn.taxon is not None:
        return [thisn.taxon.label]
    else:
        bh=[]
    for n in thisn.child_nodes():
        #if n.taxon is not None:
        #    bh.append(n.taxon.label)
        bh.extend(TaxaAtOrBelowNode(n))
    return bh

def TaxaBelowNode(thisn):
    """ Return list of taxa names BELOW and NOT INCLUDING a given node in a dendropy tree. """
    bh=[]
    for n in thisn.child_nodes():
        if n.taxon is not None:
            bh.append(n.taxon.label)
        bh.extend(TaxaBelowNode(n))
    return bh

def TupleAfy(nt):
    """
        Given nt, a newicky-nested structure that might have some lists and tuples inside, get rid of any lists that might be present and convert them to tuples

        FUTURE:  maybe should sort
    """
    g=[]
    for p in nt:
        if type(p)==type([]):
            g.append(TupleAfy(p))
        else:
            g.append(p)
    return tuple(g)

def NameMapNewick(t, tmap, sofar=[]):
    """
        Given t, in my tupleish representation WITH NO BRANCH LENGTHS, traverse it and return rebuilt version with taxa names substituted from tmap.
    """
    #print len(t)
    #print "t[0]:"
    #print t[0]
    #print "t[1]:"
    #print t[1]
    b=[]
    for p in t:
        if p in tmap:
            #print "GOT", tmap[p]
            b.append(tmap[p])
        else:
            #print p, "not in tmap, recursing"
            #print "tmap:", tmap
            if len(p) > 0:
                b.append(NameMapNewick(p, tmap))
    return tuple(b)

def InsertInNewick(nt, clique):
    """
        Helper function for converting clique representation to Newick.

        recursively split up partition in nt
        nt stored as list of lists so we can modify it
        example:
        InsertInNewick([1, 2, 3, 4, 5], [2, 3, 5]) will become:
            [1, 4 [2, 3, 5]]
        InsertInNewick([1, 4 [2, 3, 5]], [3, 5]) will become:
            [1, 4, [2, [3, 5]]
    """
    justtaxa=set([i for i in nt if type(i) is not type([])])
    justsublists=[i for i in nt if type(i) is type([])]
    if set(clique) <= justtaxa:       # some items may be taxon indices and some lists.  can't put lists into a set, have to exclude them!
        #print "before:", nt
        #print "to insert:", clique
        for p in clique:        # remove all the items present in the to-be-inserted clique from the current level, so we can subpartition it
            #print "deleting taxon", p, "at index", nt.index(p), "from", nt
            del nt[nt.index(p)]
        nt.append(clique)
    else:
        # from the sublists at our current level, look at each of them to see if our new clique can be a
        # subpartition of one of them
        for l in justsublists:
            InsertInNewick(l, clique)

def NewickACize(taxanames, cliques):
    """
        Convert cliques, a list of partitions that we think represent structure of the phylogenetic tree, to a nested newick style representation.

        Method: start with a full list of taxa, then go through the list of cliques, from largest to smallest,
            and put a nested sub-list inside the larger list, for just the taxa in that clique
            if there is a conflict (a clique is not a subset of an existing sublist), then just skip it
    """
    # we have to look at longest length subtrees first because of how NewickaCize works
    cliques.sort(key=len, reverse=True)
    #print "sorted:"
    #print cliques
    comptree=[tn for tn in taxanames]
    #print "NewickACize", comptree, "with", cliques
    for clique in cliques:
        if len(clique)>1:       # this check is new. if single taxon partitions are included in cliques, this would create extra tuple elements we didn't want
            InsertInNewick(comptree, list(clique))
    # convert from lists to tuples?
    #print "comptree before TupleAfy:"
    #print comptree
    comptree=TupleAfy(comptree)
    #print "after TupleAfy:"
    #print comptree
    return comptree

def PopulateDendropyBranchLengths(n, bl):
    """
        Set branch lengths in dendropy tree n based on what is given in bl.

        If bl is a float or int, set all edge lengths to that value

        If bl is a dict, set edge lengths to the the value in the dict using the sorted tuple of taxa below node as key

        This does not (re)label nodes at the moment. could easily do that.
    """
    lb=tuple(sorted(TaxaAtOrBelowNode(n)))
    #print "looking for", lb
    if type(bl)==type({}):
        if lb in bl:       # bl may not be a dict or default dict! may be int or float or something else even
            n.edge.length=bl[lb]
        else:
            # the partition may not have been found; default freq/edge len to 0
            #print "no edge length!"
            #exit(0)
            n.edge.length=0
    else:
        n.edge.length=bl
    for c in n.child_nodes():
        #if (len(c.child_nodes())!=0):     # if child is not a terminal node
        PopulateDendropyBranchLengths(c, bl)
        #else:
        #    print "terminal"

def IdentifyVPs(seqdata, usenumpy=True, verboselevel=0, alphabet=['C', 'G', 'A', 'T']):
    """
        Enumerate the variable positions in the sequence data seqdata.

        Returns a list where each element of the list corresponds to a variable position locus,
        and is a list of tuples of indices which correspond to 0-based internal taxa numbers for each character in the alphabet.

        e.g., if sequences for taxa T1, T2, T3 are:
        T1: AAC
        T2: AAA
        T3: CGG
        then we return [[(2,), (0, 1)], [(2,), (0, 1)], [(0, ), (2,), (1,)]]

        DOES NOT return empty tuple for characters not found (used to do that, but no longer).
        Only the occurrence of the same character matters; thus the first two VPs appear the same in the results above.

        seqdata is a dendropy DNACharacterMatrix

        FUTURE: maybe return partition frequency dict too?
    """
    # NEW: we return single-taxon partitions too!
    # FUTURE: check for X or N or x or n or - characters for "unknown"
    # FUTURE: check for upper and lower case inputs
    # FUTURE: pad short sequences
    #   maybe use DnaCharacterMatrix fill() to pad lengths?

    if (usenumpy):
        # all numpy approaches I have tried end up somewhat slower than my basic python version
        # for the basic vp indentification;
        # however we are doing more work here along the way (computing tuples for partitions) than in python version
        # so this might end up a win
        #
        # convert dendropy DnaCharacterMatrix to numpy array
        #alldata=numpy.empty((rows, cols), dtype='|S1')
        # this would leave a vector of strings:
        #alldata=numpy.array([v.symbols_as_string() for k, v in seqdata.items()])
        # print alldata.shape       # this will give (20,) if input has 20 taxa.

        # but we want a two dimensional array of characters:
        alldata=numpy.array([list(v.symbols_as_string()) for k, v in seqdata.items()])
        rows, cols=alldata.shape
        #if verboselevel > 0: print "converted to numpy array at time %.1f" % (time.time()-ts)
        if verboselevel > 0: print "seqlen %d ntaxa %d" % (cols, rows)

        # new approach that may be faster. operate on entire array?
        # think about redoing this
        #for c in alphabet:
        #    m=where(alldata == c)
        #    print type(m)
        #    print m[0]
        #    print m[1]
        #exit()

        res=[]
        #idx=0
        # this is surprisingly slow! slower than my old non-numpy method
        # but we are going to do it anyway because ... ? (explain that here)
        for i in range(cols):           # for each variable position:
            #locus=alldata[:, i]
            #print tuple(where(locus=='A'))
            # NOTE: numpy.where returns a tuple of two ndarrays; we only want the first value, which is the indices

            #thislocus=[tuple(where(alldata[:, i] == c)[0]) for c in alphabet]
            thislocus=[]
            for c in alphabet:
                # what about the issue of single-taxon partitions? do we want to track them?
                # in that context, consider cases like AAAAT vs AABBT:
                # in AAAAT case, we would get partition (0, 1, 2, 3) and could infer a single-taxon 'partition' (4)
                # in AABBT case, we would get partitions (0, 1) and (2, 3) and could not infer single-taxon 'partition' (5) unless
                # we continued to track the two partitions as a (locus) unit, which we do up to a point, but not always
                # a locus can have a partition for each item in the dictionary (ACGT for DNA data)
                tp=tuple(where(alldata[:, i] == c)[0])        # where returns 0-based indices
                # if a character is not present, we will get empty tuple! e.g. AC will return (0,) for A and (1,) for C and () for T and G
                if len(tp)>0: thislocus.append(tp)      # NEW: don't return empty tuple if char not found
                # that is not helpful. may be faster to remove later though...
                #if len(tp)>1: thislocus.append(tp)  # don't add single-taxon partitions
                #if len(tp)==1: print "complement"
                #thislocus.append(tp)               # do add single-taxon partitions
                #FIXME: check for >1 element in list, otherwise it's not a VP
            #if len(thislocus)>1: res.append(thislocus)     # NOTE! if we don't allow single-taxon partitions, then this check will prevent a size ntaxa-1 partition from being counted!

            # NEW! thought: add an artificial one-taxon partition if and only if the only other partition found is a single complementary ntaxa-1 partition (that is, not a combination of 2 or 3 other partitions!
            # this results in many fewer single-taxon partitions than otherwise, hopefully eliminating spurious single-taxon partitions
            #if len(thislocus)==1 and len(thislocus[0])==(rows-1):
            #    rempart=tuple(tn for tn in range(0, rows) if tn not in thislocus[0])
            #    #print "adding remainder single-taxon complement partition", rempart
            #    thislocus.append(rempart)

            res.append(thislocus)
            #idx+=1
            #if verboselevel>0 and (idx%10000)==0: print idx, "of", cols, time.time()-ts, "secs", idx*100/cols, "%"

        #if verboselevel > 0: print "all partitions created, time %.1f" % (time.time()-ts)
        #exit()
    else:
        print "must use numpy for now. exiting"
        exit()
    return res

def hpa(alignmentfn=None, seq=None, resampletype=0, nresamples=0, verboselevel=0, outgroup=[], dolongbranch=True, omitpartitions=[], taxareversenames={}, justreturnbigdrops=False, truetree=[], supportvalues=None, restcs=None, branchlengths=None):
    """
        Load the alignment in FASTA format from filename 'alignmentfn' OR given in 'seq' and compute and print an estimated phylogenetic tree in newick format using the HPA algorithm, optionally resampling.

        if seq is not None, instead of loading alignment from file alignmentfn, alignment will be passed in seq
            (a Dendropy DNACharacterMatrix)

        if justreturnbigdrops is True, we will simply return the bigdrops list (partition, avg, listofavgeragedvalues)

        otherwise returns tc, rnt, ct:
            tc: TUPLE of partitions in trialclique using original taxa names,
            rnt: nested tuple representation with original taxa names,
            dt: dendropy tree representation with original taxa names

        If branchlengths is a dict instead of None, estimated branch lengths will be returned in dt and also in
            dict with key being partition (FIXME: using internal taxa names? external taxa names?)

        if omitpartitions is supplied, it is a list of partitions (using real taxon names) that will be skipped
            in printout of partitions found with big drops. This is useful for filtering out the non-reticulated
            result partitions when studying reticulated data sets.
            THEN WE EXIT()! This feature is not really used now.

        if resampletype==0, we won't do any resampling/bootstrapping. Each vp is used exactly once.
        if resampletype==-1, we will do resampling/bootstrapping by randomly selecting WITH REPLACEMENT 'n' of
            the variable positions (vps) from the full set of vps, where 'n' is the number
            of vps found. Note that each vp can contain two, three, or four partitions.
            Some vps may therefore be included multiple times and some vps not examined at all.
        if resampletype>0, we will use do 'nresamples' of size n='resampletype', samples with replacement.

        repeat this 'nresamples' times.
            if nresamples > 0 and resampletype==0, then we should get the same result each time,
            but that might be useful for testing purposes

        when resampling is requested, the actual result will be based on the full data set sampled once.
            support values are returned in FIXME

        if supportvalues is a dict, then support values for the partitions will be tallied and placed in there.
            key is partitions, value is frequency.
            frequency values will range from 1 to nresamples and represent the number of resamples where the partition
            was placed in the trialclique
            remember though that the return values (tc, rnt, dt) are computed using the full set of data, non-resampled

        The dendropy tree returned in dt will use taxon_namespace from the incoming sequence

    """
    ts=time.time()

    if seq is None:
        if verboselevel > 0: print "start loading at time %.1f" % (time.time()-ts)
        seqdata=dendropy.DnaCharacterMatrix.get(file=open(alignmentfn), schema="fasta")
        if verboselevel > 1: print "seqdata:", seqdata
        if verboselevel > 0: print "done loading at time %.1f" % (time.time()-ts)
        if verboselevel > 0: print "start identify variable positions at time %.1f" % (time.time()-ts)
    else:
        seqdata=seq

    if verboselevel > 0:
        print "WARNING: currently polymorphic or ambiguous states in input sequences are simply ignored"

    taxanames={}        # key is internal number, value is 'official' name
    tn=0
    for k, v in seqdata.items():
        taxanames[tn]=k.label
        tn+=1
    taxareversenames={v: k for k, v in taxanames.items()} # key is 'official' name, value is internal taxon index number
    if verboselevel > 1:
        print "taxon mapping internal number -> official name"
        for k in taxanames.keys():
            print k, "->", taxanames[k]
    ntaxa=len(seqdata.taxon_namespace)

    # we need only compute the VPs (variable positions) once. each element in res will be for one column slice of input sequences.
    # if we are doing resampling, we will resample among these elements of res (in other words, we 
    # until we decide what to do with ambiguous characters (X or - or whatever) those are simply ignored;
    # a column that has all ambigious data will have empty list for that column
    res=IdentifyVPs(seqdata, usenumpy=True, verboselevel=verboselevel)

    outgrouptaxa=[]
    for tn in outgroup:
        origtaxanames=[taxanames[k] for k in taxanames]
        if tn in origtaxanames:
            outgrouptaxa.append(taxareversenames[tn])
        else:
            print "taxon", tn, "is requested to be in outgroup but it isn't found in the input data"
            print "valid taxa:", origtaxanames
            exit(-1)

    if verboselevel > 1:
        print "outgrouptaxa internal indices:", outgrouptaxa

    tempsupportvalues={}        # built with internal taxa names, later copied to supportvalues
    temprestcs={}                  # dict of resulting trialcliques (keys) with frequency as values

    # first do any resamples, before computing actual hpa result on full nonresampled input data set (below)
    for sampn in range(0, nresamples):
        if verboselevel > 0: print "start resample %d at time %.1f" % (sampn, time.time()-ts)
        if resampletype==-1:
            allvps=[]
            if verboselevel > 1: print "resampling with replacement %d samples (of %d)" % (len(res), len(res))
            for i in range(0, len(res)):
                allvps.append(random.choice(res))
        elif resampletype>0:
            allvps=[]
            for i in range(0, resampletype):
                allvps.append(random.choice(res))
        else:
            allvps=res
        if verboselevel > 0: print "done resample at time %.1f" % (time.time()-ts)

        # now that we know the actual loci we will be using, we can compute partition frequencies based on that
        # NOTE: if a locus is selected in the sample, ALL partitions at that locus are selected

        allfreqsd=defaultdict(int)
        # we will have 0-4 partitions at each locus (assuming ACGT data).
        # until we decide what to do with ambiguous characters (X or - or whatever) those are simply ignored;
        # a column that has all ambigious data will have empty list for that column
        # here we pull all partitions out of each locus (column slice in sequence data) into one big mass:
        for pl in allvps:
            for p in pl:
                allfreqsd[p]+=1

        if verboselevel > 0: print "start signalbuilder at %.1f" % (time.time()-ts)
        sigs=None       # use all partitions, not just significant ones
        tc, oavgbdl, allcdls=SignalBuilderClean(sigs, allfreqsd, ntaxa, truetree, verboselevel=verboselevel, outgrouptaxa=outgrouptaxa, taxanames=taxanames, dolongbranch=dolongbranch, omitpartitions=omitpartitions, taxareversenames=taxareversenames, justreturnbigdrops=justreturnbigdrops)

        # accumulate occurrences count of partitions in tempsupportvalues; moved to supportvalues for return after all resamples
        if type(supportvalues) is type({}):
            for p in tc:
                if p in tempsupportvalues: tempsupportvalues[p]+=1
                else: tempsupportvalues[p]=1

        # tc is a list of partitions using internal taxon numbers; convert to tuple so we can use it as key and keep track of frequency
        if type(restcs) is type({}):
            tco=tuple(sorted(tc))
            if tco in temprestcs: temprestcs[tco]+=1
            else: temprestcs[tco]=1

        if verboselevel > 0: print "done signalbuilder at %.1f" % (time.time()-ts)

    # We now have a frequency count of all unique tree results in temprestcs ("temporary resulting trialcliques") from
    # all the resamples we have conducted.
    # Convert result (in the frequency count dict) to external taxa names.
    # We no longer intend to return the most common of those discovered trialcliques as the overall algorithm result,
    # instead we will return the non-resampled result where each input column (vp essentially) is considered exactly
    # once. But we will track the most common resampled tree result in case we want to look at that again in the future.
    mostcommontc=None
    if restcs is not None:
        print "temprestcs:", temprestcs
        maxfreq=0
        for rstc in temprestcs.keys():
            print "resampled tc with internal taxa numbers:", rstc
            tco=tuple(sorted([   tuple(sorted([taxanames[tn] for tn in p]))    for p in rstc]))
            print "tco with original taxa names:", tco
            thisfreq=temprestcs[rstc]
            restcs[tco]=thisfreq
            if thisfreq > maxfreq:
                maxfreq=thisfreq
                mostcommontc=rstc     # will be converted to original taxanames below, just before function return
        print "most frequent trialclique result in resampling (%d occurrences):" % (maxfreq), mostcommontc

    if supportvalues is not None:
        # convert supportvalues taxa names from internal names to real names, and copy to return dict 'supportvalues'
        for p in tempsupportvalues.keys():
            newp=tuple(sorted([taxanames[tn] for tn in p]))
            supportvalues[newp]=tempsupportvalues[p]
        if verboselevel > 1: print "supportvalues:", supportvalues

    if verboselevel > 1:
        print "Resamples (if any) are complete. Starting HPA on original data set."

    ################################################################################
    # now compute the actual (non-resampled) hpa result on original input:
    allfreqsd=defaultdict(int)
    for pl in res:      # res is the full set of vps from IdentifyVPs, based on input sequence data
        # within each vp there will be from 0 to 4 partitions (if we have ACGT characters)
        # until we decide what to do with ambiguous characters (X or - or whatever) those are simply ignored;
        # a column that has all ambigious data will have empty list for that column
        #print len(pl), pl
        for p in pl:
            #if len(p)==0:
            #    print "zerolen:", pl, [len(ap) for ap in pl]
            #    exit()
            allfreqsd[p]+=1
    sigs=None       # use all partitions, not just significant ones

    if 0:
        # diagnostic code to see why certain partitions were not being found
        # see email exchange with Hoelzer 6/2017
        print "new work"
        allfreqsl=sorted([(allfreqsd[k], k) for k in allfreqsd.keys()], reverse=True)
        for i in range(0, 100):
            c, p = allfreqsl[i]
            print allfreqsl[i], "len", len(p),
            if len(p) > ((9 * ntaxa) / 10):
                # print complement if it will be hard to see what is missing:
                cp=[tn for tn in range(0, ntaxa) if tn not in p]
                print "missing:", cp
            else:
                print
        exit()

    tc, oavgbdl, allcdls=SignalBuilderClean(sigs, allfreqsd, ntaxa, truetree, verboselevel=verboselevel, outgrouptaxa=outgrouptaxa, taxanames=taxanames, dolongbranch=dolongbranch, omitpartitions=omitpartitions, taxareversenames=taxareversenames, justreturnbigdrops=justreturnbigdrops)

    if justreturnbigdrops:
        # convert to real taxa names (SignalBuilder only knows internal taxa numbers, not real names)
        pl=[]
        for p, avg, avgels in oavgbdl:
            pl.append((tuple(sorted([taxanames[i] for i in p])), avg))
        return pl

    print "tc result:", tc

    rnt=NewickACize(taxanames, tc)
    # rnt is tuple like: (((17, (1, ((13, (14, (10, (15, (9, 16))))), (19, (18, (11, 12)))))), (0, 3)), ((2, (8, (6, 7))), (4, 5)))
    rnt=NameMapNewick(rnt, taxanames)
    # rnt now is like: ((((('T6', ('T4', 'T5')), ('T8', 'T7')), ('T1', ('T2', 'T3'))), ('T11', ('T9', 'T10'))), (('T16', ('T20', ('T19', ('T17', 'T18')))), ('T12', ('T13', ('T14', 'T15')))))

    # newick requires ; at end, so we append it here:
    # we don't specify a taxon_namespace so one will be created
    dt=dendropy.Tree.get_from_string(`rnt`+';', "newick", taxon_namespace=seqdata.taxon_namespace)

    ################################################################################
    # if requested, try to estimate branch lengths.
    # FIXME: this code is pretty awkward.
    # it would be best to do this BEFORE converting tree from internal taxa numbers to external names,
    # because branch length will be determined (for now) based on raw frequency of partition occurrence in allfreqsd which uses
    # internal taxa numbers, however for annoying reasons of data conversion complexity from my internal types to dendropy tree we
    # are plugging in branch lengths AFTER, which requires a conversion of the frequency count dict into the external namespace.
    # FIXME at some point.
    # one easy thing: rather than converting ENTIRE allfreqsd into external namespace,
    # it would be more efficient to just convert the ones we need
    if branchlengths is not None:
        if verboselevel > 0: print "estimating branch lengths"
        #print "FIXME: check if single taxons are present in allfreqs"
        #print "tc:", tc
        if 0: # {{{
            # hoelzerish approach? only works if we add single-taxon partitions
            # add single-taxon partitions
            ntc=[] #[p for p in tc]
            for tn in range(0, ntaxa):
                ntc.append((tn,))
            #print "trialclique with single-taxon partitions added:", ntc
            # build up partitions from single taxon upwards:
            # sort tc by length to speed up search process
            print "building up"
            stc=sorted(list(tc), key=len)
            cnts=[]
            for p1 in ntc:
                for p2 in stc:
                    if len(set(p1).intersection(p2)) != 0:
                        print p2, "contains", p1
                        if p2 in cnts:
                            print "(already added)"
                            break
                        else:
                            cnts.append(p2)
                            print "added", p2
                            print "freq of", p1, allfreqsd[p1]
                            print "freq of", p2, allfreqsd[p2]
                            break
                #containsp1=[p in tc if set(p).intersection
                #ne=tuple(sorted(set(p1)
            #for p1 in combinations(tc, 2):
            #    print p1
            exit()
        # }}}
        if 0:       # algorithm 2: branch length is simply partition frequency
            if verboselevel > 0: print "using naive technique where branch length is simply set to the raw freqeuency of the partition"
            # convert frequency dict to external namespace. somewhat awkward. improve efficiency if the approach seems worthwhile
            allfreqsdorigname={}
            hacklengths=False    # ad hoc adjustment for underestimation based on number of taxa underneath
            for k in allfreqsd.keys():
                #print "k", k
                ok=tuple(sorted([taxanames[tn] for tn in k]))
                #if len(k)==1: print "ONETAXON", k, ok, allfreqsd[k]
                if hacklengths:
                    #allfreqsdorigname[ok]=float(allfreqsd[k])*sqrt(len(k))
                    allfreqsdorigname[ok]=float(allfreqsd[k])/sqrt(len(k))
                else:
                    allfreqsdorigname[ok]=allfreqsd[k]
            #print "converted frequency dict:"
            #for k in allfreqsdorigname.keys():
            #    print k, allfreqsdorigname[k]
            PopulateDendropyBranchLengths(dt.seed_node, allfreqsdorigname)
        if 1:       # algorithm 3. see emails. average branch length to different to each leaf in sister clade, essentially
            if verboselevel > 0: print "new branch length calculation method"
            # convert frequency dict to external namespace. somewhat awkward. improve efficiency if the approach seems worthwhile
            blens={}
            tc.sort(key=len)
            for p in tc:
                #print
                print "working on partition", [taxanames[tn] for tn in p]
                parent=[pp for pp in tc if (len(set(pp) & set(p)) > 0) and len(set(pp) - set(p)) > 0]       # larger partitions containing p
                parent.sort(key=len)
                if len(parent)==0:
                    #print "no parent found, will take other largest clade as sister"
                    # in that case, the sister clade will be identified as the largest partition that has no common elements with p
                    for sp in reversed(tc):
                        if len(set(p) & set(sp)) == 0:
                            break
                else:
                    sp=tuple(set(parent[0]) - set(p))
                print "sister partition:", [taxanames[tn] for tn in sp]
                toavg=[]
                for tn1 in p:
                    tcdl=allcdls[tn1]
                    print "correlation drop table for taxon %s:" % taxanames[tn1]
                    print [(taxanames[tn], d) for tn, d in tcdl]
                    lookfort=set(p)
                    # what we should find in the drop table are the taxa in p,
                    # followed by the taxa in sp
                    # move through drop list until we reach the full source partition, to find the 'starting frequency count'
                    print "going through drop list until source partition is found"
                    nfound=0
                    ntofind=len(sp)
                    for tn, d in tcdl:
                        if tn not in p:
                            if verboselevel > 1: print "WARNING: encounted out of order taxon %s which is not in source partition" % taxanames[tn]
                        else:
                            lookfort.remove(tn)
                            if len(lookfort)==0:
                                #print "found all taxa in source partition at frequency", d
                                #print "now looking for target taxa in sister partition"
                                break

                    # now we know the frequency at which to start the drop pct calculation (the frequency at which the entire source partition is accounted for)
                    # since the order of taxa in the drop list may not agree with the order of taxa in accepted partitions, we start again from the
                    # beginning of the drop list to look for sister taxa
                    for stn, sd in tcdl:
                        if stn in sp:
                            if d==0:
                                print "WARNING: base frequency of source taxon is 0. Maybe look into this"
                            else:
                                pd=(float(d)-float(sd))/float(d)
                                print "found sister taxon %s with frequency %d, pct drop %f" % (taxanames[stn], sd, pd)
                                toavg.append(pd)
                            nfound+=1
                            if nfound==ntofind: break

                print "the drops to average:", toavg
                a=sum(toavg)/len(toavg)
                tp=tuple(sorted([taxanames[tn] for tn in p]))
                print "the average is %f which will be stored as edge length for partition" % a, tp
                blens[tp]=a

            if verboselevel > 1:
                print
                print "the final list of partitions and branch lengths into those nodes:"
                sp=blens.items()
                for v in sp:
                    print v
            PopulateDendropyBranchLengths(dt.seed_node, blens)
        if 0:
            if verboselevel > 0: print "new approach where branch length set to avg pct drop"
            # must convert big drop list to external namespace. somewhat awkward. improve efficiency if the approach seems worthwhile
            justavgdrops={}
            lens={}
            for k, a, aels in oavgbdl:
                ok=tuple(sorted([taxanames[tn] for tn in k]))
                justavgdrops[ok]=a
                ll=len(ok)
                if ll not in lens:
                    lens[ll]=1
                else:
                    lens[ll]+=1
            print justavgdrops
            print "distribution of lengths of partitions in avgbigdropl:"
            print sorted(lens.items())
            exit()
            PopulateDendropyBranchLengths(dt.seed_node, justavgdrops)
            print dt

    else:
        # fill in dummy branch lengths so some tools (e.g. trex robinson-foulds metric) works
        PopulateDendropyBranchLengths(dt.seed_node, 1.0)

    ################################################################################
    # convert tc to original taxa names from internal indices
    # NOTE: taxa are sorted strictly lexicographically, so 'T10' will appear before 'T8'! but the sort order will be unique and consistent
    #otc=[]
    #for p in tc:
    #    otc.append(sorted([taxanames[ti] for ti in p]))
    #tc=otc
    tc=tuple(sorted([   tuple(sorted([taxanames[tn] for tn in p]))    for p in tc]))

    # returns:
    #     TUPLE of partitions in trialclique using original taxa names,
    #     nested tuple representation with original taxa names,
    #     dendropy tree representation
    return tc, rnt, dt

if __name__=="__main__":
    """
        To do a very basic test with data included in this file, try:
            python hpa.py --quicktest
        or:
            python hpa.py --quicktest --nresamples 1000 --resampletype -1 --supportvalues
        
    """
    fn=None
    resampletype=0      # by default use each vp once and only once; no resampling
    nresamples=0
    verboselevel=1      # 1 is a minimal amount of helpful information. 0 is totally quiet. 2 and higher have more descriptive info
    outgroup=None
    outtree=None
    dolongbranch=False
    dosv=False
    i=1

    testseq=None

    # support arg=value syntax?

    # sample invocations
    # do 1000 resamples of 10000 randomly selected (with replacement) positions:
    # ./hpa.py -f pyvolverun500/simulated_alignment.fasta --resampletype 10000 --nresamples 1000

    # test with outgroup:
    # ./hpa.py -f pyvolverun500/simulated_alignment.fasta --outgroup "T20 T9"

    while i < len(sys.argv):
        arg=sys.argv[i]
        if arg=='-t' or arg=='--treeout':
            i+=1
            outtree=sys.argv[i]
            print "output newick tree to:", outtree
            i+=1
        elif arg=='--lb':
            i+=1
            dolongbranch=True
            print "long branch logic turned on"
        elif arg=='-o' or arg=='--outgroup':
            i+=1
            outgroup=sys.argv[i]
            print "got outgroup", outgroup
            i+=1
        elif arg=='-f' or arg=='--infile':
            i+=1
            fn=sys.argv[i]
            i+=1
        elif arg=='--resampletype':
            i+=1
            resampletype=int(sys.argv[i])
            i+=1
        elif arg=='--nresamples':
            i+=1
            nresamples=int(sys.argv[i])
            i+=1
        elif arg=='--verboselevel':
            i+=1
            verboselevel=int(sys.argv[i])
            i+=1
        elif arg=='--supportvalues' or arg=='-sv':
            print "will save support values"
            dosv=True
            i+=1
        elif arg=='-h' or arg=='--help':
            print '--lb do long branch check (default is not to)'
            print '--verbose <verboselevel from 0 to 10'
            print '--resampletype <0 = resamples, -1 = full size resampling with replacement, n = sample n vps with replacement>'
            print '-o <space or comma separated list of taxa names defining an outgroup>'
            i+=1
        elif arg=='--quicktest':
            testdata={
                "t1": "TTCCAA",
                "t2": "TTCCAG",
                "t3": "TTCGAG",
                "t4": "TTAGAG",
            }
            testseq=dendropy.DnaCharacterMatrix.from_dict(testdata)
            i+=1
        else:
            if len(sys.argv)!=i+1:
                print "bad args"
                exit()
            # otherwise we will assume filename is remaining arg
            fn=sys.argv[1]
            i+=1

    if fn==None and testseq==None:
        print "-f <input filename in FASTA format>"
        exit(-1)

    if outgroup is not None:
        if ',' in outgroup:
            outgroup=outgroup.split(',')
        else:
            outgroup=outgroup.split(' ')
    if outgroup is None:
        outgroup=[]

    if verboselevel > 0:
        print "infile", fn
        print "outgroup", outgroup
        print "resampletype", resampletype
        print "nresamples", nresamples
        print "verboselevel", verboselevel
        print "dolongbranch", dolongbranch
        print "dosupportvalues", dosv

    #tc, newicknobranchlens, dt=hpa(fn, resampletype=-1, nresamples=100, verboselevel=1)
    bl={}
    sv=None
    if dosv:
        sv={}

    tc, newicknobranchlens, dt=hpa(fn, seq=testseq, resampletype=resampletype, nresamples=nresamples, verboselevel=verboselevel, outgroup=outgroup, branchlengths=bl, dolongbranch=dolongbranch, supportvalues=sv)
    #tc, newicknobranchlens, dt=hpa(fn, resampletype=500000, nresamples=10, verboselevel=2)
    #tc, newicknobranchlens, dt=hpa(fn, resampletype=0, verboselevel=0)
    #v=dt.as_string(schema="newick")
    ## note: dendropy includes newline at end of returned string
    #print v.rstrip()
    if dosv is not None: print "supportvalues:", sv

    if outtree is not None:
        if verboselevel > 0:
            print "saving result tree in newick format to", outtree
        dt.write_to_path(outtree, schema="newick")        # lots of kw options ...
    else:
        print "resulting estimated tree:"
        print dt
