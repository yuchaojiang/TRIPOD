phyper.test=function(peak.gr.g, sig.peak.gr.g, enhancer.db){
  # peak.gr.g: peak granges for gene g
  # sig.peak.gr.g: significant peak granges for gene g
  # enhancer.db: enhancer database granges
  
  # Hypergeomeetric test
  # m+n: Total number of peaks (overlap or dont overlap with enhancers)
  # m: Number peaks that overlap with enhancers
  # k: Number of significant peaks
  # q: Number of significant peaks that overlap with enhancers
  if(length(enhancer.db)==0) {
    pval=NA
  } else{
    m=sum(countOverlaps(peak.gr.g, enhancer.db)>0)
    n=length(peak.gr.g)-m
    k=length(sig.peak.gr.g)
    q=sum(countOverlaps(sig.peak.gr.g, enhancer.db)>0)
    #cat(m,'\t',m+n,'\t',q,'\t',k,'\n')
    if(m==0){ # no peak overlaps with enhancer databases
      pval=NA
    } else{
      pval=phyper(q, m, n, k, lower.tail=F)
    }
  }
  return(pval)
}
