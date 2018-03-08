/* compilation of functions used by sfs_code for operations on binary trees */

/* ------------------------------------------------------------------------- */



/* ------------------------------------------------------------------------- */

void countMuts(struct history *BigHead, long *numMutSeg, long min, long max)
{
  if(BigHead != NULL){
    if(BigHead->event->site >= min && BigHead->event->site <= max)
      (*numMutSeg)++;
    countMuts(BigHead->Rtree, numMutSeg, min, max);
    countMuts(BigHead->Ltree, numMutSeg, min, max);
  }
}

/* ------------------------------------------------------------------------- */

void countSeg(struct history *BigHead, long *numMutSeg)
{
  if(BigHead != NULL){
    (*numMutSeg)++;
    countSeg(BigHead->Rtree, numMutSeg);
    countSeg(BigHead->Ltree, numMutSeg);
  }
}

/* ------------------------------------------------------------------------- */

void treeSize(struct history *BigHead, long *i, int pop) 
{
  if (BigHead != NULL){
    (*i)++;
    if(BigHead->event->numCarriers[pop] == 0 || BigHead->event->numCarriers[pop]
       == gpars.P*ppars[pop].N){
      fprintf(errfile,"odd...  numCarriers = %ld in treeSize\n",
	      BigHead->event->numCarriers[pop]);
      fflush(errfile);
    }
    treeSize(BigHead->Rtree, i, pop);
    treeSize(BigHead->Ltree, i, pop);
  }
}  /* treeSize */

/* ------------------------------------------------------------------------- */

void countSegStorage(struct mutStorage *store, long *numMutSeg)
{
  if(store != NULL){
    (*numMutSeg) += store->numCarry;
    countSegStorage(store->Rtree, numMutSeg);
    countSegStorage(store->Ltree, numMutSeg);
  }
}

/* ------------------------------------------------------------------------- */

void memcpyEVENT(struct event *to, const struct event *from, const int pop)
{
  if(to == NULL)
    exitNOW("NULL event in memcpyEVENT\n");
  to->site = from->site;
  to->gen = from->gen;
  to->genFix[pop] = from->genFix[pop];
  to->genDead[pop] = from->genDead[pop];
  to->ancNuc = from->ancNuc;
  to->derNuc = from->derNuc;
  to->nonsyn = from->nonsyn;
  to->axy = from->axy;
  to->ancAA = from->ancAA;
  to->derAA = from->derAA;
  to->fit = from->fit;
  to->h = from->h;
  to->CpG = from->CpG;
  to->fiveP = from->fiveP;
  to->threeP = from->threeP;
  if(to->nonsyn != '0' && to->nonsyn != '1'){
    to->nSites = from->nSites;
    /* if(1 || to->nonsyn == 'i'){ */
      if(to->nucs == NULL)
	assert(to->nucs = malloc((to->nSites+1)*sizeof(char)));
      else
	assert(to->nucs = realloc(to->nucs, (to->nSites+1)*sizeof(char)));
#ifdef WINDOWS
      strcpy_s(to->nucs, to->nSites+1, from->nucs);
#else
      strcpy(to->nucs, from->nucs);
#endif
    /* } */
  }
}

/* ------------------------------------------------------------------------- */

void checkCarryNum(struct history *tree, long **checkCarry, long numCopies,
		   long *badSite)
{
  if(tree == NULL)
    return;
  (*checkCarry)[tree->event->index] += numCopies;
  if((*badSite)<0 && (*checkCarry)[tree->event->index] > 1){
    (*badSite) = tree->event->index;
    return;
  }
  checkCarryNum(tree->Ltree, checkCarry, numCopies, badSite);
  checkCarryNum(tree->Rtree, checkCarry, numCopies, badSite);
}

/* ------------------------------------------------------------------------- */

void PrintHistTree(struct history *BigHead, long *i, int pop) 
{
  if (BigHead != NULL){
    PrintHistTree(BigHead->Ltree,i, pop);
    (*i)++;
    if(BigHead->Parent != NULL)
      fprintf(errfile,
	      "%3ld: (%lu.%ld) %4lu.%ld %5ld F=%5ld D=%5ld %c %c %c %2i %2i FIT=%f h=%f %c (%ld) %c %c [%c:%d] (%d/%ld/%d)\n",
	      *i, BigHead->Parent->event->site, BigHead->Parent->event->index,
	      BigHead->event->site, BigHead->event->index, BigHead->event->gen,
	      BigHead->event->genFix[pop], BigHead->event->genDead[pop],
	      BigHead->event->ancNuc, BigHead->event->derNuc, 
	      BigHead->event->nonsyn, BigHead->event->ancAA, 
	      BigHead->event->derAA, BigHead->event->fit, BigHead->event->h,
	      BigHead->event->CpG, 
	      BigHead->event->numCarriers[pop],
	      BigHead->event->fiveP, BigHead->event->threeP,
	      BigHead->event->free,
	      popn[pop].polySites[0][BigHead->event->site],
	      BigHead->event->checkPop, BigHead->event->checkGen,
	      BigHead->event->checkStep);
    else{
      fprintf(errfile, "%3ld:(xx) %4lu.%ld %5ld F=%5ld D=%5ld %c %c %c %2i %2i FIT=%f h=%f %c (%ld) %c %c[%c:%d] (%d/%ld/%d)\n",
	      *i, BigHead->event->site, BigHead->event->index,
	      BigHead->event->gen,BigHead->event->genFix[pop],
	      BigHead->event->genDead[pop],
	      BigHead->event->ancNuc,
	      BigHead->event->derNuc, BigHead->event->nonsyn, 
	      BigHead->event->ancAA, BigHead->event->derAA, BigHead->event->fit,
	      BigHead->event->h,
	      BigHead->event->CpG, BigHead->event->numCarriers[pop],
	      BigHead->event->fiveP, BigHead->event->threeP,
	      BigHead->event->free,
	      popn[pop].polySites[0][BigHead->event->site],
	      BigHead->event->checkPop, BigHead->event->checkGen,
	      BigHead->event->checkStep);
      fflush(errfile);
    }
    fflush(errfile);
    PrintHistTree(BigHead->Rtree, i, pop);
  }
  fflush(errfile);
}  /* PrintHistTree */

/* ------------------------------------------------------------------------- */

void PrintErrorStats(const int pop, const long r)
{
  long i, j, k, m;
  for(i=0; i<mutAr[r]->numMuts; i++){
    fprintf(errfile,
	    "%4ld: %4lu %5ld %5ld %c %c %c %2i %2i %1.2f %c %c %c [%c%c] %d",
	    i, mutAr[r]->muts[i]->event->site,
	    mutAr[r]->muts[i]->event->gen,
	    mutAr[r]->muts[i]->event->genFix[pop],
	    mutAr[r]->muts[i]->event->ancNuc,
	    mutAr[r]->muts[i]->event->derNuc,
	    mutAr[r]->muts[i]->event->nonsyn,
	    mutAr[r]->muts[i]->event->ancAA, 
	    mutAr[r]->muts[i]->event->derAA,
	    mutAr[r]->muts[i]->event->fit, mutAr[r]->muts[i]->event->CpG,
	    mutAr[r]->muts[i]->event->fiveP,
	    mutAr[r]->muts[i]->event->threeP,
	    mutAr[r]->muts[i]->event->free,
	    mutAr[r]->muts[i]->event->fixed[pop],
	    popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
    for(k=0; k<gpars.NPOP; k++){
/*       if(mutAr[r]->muts[i]->event->free == '0' && */
/* 	 mutAr[r]->muts[i]->event->fixed[k] == '0'){ */
	fprintf(errfile," (p=%ld; mH=%ld); %ld: ", k,
		mutAr[r]->muts[i]->event->maxHaps[k],
		mutAr[r]->muts[i]->event->numCarriers[k]);
	m = 0;
	for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[k]; j++){
	  if(mutAr[r]->muts[i]->event->hapFreq[k][j] != NULL){
	    fprintf(errfile,
		    "%ld.%ld ",j,*mutAr[r]->muts[i]->event->hapFreq[k][j]);
	    m++;
	  }
	}
/*       } */
    }
    fprintf(errfile, "\n");
    fflush(errfile);
  }
  fprintf(errfile,"\n\n trees\n\n");
  fflush(errfile);
  for(i=0; i<=popn[pop].maxchr; i++){
    fprintf(errfile,"%ld.%ld = %ld: ",i,r,popn[pop].numCopies[i][r]);
    for(j=0; j<gpars.P*ppars[pop].N; j++){
      if(popn[pop].parentLoc[j][r] == i)
	fprintf(errfile,"%ld ",j);
    }
    fprintf(errfile,"\n");
    fflush(errfile);
    j=0;
    PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &j, pop);
    fflush(errfile);
  }
}

/* ------------------------------------------------------------------------- */

void getFitTree(struct history *BigHead, double *fit, const int pop,
		const long locus)
{
  if(ppars[pop].neutpop)
    return;
  if(BigHead != NULL){
    if(gpars.ADDITIVE)
      (*fit) += BigHead->event->fit;
    else
      (*fit) *= 1+BigHead->event->fit;
    if(gpars.ADDITIVE && (*fit) <= -1.0){
      (*fit) = -1.0;
      return;
    }
    else if(!gpars.ADDITIVE && (*fit) <= DBL_MIN){
      (*fit) = 0.0;
      return;
    }
    getFitTree(BigHead->Ltree, fit, pop, locus);
    getFitTree(BigHead->Rtree, fit, pop, locus);
  }
}

/* ------------------------------------------------------------------------- */

void getFitTreeDom(struct history *BigHead, struct history *BigHead2,
		   double *fit, const int pop, const int incHom)
{
  if(ppars[pop].neutpop)
    return;
  if(BigHead != NULL){
    if(fabs(BigHead->event->h - 1.0) <= FLT_MIN){ /* codominant mut */
      if(gpars.ADDITIVE)
	(*fit) += BigHead->event->fit;
      else
	(*fit) *= 1+BigHead->event->fit;
      if(gpars.ADDITIVE && (*fit) <= -1.0){
	(*fit) = -1.0;
	return;
      }
      else if(!gpars.ADDITIVE && (*fit) <= DBL_MIN){
	(*fit) = 0.0;
	return;
      }
    }
    else{ /* dom/recessive mutation */
      int foo=0;
      carryMut(BigHead2->Rtree, BigHead->event->index, &foo);
      if(foo == 0){ /* heterozygous */
	if(gpars.ADDITIVE)
	  (*fit) += BigHead->event->h*BigHead->event->fit;
	else
	  (*fit) *= 1+BigHead->event->h*BigHead->event->fit;
	if(gpars.ADDITIVE && (*fit) <= -1.0){
	  (*fit) = -1.0;
	  return;
	}
	else if(!gpars.ADDITIVE && (*fit) <= DBL_MIN){
	  (*fit) = 0.0;
	  return;
	}
      }
      else if(incHom){ /* homozygous */
	if(gpars.ADDITIVE)
	  (*fit) += BigHead->event->fit;
	else
	  (*fit) *= 1+BigHead->event->fit;
	if(gpars.ADDITIVE && (*fit) <= -1.0){
	  (*fit) = -1.0;
	  return;
	}
	else if(!gpars.ADDITIVE && (*fit) <= DBL_MIN){
	  (*fit) = 0.0;
	  return;
	}
      }
    }
    getFitTreeDom(BigHead->Ltree, BigHead2, fit, pop, incHom);
    getFitTreeDom(BigHead->Rtree, BigHead2, fit, pop, incHom);
  }
}

/* ------------------------------------------------------------------------- */

void carryMut(struct history *tree, long index, int *foo)
{
  if(tree==NULL || *foo>0)
    return;
  else if(tree->event!=NULL && tree->event->index==index){
    (*foo)++;
    return;
  }
  carryMut(tree->Ltree, index, foo);
  carryMut(tree->Rtree, index, foo);
}

/* ------------------------------------------------------------------------- */

void checkHitProbs(const double *hpTMP, const char *err, const int pop,
		   const long chr, const long locus)
{
  unsigned long i;
  long j;
  char c0,c1,c2;
  double dummy = 0;
  
  c0 = retNucHistory(0, &popn[pop].BigHead[chr][locus], locus, pop, 0);
  c1 = retNucHistory(1, &popn[pop].BigHead[chr][locus], locus, pop, 0);
  c2 = retNucHistory(2, &popn[pop].BigHead[chr][locus], locus, pop, 0);
  for(i=1;i<gpars.L[locus]-1;i++){
    for(j=0;j<4;j++)
      if(j != c1-'0')
	dummy += Q[16*(c0-'0')+4*(c1-'0')+(c2-'0')][j];
    if(fabs(dummy - hpTMP[i]) > FLT_MIN){
      fprintf(errfile,"hp NONMATCH (%ld) ind=%ld: ",INITIALSEED, chr);
      fprintf(errfile,"site=%lu (%c%c%c/%c%c%c:%d%d%d); dummy=%f; hp=%lf(%e)\n",
	      i, c0, c1, c2, popn[pop].conSeq[locus][i-1],
	      popn[pop].conSeq[locus][i], popn[pop].conSeq[locus][i+1], 
	      popn[pop].polySites[locus][i-1], popn[pop].polySites[locus][i],
	      popn[pop].polySites[locus][i+1], dummy, hpTMP[i],
	      fabs(dummy-hpTMP[i]));
      fprintf(errfile,"tree:\n");
      fflush(errfile);
      j=0;
      PrintHistTree(popn[pop].BigHead[chr][locus]->Rtree, &j, pop);
      fflush(errfile);
      fflush(errfile);
      exitNOW(err);
    }
    
    if(i < gpars.L[locus]-2){
      c0=c1;
      c1=c2;
      c2 = retNucHistory(i+2, &popn[pop].BigHead[chr][locus], locus, pop, 0);
    }
  }
  if(fabs(dummy - hpTMP[gpars.L[locus]-1]) > FLT_MIN){
    fprintf(errfile,"hp NONMATCH (%ld): site=%lu; dummy=%f; hp=%f (%e)\n",
	    INITIALSEED,gpars.L[locus]-1,dummy,hpTMP[gpars.L[locus]-1],
	    fabs(dummy-hpTMP[gpars.L[locus]-1]));
    j=0;
    PrintHistTree(popn[pop].BigHead[chr][locus]->Rtree, &j, pop);
    fflush(errfile);
    exitNOW(err);
  }
}

/* ------------------------------------------------------------------------- */

void checkForErrors(int pop, char *errMessIn)
{
  int foo2, foo3=0;
  long i, j, k, m, r, *tmpPS, *checkCarry;
  long foo;
  double tfit;
  char errMess[100], cod[3];
  if(ppars[pop].ALIVE==0)
    return;
#ifdef UNIT_TEST
  fprintf(errfile,"checking %s\n",errMessIn);
  fflush(errfile);
#endif
  for(r=0; r<gpars.R; r++){
#ifdef WINDOWS
    sprintf_s(errMess, 100, "error %s!: gen=%ld, pop=%d, r=%ld, iSEED=%ld",
	      errMessIn, popn[pop].gen, pop, r, INITIALSEED);
#else
    sprintf(errMess, "error %s!: gen=%ld, pop=%d, r=%ld, iSEED=%ld",
	    errMessIn, popn[pop].gen, pop, r, INITIALSEED);
#endif

    /* check for stop codons in sequence */
    if(gpars.ANNOTATE[r] == 'C'){
      for(i=0; i<gpars.L[r]; i+=3){
	cod[0] = popn[pop].conSeq[r][i];
	cod[1] = popn[pop].conSeq[r][i+1];
	cod[2] = popn[pop].conSeq[r][i+2];
	if(AA(cod) == 0){
	  fprintf(errfile,"found stop codon at position %ld %c%c%c\n",
		 i, cod[0], cod[1], cod[2]);
	  abort();
	}
      }
    }

    foo = 0;
    for(i=0; i<=popn[pop].maxchr; i++){
      foo += popn[pop].numCopies[i][r];
    }
    if(foo != gpars.P*ppars[pop].N){
      fprintf(errfile,"error!!  numCopies not right!  E=%ld O=%ld\n",
	     gpars.P*ppars[pop].N, foo);
      abort();
    }
    for(i=0; i<=popn[pop].maxchr; i++){
      foo = popn[pop].numCopies[i][r];
      for(j=0; j<gpars.P*ppars[pop].N; j++){
	if(popn[pop].parentLoc[j][r] == i)
	  foo--;
      }
      if(foo != 0){
	fprintf(errfile, "error %s: numCopies not right, diff=%ld: PL=%ld\n",
		errMess, foo, i);
	abort();
      }
    }

#ifdef VERBOSE_DEBUG
    fprintf(errfile,"numMuts[%ld] = %ld (next=%ld)\n", r, mutAr[r]->numMuts,
	   mutAr[r]->mutIndex);
    fflush(errfile);
#endif
    for(i=0; i<mutAr[r]->numMuts; i++){
#ifdef VERBOSE_DEBUG
      fprintf(errfile,"pop=%d; i=%ld; ",pop, i);
      fflush(errfile);
      fprintf(errfile,"site=%ld, ",mutAr[r]->muts[i]->event->site);
      fflush(errfile);
      fprintf(errfile,"r=%ld; ", r);
      fflush(errfile);
      int foo2;
      fprintf(errfile,"polysites=");
      for(foo2=0; foo2<gpars.NPOP; foo2++){
	if(ppars[foo2].ALIVE)
	  fprintf(errfile,"%d;", popn[foo2].polySites[r][mutAr[r]->muts[i]->event->site]);
      }
      fprintf(errfile," maxHaps=%ld; ", mutAr[r]->muts[i]->event->maxHaps[pop]);
      fflush(errfile);
#endif
      mutAr[r]->muts[i]->event->numCarriers[pop] = 0;
      for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
	if(mutAr[r]->muts[i]->event->hapFreq[pop][j] != NULL)
	  mutAr[r]->muts[i]->event->numCarriers[pop] +=
	    (*mutAr[r]->muts[i]->event->hapFreq[pop][j]);
      }

#ifdef VERBOSE_DEBUG
      fprintf(errfile,"numCarry[%ld]=%ld; ", mutAr[r]->muts[i]->event->index, 
	     mutAr[r]->muts[i]->event->numCarriers[pop]);
      fflush(errfile);
      fprintf(errfile,"fixed=%c\n",mutAr[r]->muts[i]->event->fixed[pop]);
      fflush(errfile);
#endif
      
      for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
	if(mutAr[r]->muts[i]->event->hapFreq[pop][j] == NULL)
	  continue;
	foo=0;
	for(k=0; k<=popn[pop].maxchr; k++){
	  if(mutAr[r]->muts[i]->event->hapFreq[pop][j] == 
	     &popn[pop].numCopies[k][r]){
	    foo=1;
	    break;
	  }
	}
	if(foo==0){
	  fprintf(errfile,"error: mut %ld (%ld) in hapFreq=%ld, not in pop %d!\n",
		 i, mutAr[r]->muts[i]->event->site, 
		 *mutAr[r]->muts[i]->event->hapFreq[pop][j], pop);
	  abort();
	}
      }
      
      foo3=0;
      foo=0;
      for(j=0; j<=popn[pop].maxchr; j++){
	if(popn[pop].numCopies>0){
	  foo2=0;
	  carryMut(popn[pop].BigHead[j][r], mutAr[r]->muts[i]->event->index,
		   &foo2);
	  if(foo2 == 1){
	    for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[pop]; k++){
	      if(mutAr[r]->muts[i]->event->hapFreq[pop][k] ==
		 &(popn[pop].numCopies[j][r])){
		foo=1;
	      }
	    }
	    if(foo==0){
	      fprintf(errfile,"chrom %ld (%ld) carries mut %ld, but not included \
in hapFreq list!\n", j, popn[pop].numCopies[j][r],
		     mutAr[r]->muts[i]->event->index);
	      fflush(errfile);
	      foo3=1;
	    }
	  }
	}
      }

      foo = 0;
      for(j=0; j<=popn[pop].maxchr; j++){
	foo2 = 0;
	carryMut(popn[pop].BigHead[j][r]->Rtree, i, &foo2);
	if(foo2 > 1){
	  fprintf(errfile,"error: %s! ind %ld carries mutation %ld multiple times!\n", 
		 errMess, j, i);
	  fflush(errfile);
	  i=0;
	  PrintHistTree(popn[pop].BigHead[j][r]->Rtree, &i, pop);
	  fflush(errfile);
	  abort();
	}
	foo += foo2;
	if(foo == 1 && (mutAr[r]->muts[i]->event->numCarriers[pop] == 0 ||
			mutAr[r]->muts[i]->event->fixed[pop] == '1' ||
			mutAr[r]->muts[i]->event->free == '1')){
	  fprintf(errfile,
		  "error %s:mut %ld (%ld) nc=0, but BH[%d][%ld][%ld] (%ld)!!\n",
		  errMess, i, mutAr[r]->muts[i]->event->site, pop, j, r, 
		  popn[pop].numCopies[j][r]);
	  fprintf(errfile, "%ld\t%c\t%c\n",mutAr[r]->muts[i]->event->numCarriers[pop], mutAr[r]->muts[i]->event->fixed[pop], mutAr[r]->muts[i]->event->free);
	  fprintf(errfile,"\n\n");
	  PrintErrorStats(pop, r);
	  foo=0;
	  PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &foo, pop);
	  abort();
	}
      }
      if(foo == 0 && mutAr[r]->muts[i]->event->numCarriers[pop] > 0){
	fprintf(errfile,
		"error: mutation %ld (%ld), nc>0 but has no carriers!\n",i,
		mutAr[r]->muts[i]->event->site);
	PrintErrorStats(pop, r);
	fflush(errfile);
	abort();
      }
    }
    m = 0;
    for(i=0; i<=popn[pop].maxchr; i++){
      m += popn[pop].numCopies[i][r];
      if(popn[pop].numCopies[i][r] < 0){
	for(k=0; k<gpars.P*ppars[pop].N; k++){
	  if(popn[pop].parentLoc[k][r] == i){
	    fprintf(errfile,"error %s (%ld): parentLoc[%d][%ld][%ld]= %ld\n",
		    errMess, INITIALSEED,pop,k,r,i);
	    fprintf(errfile,"But numCopies[%d][%ld][%ld] = %ld!!\n",
		    pop,i,r,popn[pop].numCopies[i][r]);
	    fflush(errfile);
	    abort();			    
	  }
	}
	fprintf(errfile,"error %s(%ld): popn[%d].numCopies[%ld][%ld] = %ld!!\n",
		errMess, INITIALSEED, pop, i, r, popn[pop].numCopies[i][r]);
	abort();
      }
      else if(popn[pop].numCopies[i][r] == 0)
	continue;

      if(gpars.substMod == 5 && strcmp(errMessIn, "After regeneratePOP") != 0)
	checkHitProbs(popn[pop].hpSUM[i][r], errMess, pop, i, r);
      tfit = (gpars.ADDITIVE ? 0.0 : 1.0);
      getFitTree(popn[pop].BigHead[i][r]->Rtree, &tfit, pop, r);
      if(fabs(tfit-popn[pop].indfit[i][r])>1e-12){
	fprintf(errfile, "%s: fit=%.12e; indfit[%d][%ld][%ld] = %.12e; delta=%.9e; err_tol=%.9e\n",errMess,
		tfit, pop, i, r, popn[pop].indfit[i][r], fabs(tfit-popn[pop].indfit[i][r]), 10*FLT_MIN);
	foo=0;
	PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &foo, pop);
	abort();
      }
      checkHistTreeOrder(popn[pop].BigHead[i][r]->Rtree, 0, gpars.L[r], pop);
    }
    if(m != gpars.P*ppars[pop].N){
      fprintf(errfile,"error %s (%ld): numCopies does not sum to \
popsize!  sum = %ld; popsize = %ld for locus %ld\n",errMess, INITIALSEED, m, 
	      gpars.P*ppars[pop].N, r);
      abort();
    }
    
    assert(tmpPS = malloc(gpars.L[r]*sizeof(*tmpPS)));
    for(i=0; i<gpars.L[r]; i++)
      tmpPS[i] = 0;
    for(i=0; i<mutAr[r]->numMuts; i++){
      if(mutAr[r]->muts[i]->event->free == '0' &&
	 mutAr[r]->muts[i]->event->fixed[pop] == '0' && 
	 mutAr[r]->muts[i]->event->numCarriers[pop] > 0 && 
	 (mutAr[r]->muts[i]->event->nonsyn == '0' ||
	  mutAr[r]->muts[i]->event->nonsyn == '1')){
	foo = 0;
	for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
	  if(mutAr[r]->muts[i]->event->hapFreq[pop][j] != NULL && 
	     (*mutAr[r]->muts[i]->event->hapFreq[pop][j]) > 0){
	    foo = 1;
	    break;
	  }
	}
	if(foo == 1)
	  tmpPS[mutAr[r]->muts[i]->event->site]++;
      }
    }
    for(i=0; i<gpars.L[r]; i++){
      if(tmpPS[i] != popn[pop].polySites[r][i]){
	fprintf(errfile,
		"error %s! PS not right, site %ld in locus %ld E=%ld, O=%d\n",
		errMess, i, r, tmpPS[i], popn[pop].polySites[r][i]);
	for(i=0; i<mutAr[r]->numMuts; i++){
	  fprintf(errfile,
		  "%4ld: %4lu %5ld %5ld %c %c %c %2i %2i %1.2f %c %c %c [%c%c]",
		  i, mutAr[r]->muts[i]->event->site,
		  mutAr[r]->muts[i]->event->gen,
		  mutAr[r]->muts[i]->event->genFix[pop],
		  mutAr[r]->muts[i]->event->ancNuc,
		  mutAr[r]->muts[i]->event->derNuc,
		  mutAr[r]->muts[i]->event->nonsyn,
		  mutAr[r]->muts[i]->event->ancAA, 
		  mutAr[r]->muts[i]->event->derAA,
		  mutAr[r]->muts[i]->event->fit, mutAr[r]->muts[i]->event->CpG,
		  mutAr[r]->muts[i]->event->fiveP,
		  mutAr[r]->muts[i]->event->threeP,
		  mutAr[r]->muts[i]->event->free,
		  mutAr[r]->muts[i]->event->fixed[pop]);
	  for(k=0; k<gpars.NPOP; k++){
	    if(ppars[k].ALIVE == 0)
	      continue;
	    if(mutAr[r]->muts[i]->event->free == '0' &&
	       mutAr[r]->muts[i]->event->fixed[k] == '0'){
	      fprintf(errfile," (p=%ld; mH=%ld); %ld: ", k,
		      mutAr[r]->muts[i]->event->maxHaps[k],
		      mutAr[r]->muts[i]->event->numCarriers[k]);
	      m = 0;
	      for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[k]; j++){
		if(mutAr[r]->muts[i]->event->hapFreq[k][j] != NULL){
		  fprintf(errfile,
			  "%ld=%ld ",j,*mutAr[r]->muts[i]->event->hapFreq[k][j]);
		  m++;
		}
	      }
	    }
	  }
	  fprintf(errfile, "\n");
	  fflush(errfile);
	}
	fprintf(errfile,"\n\n trees\n\n");
	fflush(errfile);
	for(i=0; i<popn[pop].maxchr; i++){
	  fprintf(errfile,"%ld.%ld = %ld:\n",i,r,popn[pop].numCopies[i][r]);
	  j=0;
	  PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &j, pop);
	  fflush(errfile);
	}
	abort();
      }
    }
    free(tmpPS);
    
    assert(checkCarry = malloc((mutAr[r]->numMuts)*sizeof(*checkCarry)));
    for(i=0; i<mutAr[r]->numMuts; i++)
      checkCarry[i] = 0;
    for(i=0; i<gpars.P*ppars[pop].N; i++){
      for(j=0; j<mutAr[r]->numMuts; j++)
	checkCarry[j] = 0;
      k = -1;
      checkCarryNum(popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
		    &checkCarry, 1, &k);
      if(k>=0){
	fprintf(errfile, "site %ld is carried multiple times by chrom %ld (%ld)\
; tree:\n", k, i, popn[pop].parentLoc[i][r]);
	j=0;
	PrintHistTree(popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
		      &j, pop);
	abort();
      }
    }
    for(i=0; i<mutAr[r]->numMuts; i++)
      checkCarry[i] = 0;
    for(i=0; i<gpars.P*ppars[pop].N; i++){
      j = 0;
      checkCarryNum(popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
		    &checkCarry, 1, &j);
    }
    for(i=0; i<mutAr[r]->numMuts; i++){
      if(mutAr[r]->muts[i]->event->free == '1')
	continue;
/* #ifdef VERBOSE_DEBUG */
/*       fprintf(errfile,"%ld::%ld (site=%ld)\n",i,checkCarry[i], */
/* 	     mutAr[r]->muts[i]->event->site); */
/*       fflush(errfile); */
/* #endif */
      if(mutAr[r]->muts[i]->event->free == '0' &&
	 mutAr[r]->muts[i]->event->fixed[pop] == '0'){
	if(mutAr[r]->muts[i]->event->numCarriers[pop] != checkCarry[i]){
	  fprintf(errfile,"not the right number of carriers at mut %ld(%ld)! %ldvs%ld\n",
		 i,mutAr[r]->muts[i]->event->site,
		 mutAr[r]->muts[i]->event->numCarriers[pop], checkCarry[i]);
	  fflush(errfile);
	  foo3=1;
	  for(j=0; j<gpars.P*ppars[pop].N; j++){
	    foo2=0;
	    carryMut(popn[pop].BigHead[popn[pop].parentLoc[j][r]][r]->Rtree,
		     mutAr[r]->muts[i]->event->index, &foo2);
	    if(foo2==0)
	      continue;
	    foo=0;
	    for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[pop]; k++){
	      if(mutAr[r]->muts[i]->event->hapFreq[pop][k] ==
		 &popn[pop].numCopies[popn[pop].parentLoc[j][r]][r]){
		foo=1;
		break;
	      }
	    }
	    if(foo==0){
	      fprintf(errfile,"chromosome %ld (pl=%ld) not in hapFreq!\n",
		     j, popn[pop].parentLoc[j][r]);
	      fflush(errfile);
	    }
	  }
	}
	for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
	  if(mutAr[r]->muts[i]->event->hapFreq[pop][j] != NULL)
	    checkCarry[i] -= (*mutAr[r]->muts[i]->event->hapFreq[pop][j]);
	}
	if(checkCarry[i] != 0){
	  fprintf(errfile, "error %s: carriers of mutation %ld (site=%ld) \
locus %ld (diff=%ld)\n", errMess, i, mutAr[r]->muts[i]->event->site, r,
		  checkCarry[i]);
	  fflush(errfile);
	  foo3=1;
	  foo=0;
	  for(j=0; j<=popn[pop].maxchr; j++){
	    if(popn[pop].numCopies>0){
	      foo2=0;
	      carryMut(popn[pop].BigHead[j][r], mutAr[r]->muts[i]->event->index,
		       &foo2);
	      if(foo2 == 1){
		for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[pop]; k++){
		  if(mutAr[r]->muts[i]->event->hapFreq[pop][k] ==
		     &(popn[pop].numCopies[j][r])){
		    foo=1;
		  }
		}
		if(foo==0){
		  fprintf(errfile,"chrom %ld (%ld) carries mut %ld, but not included \
in hapFreq list!\n", j, popn[pop].numCopies[j][r],
			 mutAr[r]->muts[i]->event->index);
		  fflush(errfile);
		  abort();
		}
	      }
	    }
	  }
	  
	  fflush(errfile);
	  PrintErrorStats(pop, r);

	  foo = 0;
	  for(k=0; k<popn[pop].maxchr; k++){
	    if(popn[pop].numCopies[k][r] < 0){
	      fprintf(errfile,"numCopies[%ld] == %ld!!!\n",k,popn[pop].numCopies[k][r]);
	      fflush(errfile);
	    }
	    if(popn[pop].numCopies[k][r] == 0)
	      continue;
	    foo += popn[pop].numCopies[k][r];
	    foo2=0;
	    carryMut(popn[pop].BigHead[k][r]->Rtree,
		     mutAr[r]->muts[i]->event->index, &foo2);
	    if(foo2 == 0){
	      fprintf(errfile,"chr %ld does not carry mutation (%ld)!!!\n", k,
		     popn[pop].numCopies[k][r]);
	      fflush(errfile);
	      for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
		if(mutAr[r]->muts[i]->event->hapFreq[pop][j] ==
		   &popn[pop].numCopies[k][r]){
		  fprintf(errfile,"mutation linked to chr %ld!!\n",k);
		  fflush(errfile);
		}
	      }
	    }
	  }
	  fprintf(errfile,"total chrms:  %ld\n",foo);
	  fflush(errfile);
	  foo2 = 0;
	  for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
	    if(mutAr[r]->muts[i]->event->hapFreq[pop][j] != NULL){
	      foo -= *mutAr[r]->muts[i]->event->hapFreq[pop][j];
	      for(k=0; k<=popn[pop].maxchr; k++){
		if(mutAr[r]->muts[i]->event->hapFreq[pop][j] ==
		   &popn[pop].numCopies[k][r]){
		  fprintf(errfile,"hap[%ld] == chr %ld (%ld)\n",j, k,
			 popn[pop].numCopies[k][r]);
		  foo2 += popn[pop].numCopies[k][r];
		  break;
		}
	      }
	    }
	  }
	  fprintf(errfile,"sum = %d\n", foo2);
	  fprintf(errfile,"remainder = %ld\n",foo);
	  fflush(errfile);
	}
      }
    }
    if(foo3==1)
      abort();
    free(checkCarry);
  }
  fprintf(errfile,"passed checks %s\n", errMessIn);
}

/* ------------------------------------------------------------------------- */

void newDead(long chr, long loc)
{
  struct deadList *tmp=NULL;
  assert(tmp = malloc(sizeof(*tmp)));
  tmp->prev = NULL;
  if(dead->next == NULL)
    tmp->next = NULL;
  else{
    tmp->next = dead->next;
    dead->next->prev = tmp;
  }
  dead->next = tmp;
  tmp->chr = chr;
  tmp->loc = loc;
}

/* ------------------------------------------------------------------------- */

struct history* popTrash()
{
  struct history *foo;
  if(trash->Rtree == NULL){
    assert(foo = malloc(sizeof(*foo)));
  }
  else if(trash->Rtree->Ltree != NULL){
    foo = trash->Rtree;
    trash->Rtree = foo->Ltree;
    if(foo->Rtree != NULL){
      trash->Rtree->Parent = foo->Rtree;
      foo->Rtree->Parent = foo->Parent;
    }
    else
      trash->Rtree->Parent = foo->Parent;
  }
  else if(trash->Rtree->Rtree != NULL){
    foo = trash->Rtree;
    trash->Rtree = foo->Rtree;
    trash->Rtree->Parent = foo->Parent;
  }
  else{ /* both subtrees are NULL */
    foo = trash->Rtree;
    if(foo->Parent != NULL)
      trash->Rtree = foo->Parent;
    else
      trash->Rtree = NULL;
  }
  foo->event = NULL;
  foo->Rtree = NULL;
  foo->Ltree = NULL;
  foo->Parent = NULL;
  return foo;
}

/* ------------------------------------------------------------------------- */

void copyHistory(struct history **to, struct history **from, const long locus,
		 const int popTo, const long chrTo, const int popFrom, 
		 const int updatePS, const int fixed, const int cpNC)
{
  int k;
  long i;
  if((*from) == NULL)
    return;
  if((*from)->Rtree == NULL)
    (*to)->Rtree = NULL;
  else{
    (*to)->Rtree = popTrash();
    (*to)->Rtree->event = (*from)->Rtree->event;
    if(cpNC)
      (*to)->Rtree->event->numCarriers[popTo] = 
	(*to)->Rtree->event->numCarriers[popFrom];
    else
      (*to)->Rtree->event->numCarriers[popTo] +=
	popn[popTo].numCopies[chrTo][locus];

    if(chrTo != -1){ /* not a substitution, set hapFreq */
      k = 0;
      for(i=(*to)->Rtree->event->nextHap[popTo]; 
	  i<(*to)->Rtree->event->maxHaps[popTo]; i++){
	if((*to)->Rtree->event->hapFreq[popTo][i] == NULL)
	  break;
	else if(updatePS == 1 && (*(*to)->Rtree->event->hapFreq[popTo][i]) > 0)
	  k = 1;
      }
      (*to)->Rtree->event->nextHap[popTo] = i;
      
      if(updatePS == 1 && k == 0){  /* adding mutation to new pop */
	for(i=0; i<(*to)->Rtree->event->maxHaps[popTo]; i++){
	  if((*to)->Rtree->event->hapFreq[popTo][i] != NULL &&
	     (*(*to)->Rtree->event->hapFreq[popTo][i]) > 0){
	    k = 1;  /* mutation already segregating in population */
	    break;
	  }
	}
	if(k == 0){
	  if((*to)->Rtree->event->nonsyn == '0' ||
	     (*to)->Rtree->event->nonsyn == '1'){
	    popn[popTo].polySites[locus][(*to)->Rtree->event->site]++;
	    (*to)->Rtree->event->genDead[popTo] = -LONG_MAX;
#ifdef VERBOSE_DEBUG
	    fprintf(errfile,"update PS[%ld]=%d\n",(*to)->Rtree->event->site, 
		   popn[popTo].polySites[locus][(*to)->Rtree->event->site]);
	    if(popn[popTo].polySites[locus][(*to)->Rtree->event->site]>1){
	      //getchar();
	    }
#endif
	  }
	}
      }
      if((*to)->Rtree->event->nextHap[popTo] ==
	 (*to)->Rtree->event->maxHaps[popTo]){ /* allocate new haplotype */
	(*to)->Rtree->event->maxHaps[popTo]++;
	assert((*to)->Rtree->event->hapFreq[popTo] = 
	       realloc((*to)->Rtree->event->hapFreq[popTo],
		       (*to)->Rtree->event->maxHaps[popTo]*sizeof(long*)));
      }
      (*to)->Rtree->event->hapFreq[popTo][(*to)->Rtree->event->nextHap[popTo]] =
	&popn[popTo].numCopies[chrTo][locus];
    }
    (*to)->Rtree->Rtree = NULL;
    (*to)->Rtree->Ltree = NULL;
    if((*from)->event == NULL) /* root of tree */
      (*to)->Rtree->Parent = NULL;
    else
      (*to)->Rtree->Parent = (*to);
    
    if(fixed == 1)  (*to)->Rtree->event->fixed[popTo] = '1';
    copyHistory(&(*to)->Rtree, &(*from)->Rtree, locus, popTo, chrTo, popFrom,
		updatePS, fixed, cpNC);
  }
  
  if((*from)->Ltree == NULL)
    (*to)->Ltree = NULL;
  else{
    (*to)->Ltree = popTrash();
    (*to)->Ltree->event = (*from)->Ltree->event;
    if(cpNC)
      (*to)->Ltree->event->numCarriers[popTo] = 
	(*to)->Ltree->event->numCarriers[popFrom];
    else
      (*to)->Ltree->event->numCarriers[popTo] += 
	popn[popTo].numCopies[chrTo][locus];

    if(chrTo != -1){ /* not a substitution, set hapFreq */
      k = 0;
      for(i=(*to)->Ltree->event->nextHap[popTo]; 
	  i<(*to)->Ltree->event->maxHaps[popTo]; i++){
	if((*to)->Ltree->event->hapFreq[popTo][i] == NULL)
	  break;
	else if(updatePS == 1 && (*(*to)->Ltree->event->hapFreq[popTo][i]) > 0)
	  k = 1;
      }
      (*to)->Ltree->event->nextHap[popTo] = i;
      if(updatePS == 1 && k == 0){  /* adding mutation to new pop */
	for(i=0; i<(*to)->Ltree->event->maxHaps[popTo]; i++){
	  if((*to)->Ltree->event->hapFreq[popTo][i] != NULL &&
	     (*(*to)->Ltree->event->hapFreq[popTo][i]) > 0){
	    k = 1;
	    break;
	  }
	}
	if(k == 0){
	  if((*to)->Ltree->event->nonsyn == '0' || 
	     (*to)->Ltree->event->nonsyn == '1'){
	    popn[popTo].polySites[locus][(*to)->Ltree->event->site]++;
	    (*to)->Ltree->event->genDead[popTo] = -LONG_MAX;

#ifdef VERBOSE_DEBUG
	    fprintf(errfile,"update PS[%ld]=%d\n",(*to)->Ltree->event->site, 
		   popn[popTo].polySites[locus][(*to)->Ltree->event->site]);

	    if(popn[popTo].polySites[locus][(*to)->Ltree->event->site]>1){
	      //getchar();
	    }
#endif
	  }
	}
      }
      if((*to)->Ltree->event->nextHap[popTo] ==
	 (*to)->Ltree->event->maxHaps[popTo]){ /* allocate new haplotype */
	(*to)->Ltree->event->maxHaps[popTo]++;
	assert((*to)->Ltree->event->hapFreq[popTo] = 
	       realloc((*to)->Ltree->event->hapFreq[popTo],
		       (*to)->Ltree->event->maxHaps[popTo]*sizeof(long*)));
      }
      (*to)->Ltree->event->hapFreq[popTo][(*to)->Ltree->event->nextHap[popTo]] =
	&popn[popTo].numCopies[chrTo][locus];
    }
    (*to)->Ltree->Rtree = NULL;
    (*to)->Ltree->Ltree = NULL;
    if((*from)->event == NULL) /* root of tree */
      (*to)->Ltree->Parent = NULL;
    else
      (*to)->Ltree->Parent = (*to);
    if(fixed == 1)  (*to)->Ltree->event->fixed[popTo] = '1';
    copyHistory(&(*to)->Ltree, &(*from)->Ltree, locus, popTo, chrTo, popFrom,
		updatePS, fixed, cpNC);
  }
}

/* ------------------------------------------------------------------------- */

void copyHistoryNode(struct history **to, struct history **from,
		     struct history **MOM, const int pop, const long chr,
		     const long locus)
{
  long i;
  if((*to) == NULL){ /* reached leaf node */
    (*to) = popTrash();
    (*to)->event = (*from)->event;
    for(i=(*to)->event->nextHap[pop]; i<(*to)->event->maxHaps[pop]; i++){
      if((*to)->event->hapFreq[pop][i] == NULL)
	break;
    }
    if(i == (*to)->event->maxHaps[pop]){ /* allocate a new haplotype */
      (*to)->event->maxHaps[pop]++;
      assert((*to)->event->hapFreq[pop] = 
	     realloc((*to)->event->hapFreq[pop], (i+1)*sizeof(long*)));
      (*to)->event->nextHap[pop] = 0;
    }
    else
      (*to)->event->nextHap[pop] = i;
    (*to)->event->hapFreq[pop][i] = &popn[pop].numCopies[chr][locus];
    (*to)->event->numCarriers[pop] += popn[pop].numCopies[chr][locus];
    (*to)->Rtree = NULL;
    (*to)->Ltree = NULL;
    if(MOM != NULL)
      (*to)->Parent = *MOM;
    else
      (*to)->Parent = NULL;
    splayHistory(*to, &popn[pop].BigHead[chr][locus], pop);
    return;
  }
  else{
    if((*to)->event->site == (*from)->event->site){
      if((*to)->event->gen < (*from)->event->gen) /* new data from later gen */ 
	copyHistoryNode(&(*to)->Rtree, from, to, pop, chr, locus);
      else if((*to)->event->gen > (*from)->event->gen) /* event from prev gen */
	copyHistoryNode(&(*to)->Ltree, from, to, pop, chr, locus);
      /* else event already exists! */
    }
    else if((*to)->event->site < (*from)->event->site)
      copyHistoryNode(&(*to)->Rtree, from, to, pop, chr, locus);
    else
      copyHistoryNode(&(*to)->Ltree, from, to, pop, chr, locus);
  }
}

/* ------------------------------------------------------------------------- */

void copyPartialHistory(struct history **to, struct history **from,
			const long min, const long max, const int pop,
			const long chr, const long locus)
{
  if((*from) == NULL)
    return;
  else{
    /* add in post-order */
    if((*from)->event->site >= min && (*from)->Ltree != NULL)
      copyPartialHistory(to, &(*from)->Ltree, min, max, pop, chr, locus);
    if((*from)->event->site <= max && (*from)->Rtree != NULL)
      copyPartialHistory(to, &(*from)->Rtree, min, max, pop, chr, locus);
    
    if((*from)->event->site >= min && (*from)->event->site <= max)
      copyHistoryNode(to, from, NULL, pop, chr, locus);
  }
}

/* ------------------------------------------------------------------------- */

void convertTract(struct history **to, struct history **from, const long min,
		  const long max, const int pop, const long chr,
		  const long locus, const long altPar)
{
  int foo;
  if((*from) == NULL)
    return;
  /* add in post-order */
  if((*from)->event->site >= min && (*from)->Ltree != NULL)
    convertTract(to, &(*from)->Ltree, min, max, pop, chr, locus, altPar);
  if((*from)->event->site <= max && (*from)->Rtree != NULL)
    convertTract(to, &(*from)->Rtree, min, max, pop, chr, locus, altPar);
  
  if((*from)->event->site >= min && (*from)->event->site <= max){
#ifdef VERBOSE_DEBUG
    fprintf(errfile,"looking at site %ld(%ld)\n",(*from)->event->site, 
	   (*from)->event->index);
    fflush(errfile);
#endif
    if(fabs(ppars[pop].BGC - 0.5) <= FLT_MIN || /* no BGC, copy directly */
       (((*from)->event->ancNuc-'0') < 2 && ((*from)->event->derNuc-'0') < 2) ||
       (((*from)->event->ancNuc-'0') > 1 && ((*from)->event->derNuc-'0') > 1)){
      /* or C<->G or T<->A mutations not affected by BGC, copy! */
      copyHistoryNode(to, from, NULL, pop, chr, locus);
#ifdef VERBOSE_DEBUG
      fprintf(errfile,"copied site %ld(%ld)\n",(*from)->event->site, 
	     (*from)->event->index);
      fflush(errfile);
#endif
    }
    else{/* GC<->AT need to see if other parent chromosome carries mutation */
      foo = 0;
      carryMut(popn[pop].BigHead[altPar][locus]->Rtree, (*from)->event->index,
	       &foo);
      if(foo == 1 || 
	 (((*from)->event->derNuc-'0') < 2 &&
	  ran1(&gpars.seed) < ppars[pop].BGC) ||
	 (((*from)->event->derNuc-'0') > 1 && 
	  ran1(&gpars.seed)>ppars[pop].BGC)){
	/* mutation carried by both chromosomal parents, AT->GC mutation 
	   converted, or a lucky escape from BGC!!  copy! */
	copyHistoryNode(to, from, NULL, pop, chr, locus);
#ifdef VERBOSE_DEBUG
	fprintf(errfile,"copied site %ld(%ld)\n",(*from)->event->site, 
	       (*from)->event->index);
	fflush(errfile);
#endif
      }
#ifdef VERBOSE_DEBUG
      else{
	fprintf(errfile,"did NOT copy site %ld(%ld)\n",(*from)->event->site, 
	       (*from)->event->index);
	fflush(errfile);
      }
#endif
    }
  }
}

/* ------------------------------------------------------------------------- */

void addToTract(struct history **to, struct history **from, const long min,
		const long max, const int pop, const long chr,
		const long locus, const long altPar)
{
  int foo;
  if((*from) == NULL)
    return;
  if((*from)->event->site >= min && (*from)->Ltree != NULL)
    addToTract(to, &(*from)->Ltree, min, max, pop, chr, locus, altPar);
  if((*from)->event->site <= max && (*from)->Rtree != NULL)
    addToTract(to, &(*from)->Rtree, min, max, pop, chr, locus, altPar);
  if((*from)->event->site >= min && (*from)->event->site <= max){
    if((((*from)->event->ancNuc-'0') < 2 && ((*from)->event->derNuc-'0') < 2) ||
       (((*from)->event->ancNuc-'0') > 1 && ((*from)->event->derNuc-'0') > 1)){
      /* or C<->G or T<->A mutations not affected by BGC, skip! */
      return;
    }
    else{/* GC<->AT need to see if other parent chromosome carries mutation */
      foo = 0;
      carryMut(popn[pop].BigHead[altPar][locus]->Rtree, (*from)->event->index,
	       &foo);
      if(foo == 0 &&  
	 ((((*from)->event->derNuc-'0')<2&&ran1(&gpars.seed)<ppars[pop].BGC) ||
	  (((*from)->event->derNuc-'0')>1&&ran1(&gpars.seed)>ppars[pop].BGC))){
	/* AT->GC mutation only carried by this parent and it is converted
	   or GC->AT mutation that escapes BGC */
	copyHistoryNode(to, from, NULL, pop, chr, locus);
#ifdef VERBOSE_DEBUG
	fprintf(errfile,"copied site %ld(%ld)\n",(*from)->event->site, 
	       (*from)->event->index);
	fflush(errfile);
#endif
      }
#ifdef VERBOSE_DEBUG
      else{
	fprintf(errfile,"did NOT copy site %ld(%ld) (addToTract)\n",(*from)->event->site,
	       (*from)->event->index);
	fflush(errfile);
      }
#endif
    }
  }
}

/* ------------------------------------------------------------------------- */

void getMutSite(unsigned long *site, double *hpSUM, const long locus)
{
  double rn;
  double foo;
  int steps = 0;

  /* use an inverse-CDF method to get site */
  rn=ran1(&gpars.seed);
  
  (*site) = (long)(rn*(gpars.L[locus]-2)+1);   /* initial guess in [1,L-2] */

#ifdef VERBOSE_DEBUG
  fprintf(errfile,"rn = %f; hp = %f;  ",rn,hpSUM[gpars.L[locus]-1]);
  fflush(errfile);
#endif

  rn *= hpSUM[gpars.L[locus]-1];

#ifdef VERBOSE_DEBUG
  fprintf(errfile,"rn = %f;  ",rn);
  fflush(errfile);
  if(rn <= DBL_MIN){
    fprintf(errfile,"rn(2) == zero!\n");
    assert(1);
  }
#endif
  while(1){
    steps++;
    if(*site >= gpars.L[locus]-1)
      *site = gpars.L[locus]-2;
    else if(*site <= 0)
      *site = 0;
    
    if(hpSUM[(*site)] < rn && hpSUM[(*site)+1]>=rn)
      break;
    /* take large steps for the first 5 rounds, then take small steps... */
    if(steps<5 && fabs(foo=(rn-hpSUM[(*site)])) > 3.0){
      foo /= 3.0;
      (*site) += (long)foo;
    }
    else{
      if(rn>hpSUM[(*site)]){
	(*site)++;
      }
      else{
	(*site)--;
      }
    }
  }
  (*site)++;
  
#ifdef VERBOSE_DEBUG
  fprintf(errfile,"%f - %f (steps = %d, site=%ld)\n",hpSUM[(*site)-1],hpSUM[(*site)],steps,*site);
  fflush(errfile);
#endif

  if((*site)>=gpars.L[locus]-1 || *site <= 0){
    fprintf(errfile,"damn site(%ld):%lu; rn=%f; hpSUM[%lu]=%f; hpSUM[%lu]=%f\n",
	    INITIALSEED,*site,rn,gpars.L[locus]-2,hpSUM[gpars.L[locus]-2],
	    gpars.L[locus]-1,hpSUM[gpars.L[locus]-1]);
    abort();
  }
}

/* ------------------------------------------------------------------------- */

long invCDF(float rn, double *dist, long beg, long end)
{
  long steps = 0, MAXSTEPS=end-beg+1;
  long site;
  long numBigSteps;
  long toobig=end, toosmall=beg;
  double avgStep = (dist[end]-dist[beg])/MAXSTEPS;

  if(rn <= dist[beg]+DBL_MIN)
    return beg;
  else if(rn >= dist[end]-DBL_MIN)
    return end;
  else if((beg > 0 && rn < dist[beg-1]) || rn > dist[end] || dist[beg] < 0){
    fprintf(errfile,"sfs_code error(%ld): invCDF requesting invalid value\n",
	    INITIALSEED);
    fprintf(errfile,"rn=%f; dist[%ld]=%f; dist[%ld]=%1.9f\n",rn,beg,
	    dist[beg], end, dist[end]);
    fflush(errfile);
    abort();
  }
  
  /* use an inverse-CDF method to get index */
  (site) = (long)floor(rn/dist[end]*end+0.5);   /* initial guess in [0,end] */
  if(site > end){
    site = end;
  }
  else if(site < beg){
    site = beg;
  }

  if(dist[site]>=rn){
    if(site == beg || dist[site-1] < rn){
      return(site);
    }
    else if(site > beg+1 && dist[site-2] < rn && dist[site-1]>=rn){
      return(site-1);
    }
  }
  else if(site < end && dist[site+1]>=rn){
    return(site+1);
  }

  /* take large steps for the first n rounds, dividing interval in half each time,
     then take small steps... */
  numBigSteps = log((end-beg+1)/4)/log(2); /* after this many, interval is ~4 steps */
  while(steps < MAXSTEPS){
    steps++;
    if(site > end){
      site = end;
    }
    else if(site < beg){
      site = beg;
    }
    if((site == beg && dist[site] >= rn) ||
       (site > beg && dist[site-1] < rn && dist[site]>=rn))
      break;
    else if(steps < numBigSteps && fabs(rn-dist[site]) > 2*avgStep){
      if(dist[site] < rn){
	toosmall = site;
	site += (long)(0.5*(toobig-site));
      }
      else{
	toobig = site;
	site -= (long)(0.5*(site-toosmall));
      }
    }
    else{
      if(rn > dist[site]){
	site++;
      }
      else{
	site--;
      }
    }
  }

  if(steps<MAXSTEPS)
    return(site);
  else{
    fprintf(errfile,"sfs_code error(%ld): too many steps in invCDF\n",
	    INITIALSEED);
    fprintf(errfile,"rn=%f; dist[%ld]=%f; dist[%ld]=%f; site=%ld:(%f,%f)\n",
	    rn,beg, dist[beg], end, dist[end], site, dist[site-1], dist[site]);
    fflush(errfile);
    abort();
    return -1;
  }
}

/* ------------------------------------------------------------------------- */

void getNucHistory(const unsigned long site, char *retNuc, 
		   struct history *BigHead, struct history **root, 
		   const long locus, int pop, int SPLAY)
{
#ifdef UNIT_TEST
  if(site>=gpars.L[locus]){
    fprintf(errfile,"what kind of site is %lu??? max[%ld] = %lu  (%ld)\n",
	    site,locus,gpars.L[locus],INITIALSEED);
    (*retNuc)='4';
    abort();
  }
#endif  
  if(BigHead == NULL || popn[pop].polySites[locus][site] == 0){
    if((*retNuc) == 'x')
      (*retNuc) = popn[pop].conSeq[locus][site];
  }
  else if(BigHead->event->site < site)
    getNucHistory(site, retNuc, BigHead->Rtree, root, locus, pop, SPLAY);
  else if(BigHead->event->site > site)
    getNucHistory(site, retNuc, BigHead->Ltree, root, locus, pop, SPLAY);
  else{
    if(BigHead->event->nonsyn == '0' || BigHead->event->nonsyn == '1'){
      (*retNuc)=BigHead->event->derNuc; /*search Rtree for more recent mutant*/
      if(SPLAY==1)
	splayHistory(BigHead, root, pop);
    }
    else{ /* indels don't provide nucleotide info at site */
      if((*retNuc) == 'x')
	(*retNuc) = popn[pop].conSeq[locus][site];
    }
    if(BigHead->Rtree != NULL)
      getNucHistory(site, retNuc, BigHead->Rtree, root, locus, pop, SPLAY);
  }
}

/* ------------------------------------------------------------------------- */

char retNucHistory(const unsigned long site, struct history **BigHead, 
		   const long locus, int pop, int SPLAY)
{
  char nuc = 'x';
  getNucHistory(site, &nuc, (*BigHead)->Rtree, BigHead, locus, pop, SPLAY);
  return nuc;
}

/* ------------------------------------------------------------------------- */

void getNewHitProbsReg(double **hpTMP, struct history **BigHead,
		       const long locus, const int pop, const long first,
		       const long last)
{
  int j;
  unsigned long i;
  char c0,c1,c2;
  
  if(first>=gpars.L[locus]-1)
    return; /* nothing to do */

  if(first<=0)
    i = 1;
  else if(first==1)
    i = 1;
  else
    i = first-1;
  
  c0 = retNucHistory(i-1, BigHead, locus, pop, 1);
  c1 = retNucHistory(i, BigHead, locus, pop, 0);
  c2 = retNucHistory(i+1, BigHead, locus, pop, 0);
  for(; i<=last+2 && i<gpars.L[locus]-1; i++){
    (*hpTMP)[i] = (*hpTMP)[i-1];
    for(j=0;j<4;j++)
      if(j!=c1-'0')
	(*hpTMP)[i] += (Q[16*(c0-'0')+4*(c1-'0')+(c2-'0')][j]);
    c0=c1;
    c1=c2;
    if(i<last+2 && i<gpars.L[locus]-2)
      c2 = retNucHistory(i+2, BigHead, locus, pop, i%100);
  }
  if(i >= gpars.L[locus]-1){
    (*hpTMP)[gpars.L[locus]-1] = (*hpTMP)[gpars.L[locus]-2];
  }
}

/* ------------------------------------------------------------------------- */

void rotateHistoryNode(struct history *node, struct history **tree, int pop)
{
  struct history *parent, *grandparent;
  long x=0;
  
  if (node->Parent == NULL) return;
  
  parent = node->Parent;
  grandparent = parent->Parent;
  
  if (parent->Ltree == node){
    parent->Ltree = node->Rtree;
    if (parent->Ltree != NULL) 
      parent->Ltree->Parent = parent;
    node->Rtree = parent;
  } 
  else if (parent->Rtree == node){
    parent->Rtree = node->Ltree;
    if (parent->Rtree != NULL) 
      parent->Rtree->Parent = parent;
    node->Ltree = parent;
  } 
  else{
    fprintf(errfile,"rotateHistoryNode error(%ld):parent children not right\n",
	    INITIALSEED);
    PrintHistTree((*tree)->Rtree, &x, pop);
    fprintf(errfile,"\n");
    x=0;
    PrintHistTree(node, &x, pop);
    fprintf(errfile,"\n\n");
    abort();
  }
  
  parent->Parent = node;
  node->Parent = grandparent;
  
  if (grandparent == NULL){
    (*tree)->Rtree = node;
    node->Parent = NULL;
  }
  else if (grandparent->Ltree == parent) 
    grandparent->Ltree = node;
  else if (grandparent->Rtree == parent) 
    grandparent->Rtree = node;
  else{
    fprintf(errfile,"rotateHistoryNode error(%ld):grandparent kids not right\n",
	    INITIALSEED);
    x=0;
    PrintHistTree((*tree)->Rtree, &x, pop);
    fprintf(errfile,"\n");
    x=0;
    PrintHistTree(node, &x, pop);
    fprintf(errfile,"\n\n");
    abort();
  }
}

/* ------------------------------------------------------------------------- */

void splayHistory(struct history *node, struct history **tree, int pop)
{
  /* this function actually runs with probability RandSplayP */
  if(ran1(&gpars.seed) < gpars.RandSplayP){
    struct history *parent, *grandparent;
    
    if (node == NULL) return;
    while(1){
      if (node->Parent == NULL) break;
      
      parent = node->Parent;
      grandparent = parent->Parent;
      
      /* If the node's parent is the root of the tree, do one rotation */
      if (grandparent == NULL) 
	rotateHistoryNode(node, tree, pop);
      
      /* If we have a zig-zig, then rotate my parent, then rotate me */
      else if ((parent->Ltree == node && grandparent->Ltree == parent) ||
	       (parent->Rtree == node && grandparent->Rtree == parent)){
	rotateHistoryNode(parent, tree, pop);
	rotateHistoryNode(node, tree, pop);
	
	/* If we have a zig-zag, then rotate me twice */
      }
      else{
	rotateHistoryNode(node, tree, pop);
	rotateHistoryNode(node, tree, pop);
      }
    }
  }
}

/* ------------------------------------------------------------------------- */

void getLargestHist(struct history *tree, struct history **node)
{
  if(tree != NULL){
    (*node)->Rtree = tree;
    getLargestHist(tree->Rtree,node);
  }
}

/* ------------------------------------------------------------------------- */

void getSmallestHist(struct history *tree, struct history **node)
{
  if(tree != NULL){
    (*node)->Rtree = tree;
    getSmallestHist(tree->Ltree,node);
  }
}

/* ------------------------------------------------------------------------- */

void swapHapFreq(struct history *tree, const int pop, const long xreg,
		 const long chrTo, const long chrFrom)
{
  int check = 0;
  long i;
  if(tree == NULL)
    return;
  for(i=tree->event->nextHap[pop]; i>=0; i--){
    if(tree->event->hapFreq[pop][i]==&popn[pop].numCopies[chrFrom][xreg]){
      tree->event->hapFreq[pop][i] = &popn[pop].numCopies[chrTo][xreg];
      check = 1;
      break;
    }
  }
  if(check == 0){
    for(i=tree->event->nextHap[pop]+1; i<tree->event->maxHaps[pop]; i++){
      if(tree->event->hapFreq[pop][i]==&popn[pop].numCopies[chrFrom][xreg]){
	tree->event->hapFreq[pop][i] = &popn[pop].numCopies[chrTo][xreg];
	check = 1;
	break;
      }
    }
  }
#ifdef UNIT_TEST
  if(check == 0){
    fprintf(errfile, "error in swapHapFreq, could not find nc[%d][%ld][%ld] in \
hapFreq for ind %ld\n",pop, chrFrom, xreg, chrTo);
    abort();
  }
#endif
  swapHapFreq(tree->Ltree, pop, xreg, chrTo, chrFrom);
  swapHapFreq(tree->Rtree, pop, xreg, chrTo, chrFrom);
}

/* ------------------------------------------------------------------------- */

void setHistoryParents(struct history *tree)
{
  if(tree != NULL){
    if(tree->Ltree != NULL){
      tree->Ltree->Parent = tree;
      setHistoryParents(tree->Ltree);
    }
    if(tree->Rtree != NULL){
      tree->Rtree->Parent = tree;
      setHistoryParents(tree->Rtree);
    }
  }
}

/* ------------------------------------------------------------------------- */

void emptyMut(const long mutreg, const int pop)
{
  long i, j;
  for(i=mutAr[mutreg]->mutIndex; i<mutAr[mutreg]->numMuts; i++){
    if(mutAr[mutreg]->muts[i]->event->free == '1'){ /* found one! */
      mutAr[mutreg]->mutIndex = i;
      mutAr[mutreg]->muts[i]->event->free = '0';
#ifdef UNIT_TEST
      for(j=0; j<mutAr[mutreg]->muts[i]->event->maxHaps[pop]; j++){
	if(mutAr[mutreg]->muts[i]->event->hapFreq[pop][j] != NULL){
	  fprintf(errfile,"error in emptyMut, not actually free!!\n");
	  abort();
	}
      }
      for(j=0; j<=popn[pop].maxchr; j++){
	int foo = 0;
	carryMut(popn[pop].BigHead[j][mutreg]->Rtree, i, &foo);
	if(foo > 0){
	  fprintf(errfile,"error in emptyMut, carryMut!!!\n");
	  abort();
	}
      }
#endif
      mutAr[mutreg]->muts[i]->event->checkPop = -1;
      mutAr[mutreg]->muts[i]->event->checkGen = -1;
      mutAr[mutreg]->muts[i]->event->checkStep = -1;
      return;
    }
  }
  /* reaching this point means there were no free mutations, create new space!*/
  mutAr[mutreg]->mutIndex = mutAr[mutreg]->numMuts;
  mutAr[mutreg]->numMuts++; /* arb. allocate space for 10 more mutations */
  assert(mutAr[mutreg]->muts = 
	 realloc(mutAr[mutreg]->muts, 
		 mutAr[mutreg]->numMuts*sizeof(struct history*)));
  while(i<mutAr[mutreg]->numMuts){
    mutAr[mutreg]->muts[i] = popTrash();
    assert(mutAr[mutreg]->muts[i]->event = malloc(sizeof(struct event)));
    assert(mutAr[mutreg]->muts[i]->event->genFix =
	   malloc(gpars.NPOP*sizeof(long)));
    assert(mutAr[mutreg]->muts[i]->event->genDead =
	   malloc(gpars.NPOP*sizeof(long)));
    mutAr[mutreg]->muts[i]->event->nSites = 0;
    mutAr[mutreg]->muts[i]->event->nucs = NULL;
    assert(mutAr[mutreg]->muts[i]->event->fixed =
	   malloc(gpars.NPOP*sizeof(char)));
    mutAr[mutreg]->muts[i]->event->site = 0;
    mutAr[mutreg]->muts[i]->event->free = '1';
    mutAr[mutreg]->muts[i]->event->index = i;
    mutAr[mutreg]->muts[i]->event->checkPop = -1;
    mutAr[mutreg]->muts[i]->event->checkGen = -1;
    mutAr[mutreg]->muts[i]->event->checkStep = -1;
    assert(mutAr[mutreg]->muts[i]->event->nextHap =
	   malloc(gpars.NPOP*sizeof(long)));
    assert(mutAr[mutreg]->muts[i]->event->numCarriers =
	   malloc(gpars.NPOP*sizeof(long)));
    assert(mutAr[mutreg]->muts[i]->event->maxHaps =
	   malloc(gpars.NPOP*sizeof(long)));
    assert(mutAr[mutreg]->muts[i]->event->hapFreq =
	   malloc(gpars.NPOP*sizeof(long**)));
    for(j=0; j<gpars.NPOP; j++){
      mutAr[mutreg]->muts[i]->event->genFix[j] = 0;
      mutAr[mutreg]->muts[i]->event->genDead[j] = -LONG_MAX;
      mutAr[mutreg]->muts[i]->event->fixed[j] = '0';
      mutAr[mutreg]->muts[i]->event->nextHap[j] = 0;
      mutAr[mutreg]->muts[i]->event->numCarriers[j] = 0;
      mutAr[mutreg]->muts[i]->event->maxHaps[j] = 1;
      assert(mutAr[mutreg]->muts[i]->event->hapFreq[j] = malloc(sizeof(long*)));
      mutAr[mutreg]->muts[i]->event->hapFreq[j][0] = NULL;/* no carriers yet */
    }
    i++;
  }
  mutAr[mutreg]->muts[mutAr[mutreg]->mutIndex]->event->free = '0';
  mutAr[mutreg]->muts[mutAr[mutreg]->mutIndex]->event->checkPop = -1;
  mutAr[mutreg]->muts[mutAr[mutreg]->mutIndex]->event->checkGen = -1;
  mutAr[mutreg]->muts[mutAr[mutreg]->mutIndex]->event->checkStep = -1;
}

/* ------------------------------------------------------------------------- */

void addHistoryNode(const struct event *data, struct history **BigHead, 
		    struct history **MOM, const int pop, const long chr,
		    const long locus,  const int k)
{
  long i;
#ifdef UNIT_TEST
  long x;
  if(data->site >= gpars.L[locus]){
    fprintf(errfile,"site: %ld, nuc: %c, fit: %f\n",
	    data->site,data->ancNuc,data->fit);
    x=0;
    PrintHistTree((*BigHead), &x, pop);
    fprintf(errfile,"\nsomething screwy in addHistoryNode (%ld)\n",INITIALSEED);
    abort();
  }
#endif
  if((*BigHead) == NULL){ /* reached leaf node */
    (*BigHead) = popTrash();

    if(k == 0){ /* find empty mutation in mutAr[locus] (new or dead) */
      emptyMut(locus, pop);
      (*BigHead)->event = mutAr[locus]->muts[mutAr[locus]->mutIndex]->event;
      memcpyEVENT((*BigHead)->event, data, pop);
      (*BigHead)->event->hapFreq[pop][0] = &popn[pop].numCopies[chr][locus];
      (*BigHead)->event->numCarriers[pop] += popn[pop].numCopies[chr][locus];
      
#ifdef VERBOSE_DEBUG
      fprintf(errfile,"mutation: index=%ld; chr=%ld; freq=%ld (%ld)\n",
	      (*BigHead)->event->index, chr, *(*BigHead)->event->hapFreq[pop][0],
	      (*BigHead)->event->numCarriers[pop]);
      fflush(errfile);
#endif
    }
    else{/* otherwise have already found location for mutation */
      (*BigHead)->event = mutAr[locus]->muts[mutAr[locus]->mutIndex]->event;
      for(i=(*BigHead)->event->nextHap[pop]; i<(*BigHead)->event->maxHaps[pop];
	  i++){
	if((*BigHead)->event->hapFreq[pop][i] == NULL){
	  break;
	}
      }
      if(i == (*BigHead)->event->maxHaps[pop]){ /* allocate a new haplotype */
	(*BigHead)->event->maxHaps[pop]++;
	assert((*BigHead)->event->hapFreq[pop] = 
	       realloc((*BigHead)->event->hapFreq[pop], (i+1)*sizeof(long*)));
      }
      (*BigHead)->event->nextHap[pop] = i;
      (*BigHead)->event->hapFreq[pop][i] = &popn[pop].numCopies[chr][locus];
      (*BigHead)->event->numCarriers[pop] += popn[pop].numCopies[chr][locus];
#ifdef UNIT_TEST
      if((*BigHead)->event->numCarriers[pop] == 0){
	fprintf(errfile,"error addHistoryNode: numCopies[%d][%ld][%ldl]=0!\n",
		pop, chr, locus);
	abort();
      }
#endif
    }
    if(gpars.CONDTRAJ != NULL && gpars.condStarted==0 && 
       locus == gpars.condLoc && (*BigHead)->event->site == gpars.condSite){
      gpars.condMutIndex = (*BigHead)->event->index;
      gpars.condStarted = 1;
    }
    (*BigHead)->Rtree = NULL;
    (*BigHead)->Ltree = NULL;
    if(MOM != NULL)
      (*BigHead)->Parent = *MOM;
    else
      (*BigHead)->Parent = NULL;
    splayHistory(*BigHead, &popn[pop].BigHead[chr][locus], pop);
    return;
  }
  else{
    if((*BigHead)->event->site == data->site){
      if((*BigHead)->event->gen == data->gen && 
	 (*BigHead)->event->ancNuc == data->ancNuc &&
	 (*BigHead)->event->derNuc == data->derNuc &&
	 (*BigHead)->event->genFix[pop] == data->genFix[pop]){
	for(i=(*BigHead)->event->nextHap[pop];
	    i<(*BigHead)->event->maxHaps[pop]; i++){
	  if((*BigHead)->event->hapFreq[pop][i] == NULL){
	    break;
	  }
	}
	if(i == (*BigHead)->event->maxHaps[pop]){ /* allocate a new haplotype */
	  (*BigHead)->event->maxHaps[pop]++;
	  assert((*BigHead)->event->hapFreq[pop] = 
		 realloc((*BigHead)->event->hapFreq[pop], (i+1)*sizeof(long*)));
	}
	(*BigHead)->event->hapFreq[pop][i] = &popn[pop].numCopies[chr][locus];
	(*BigHead)->event->numCarriers[pop] += popn[pop].numCopies[chr][locus];
	splayHistory(*BigHead, &popn[pop].BigHead[chr][locus], pop);
      }
      else if((*BigHead)->event->gen <= data->gen) /* new data from later 
						      generation */ 
	addHistoryNode(data, &(*BigHead)->Rtree, BigHead, pop, chr, locus, k);
      else /* new event from earlier generation */
	addHistoryNode(data, &(*BigHead)->Ltree, BigHead, pop, chr, locus, k);
    }
    else if((*BigHead)->event->site < data->site)
      addHistoryNode(data, &(*BigHead)->Rtree, BigHead, pop, chr, locus, k);
    else
      addHistoryNode(data, &(*BigHead)->Ltree, BigHead, pop, chr, locus, k);
  }
}

/* ------------------------------------------------------------------------- */

void freeStorage(struct mutStorage *storage)
{
  if(storage != NULL){
    freeStorage(storage->Ltree);
    freeStorage(storage->Rtree);
    free(storage->pops);
    free(storage->chrs);
    free(storage->event->genFix);
    free(storage->event->genDead);
    if(storage->event->nucs != NULL)
      free(storage->event->nucs);
    free(storage->event);
    storage->event = NULL;
    storage->Ltree = NULL;
    storage->Rtree = NULL;
    storage->Parent = NULL;
    free(storage);
    storage = NULL;
  }
}

/* ------------------------------------------------------------------------- */

void printStorage(struct mutStorage *storage, long *x, FILE *out, long sampTime)
{
  long i;
  if(storage != NULL && sampTime == -1){
    printStorage(storage->Ltree, x, out, sampTime);
    if(gpars.PRINTLOC[storage->locus] == '1'){
      if(fabs(storage->event->fit) <= DBL_MIN){
	if(storage->event->nonsyn == '0' || storage->event->nonsyn == '1'){
	  fprintf(out,"%ld,%c,%lu,%ld,%ld,%c%c%c,%c,%c,%c,%c,%1.1f,",
		  storage->locus, storage->event->axy, storage->event->site, 
		  storage->event->gen, storage->event->genFix[0], 
		  ToCGTA(storage->event->fiveP),ToCGTA(storage->event->ancNuc),
		  ToCGTA(storage->event->threeP),ToCGTA(storage->event->derNuc),
		  storage->event->nonsyn, NumToAA(storage->event->ancAA), 
		  NumToAA(storage->event->derAA), storage->event->fit);
	  if(!gpars.allCodom)
	    fprintf(out,"%1.1f,",storage->event->h);
	  fprintf(out,"%ld", storage->numCarry);
	}
	else{
	  if(storage->event->nonsyn == 'i'){
	    for(i=0; i<storage->event->nSites; i++)
	      storage->event->nucs[i] = ToCGTA(storage->event->nucs[i]);
	    fprintf(out,"%ld,%c,%lu,%ld,%ld,%c-%c,%c,%c,%d,%s,%1.1f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], 
		    ToCGTA(storage->event->fiveP), 
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites, 
		    storage->event->nucs, storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.1f,",storage->event->h);
	    fprintf(out,"%ld", storage->numCarry);
	  }
	  else{
	    fprintf(out,"%ld,%c,%lu,%ld,%ld,%c+%c,%c,%c,%d,%1.1f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], 
		    ToCGTA(storage->event->fiveP),  
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites,
		    storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.1f,",storage->event->h);
	    fprintf(out,"%ld", storage->numCarry);
	  }
	}
      }
      else{
	if(storage->event->nonsyn == '0' || storage->event->nonsyn == '1'){
	  fprintf(out,"%ld,%c,%lu,%ld,%ld,%c%c%c,%c,%c,%c,%c,%e,",
		  storage->locus, storage->event->axy, storage->event->site, 
		  storage->event->gen, storage->event->genFix[0], 
		  ToCGTA(storage->event->fiveP), ToCGTA(storage->event->ancNuc),
		  ToCGTA(storage->event->threeP),ToCGTA(storage->event->derNuc),
		  storage->event->nonsyn, NumToAA(storage->event->ancAA), 
		  NumToAA(storage->event->derAA), storage->event->fit);
	  if(!gpars.allCodom)
	    fprintf(out,"%1.4f,",storage->event->h);
	  fprintf(out,"%ld",storage->numCarry);
	}
	else{
	  if(storage->event->nonsyn == 'i'){
	    for(i=0; i<storage->event->nSites; i++)
	      storage->event->nucs[i] = ToCGTA(storage->event->nucs[i]);
	    fprintf(out,"%ld,%c,%lu,%ld,%ld,%c-%c,%c,%c,%d,%s,%1.1f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], 
		    ToCGTA(storage->event->fiveP),  
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites, 
		    storage->event->nucs, storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.1f,",storage->event->h);
	    fprintf(out,"%ld",storage->numCarry);
	  }
	  else{
	    fprintf(out,"%ld,%c,%lu,%ld,%ld,%c+%c,%c,%c,%d,%.9f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], 
		    ToCGTA(storage->event->fiveP), 
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites,
		    storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%.9f,",storage->event->h);
	    fprintf(out, "%ld", storage->numCarry);
	  }
	}
      }
      for(i=0; i<storage->numCarry; i++)
	fprintf(out,",%d.%ld",storage->pops[i], storage->chrs[i]);
      fprintf(out,";");
      if((++(*x))%20 == 0 && storage->Rtree != NULL)
	fprintf(out,"\n");
    }
    printStorage(storage->Rtree, x, out, sampTime);
  }
  else if(storage != NULL && sampTime != -1){
    printStorage(storage->Ltree, x, out, sampTime);
    if(gpars.PRINTLOC[storage->locus] == '1'){
      if(fabs(storage->event->fit) <= DBL_MIN){
	if(storage->event->nonsyn == '0' || storage->event->nonsyn == '1'){
	  fprintf(out,"%ld,%c,%lu,%ld,%ld.%ld,%c%c%c,%c,%c,%c,%c,%1.1f,",
		  storage->locus, storage->event->axy, storage->event->site, 
		  storage->event->gen, storage->event->genFix[0], sampTime,
		  ToCGTA(storage->event->fiveP), ToCGTA(storage->event->ancNuc), 
		  ToCGTA(storage->event->threeP), ToCGTA(storage->event->derNuc),
		  storage->event->nonsyn, NumToAA(storage->event->ancAA), 
		  NumToAA(storage->event->derAA), storage->event->fit); 
	  if(!gpars.allCodom)
	    fprintf(out,"%1.1f,",storage->event->h);
	  fprintf(out,"%ld",storage->numCarry);
	}
	else{
	  if(storage->event->nonsyn == 'i'){
	    for(i=0; i<storage->event->nSites; i++)
	      storage->event->nucs[i] = ToCGTA(storage->event->nucs[i]);
	    fprintf(out,"%ld,%c,%lu,%ld,%ld.%ld,%c-%c,%c,%c,%d,%s,%1.1f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], sampTime,
		    ToCGTA(storage->event->fiveP),
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites, 
		    storage->event->nucs, storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.1f,",storage->event->h);
	    fprintf(out,"%ld",storage->numCarry);
	  }
	  else{
	    fprintf(out,"%ld,%c,%lu,%ld,%ld.%ld,%c+%c,%c,%c,%d,%1.1f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], sampTime,
		    ToCGTA(storage->event->fiveP),
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites,
		    storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.1f,",storage->event->h);
	    fprintf(out,"%ld",storage->numCarry);
	  }
	}
      }
      else{
	if(storage->event->nonsyn == '0' || storage->event->nonsyn == '1'){
	  fprintf(out,"%ld,%c,%lu,%ld,%ld.%ld,%c%c%c,%c,%c,%c,%c,%e,",
		  storage->locus, storage->event->axy, storage->event->site, 
		  storage->event->gen, storage->event->genFix[0], sampTime,
		  ToCGTA(storage->event->fiveP), ToCGTA(storage->event->ancNuc), 
		  ToCGTA(storage->event->threeP), ToCGTA(storage->event->derNuc),
		  storage->event->nonsyn, NumToAA(storage->event->ancAA), 
		  NumToAA(storage->event->derAA), storage->event->fit); 
	  if(!gpars.allCodom)
	    fprintf(out,"%1.4f,",storage->event->h);
	  fprintf(out,"%ld",storage->numCarry);
	}
	else{
	  if(storage->event->nonsyn == 'i'){
	    for(i=0; i<storage->event->nSites; i++)
	      storage->event->nucs[i] = ToCGTA(storage->event->nucs[i]);
	    fprintf(out,"%ld,%c,%lu,%ld,%ld.%ld,%c-%c,%c,%c,%d,%s,%1.1f,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], sampTime,
		    ToCGTA(storage->event->fiveP),
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites, 
		    storage->event->nucs, storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.1f,",storage->event->h);
	    fprintf(out,"%ld",storage->numCarry);
	  }
	  else{
	    fprintf(out,"%ld,%c,%lu,%ld,%ld.%ld,%c+%c,%c,%c,%d,%e,",
		    storage->locus, storage->event->axy, storage->event->site, 
		    storage->event->gen, storage->event->genFix[0], sampTime,
		    ToCGTA(storage->event->fiveP),
		    ToCGTA(storage->event->threeP), storage->event->derNuc, 
		    storage->event->nonsyn, storage->event->nSites,
		    storage->event->fit);
	    if(!gpars.allCodom)
	      fprintf(out,"%1.4f,",storage->event->h);
	    fprintf(out,"%ld",storage->numCarry);
	  }
	}
      }
      for(i=0; i<storage->numCarry; i++)
	fprintf(out,",%d.%ld",storage->pops[i], storage->chrs[i]);
      fprintf(out,";");
      if((++(*x))%20 == 0 && storage->Rtree != NULL)
	fprintf(out,"\n");
    }
    printStorage(storage->Rtree, x, out, sampTime);
  }
}

/* ------------------------------------------------------------------------- */

void printStorageVCF(struct mutStorage *storage,long *x,FILE *out,long sampTime)
{
  long i,j,AF=0;
  int pop,p,PNT;
  float TOTSS = 0.0;
  if(storage == NULL)
    return;
  
  printStorageVCF(storage->Ltree, x, out, sampTime);
  if(gpars.PRINTLOC[storage->locus] == '1'){
    for(i=0; i<gpars.NPOP; i++){
      TOTSS += gpars.P*ppars[i].SS;
    }
    for(i=0; i<storage->numCarry; i++){
      if(storage->chrs[i] == -1)
	AF += gpars.P*ppars[storage->pops[i]].SS;
      else
	AF++;
    }
    if(storage->event->nonsyn=='0' || storage->event->nonsyn=='1'){
      fprintf(out, "%ld\t%ld\t.\t%c\t%c\t.\tPASS\t", storage->locus,
	      storage->event->site, ToCGTA(storage->event->ancNuc),
	      ToCGTA(storage->event->derNuc));
      fprintf(out, "NS=%ld;",(long)TOTSS);
      
      fprintf(out,"AF=%.3f;",AF/TOTSS);
      fprintf(out,"AXY=%c;GA=%ld;GF=%ld;5P=%c;3P=%c;NSYN=%c;",
	      storage->event->axy, storage->event->gen,
	      storage->event->genFix[0], ToCGTA(storage->event->fiveP),
	      ToCGTA(storage->event->threeP), storage->event->nonsyn);
      if(gpars.ANNOTATE[storage->locus] == 'C')
	fprintf(out,"AA=%c;DA=%c;",
		NumToAA(storage->event->ancAA), NumToAA(storage->event->derAA));
    }
    else{ /* indel... */
      for(i=0; i<storage->event->nSites; i++)
	storage->event->nucs[i] = ToCGTA(storage->event->nucs[i]);
      fprintf(out, "%ld\t%ld\t.\t", storage->locus,
	      storage->event->site-1);
      if(storage->event->nonsyn=='i'){
	fprintf(out,"%c\t", ToCGTA(storage->event->fiveP));
	fprintf(out,"%c%s\t.\tPASS\t",ToCGTA(storage->event->fiveP),
		storage->event->nucs);
      }
      else{
	fprintf(out,"%c%s\t", ToCGTA(storage->event->ancNuc),
		storage->event->nucs);
	fprintf(out,"%c\t.\tPASS\t",ToCGTA(storage->event->fiveP));
      }
      fprintf(out, "NS=%ld;",(long)TOTSS);
      fprintf(out,"AF=%.3f;",AF/TOTSS);
      fprintf(out,"AXY=%c;GA=%ld;GF=%ld;",storage->event->axy,
	      storage->event->gen,storage->event->genFix[0]);
    }

    if(fabs(storage->event->fit) > FLT_MIN)
      fprintf(out,"H=%.4f;S=%.4f", storage->event->h, storage->event->fit);
    else
      fprintf(out,"H=1.0;S=0.0");
    fprintf(out,"\tGT");
    j=0;
    for(pop=0; pop<gpars.NPOP; pop++){
      for(i=0; i<ppars[pop].SS; i++){
	fprintf(out,"\t");
	for(p=0; p<gpars.P; p++){
	  PNT=0;
	  do{
	    if(j<storage->numCarry){
	      if(storage->pops[j] == pop &&
		 (storage->chrs[j] == -1 || storage->chrs[j] == gpars.P*i+p)){
		fprintf(out,"1");
		PNT=1;
		break;
	      }
	      else if(storage->pops[j] > pop || storage->chrs[j] > gpars.P*i+p){
		fprintf(out,"0");
		PNT=1;
		break;
	      }
	      else if(storage->pops[j] < pop || (storage->pops[j]==pop && storage->chrs[j] < gpars.P*i+p)){
		j++;
	      }
	    }
	    else{
	      fprintf(out,"0");
	      PNT=1;
	    }
	  }while(j<storage->numCarry);
	  if(PNT == 0){
	    fprintf(out,"0");
	  }
	  if(p<gpars.P-1){
	    fprintf(out, "|");
	  }
	}
      }
    }
    fprintf(out,"\n");
  }
  printStorageVCF(storage->Rtree, x, out, sampTime);
}

/* ------------------------------------------------------------------------- */

void storeInfo(struct mutStorage **storage, struct history *data, const int pop,
	       const long chr, const long locus, const long offset,
	       const long genFix)
{
  if(data == NULL)
    return;
  storeInfo(storage, data->Ltree, pop, chr, locus, offset, genFix);
  storeInfo(storage, data->Rtree, pop, chr, locus, offset, genFix);
  if(data->event != NULL){
    addToStorage(storage, data->event, pop, chr, locus, offset, 
		 (data->event->genFix[pop] == 0 ? 
		  genFix : data->event->genFix[pop]));
  }
}

/* ------------------------------------------------------------------------- */

void addToStorage(struct mutStorage **storage, struct event *data,
		  const int pop, const long chr, const long locus,
		  const long offset, const long genFix)
{
  int i;
  if((*storage) == NULL){ /* reached leaf node */
    assert((*storage) = malloc(sizeof(struct mutStorage)));
    assert((*storage)->event = malloc(sizeof(struct event)));
    assert((*storage)->event->genFix = malloc(gpars.NPOP*sizeof(long)));
    assert((*storage)->event->genDead = malloc(gpars.NPOP*sizeof(long)));
    (*storage)->event->nucs = NULL;
    for(i=0; i<gpars.NPOP; i++){
      (*storage)->event->genFix[i] = 0;
      (*storage)->event->genDead[i] = -LONG_MAX;
    }
    (*storage)->numCarry = 1;
    assert((*storage)->pops = malloc(sizeof(int)));
    assert((*storage)->chrs = malloc(sizeof(long)));
    (*storage)->pops[0] = pop;
    (*storage)->chrs[0] = chr;
    (*storage)->locus = locus;
    memcpyEVENT((*storage)->event, data, pop);
    (*storage)->event->gen -= offset;
    if(genFix != 0)
      (*storage)->event->genFix[0] = genFix;
    else  /* substitutions at different generations in different structures */
      (*storage)->event->genFix[0] = data->genFix[pop];
    (*storage)->Rtree = NULL;
    (*storage)->Ltree = NULL;
  }
  else{
    if((*storage)->locus == locus && (*storage)->event->site == data->site &&
       (*storage)->event->gen+offset == data->gen && 
       ((*storage)->event->genFix[0] == data->genFix[pop] ||
	(*storage)->event->genFix[0] == genFix) &&
       fabs((*storage)->event->fit - data->fit) <= DBL_MIN &&
       (*storage)->event->derNuc == data->derNuc && 
       (*storage)->event->ancNuc == data->ancNuc){ /* match */
      (*storage)->numCarry++;
      assert((*storage)->pops = 
	     realloc((*storage)->pops,
		     (*storage)->numCarry*sizeof(*(*storage)->pops)));
      (*storage)->pops[(*storage)->numCarry-1] = pop;
      assert((*storage)->chrs = 
	     realloc((*storage)->chrs,
		     (*storage)->numCarry*sizeof(*(*storage)->chrs)));
      (*storage)->chrs[(*storage)->numCarry-1] = chr;
    }
    else if((*storage)->locus < locus)
      addToStorage(&(*storage)->Rtree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->locus > locus)
      addToStorage(&(*storage)->Ltree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->site < data->site)
      addToStorage(&(*storage)->Rtree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->site > data->site)
      addToStorage(&(*storage)->Ltree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->gen+offset < data->gen)
      addToStorage(&(*storage)->Rtree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->gen+offset > data->gen)
      addToStorage(&(*storage)->Ltree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->genFix[0] < data->genFix[pop])
      addToStorage(&(*storage)->Rtree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->genFix[0] > data->genFix[pop])
      addToStorage(&(*storage)->Ltree, data, pop, chr, locus, offset, genFix);
    else if((*storage)->event->fit < data->fit)
      addToStorage(&(*storage)->Rtree, data, pop, chr, locus, offset, genFix);
    else /* new data into Ltree */
      addToStorage(&(*storage)->Ltree, data, pop, chr, locus, offset, genFix);
  }
}

/* ------------------------------------------------------------------------- */

void outputSample(long *SS){
  struct mutStorage *tsampStore = NULL;
  long *samps=NULL;
  long nF=0, nM=0, samp, i, j, r, pop;
  int cnt;
  
  if(tmpOUT != NULL){
    fprintf(tmpOUT, "MALES:");
  }
  else{
    fprintf(outfile, "MALES:");
  }

  for(pop=0; pop<gpars.NPOP; pop++){
    if(pop >= gpars.NPOPDEF){
      break;
    }

    nF=nM=0;
    /* pick random samples */
    assert(samps = malloc(SS[pop]*sizeof(*samps)));
    for(i=0; i<SS[pop]; i++){
      cnt=1;
      while(cnt==1){
	cnt=0;
	samp = (long)(ran1(&gpars.seed)*ppars[pop].N);
	for(j=0; j<nF; j++){
	  if(samps[j] == samp){
	    cnt = 1;
	    break;
	  }
	}
	for(j=(long)(ppars[pop].pFEMALES*SS[pop]+0.5); 
	    j<(long)(ppars[pop].pFEMALES*SS[pop]+0.5)+nM; j++){
	  if(samps[j] == samp){
	    cnt = 1;
	    break;
	  }
	}
	if(cnt == 0){
	  if(samp < ppars[pop].MALES){
	    if(nF >= ppars[pop].pFEMALES*SS[pop]){
	      cnt = 1;
	      continue;
	    }
	    else{
	      samps[nF] = samp;
	      nF++;
	      cnt=0;
	      break;
	    }
	  }
	  else{
	    if(nM >= (1-ppars[pop].pFEMALES)*SS[pop]){
	      cnt = 1;
	      continue;
	    }
	    else{
	      samps[(long)(ppars[pop].pFEMALES*SS[pop]+0.5)+nM] = samp;
	      nM++;
	      cnt=0;
	      break;
	    }
	  }
	}
      }
    }

    if(tmpOUT != NULL){
      if(pop < gpars.NPOP-1)
	fprintf(tmpOUT, "%ld,", gpars.P*nF);
      else
	fprintf(tmpOUT, "%ld", gpars.P*nF);
    }
    else{
      if(pop < gpars.NPOP-1)
	fprintf(outfile, "%ld,", gpars.P*nF);
      else
	fprintf(outfile, "%ld", gpars.P*nF);
    }
    
    for(r=0; r<gpars.R; r++){
      /* get fixed sites */
      storeInfo(&tsampStore, popn[pop].fixed[r]->Rtree, pop, -1, r, 0,
		popn[pop].gen);
      /* now get variants carried by each sample */
      for(i=0; i<SS[pop]; i++){
	storeInfo(&tsampStore,
		  popn[pop].BigHead[popn[pop].parentLoc[gpars.P*samps[i]][r]][r]->Rtree,
		  pop, gpars.P*i, r, 0, popn[pop].gen);
	storeInfo(&tsampStore,
		  popn[pop].BigHead[popn[pop].parentLoc[gpars.P*samps[i]+1][r]][r]->Rtree,
		  pop, gpars.P*i+1, r, 0, popn[pop].gen);
      }
    }
    free(samps);
  }
  if(tmpOUT != NULL){
    fprintf(tmpOUT, "\n");
  }
  else{
    fprintf(outfile, "\n");
  }
  
  i=0;
  if(tmpOUT != NULL){
    printStorage(tsampStore, &i, tmpOUT, -1);
    fprintf(tmpOUT, "\n");
  }
  else{
    printStorage(tsampStore, &i, outfile, -1);
    fprintf(outfile, "\n");
  }
  freeStorage(tsampStore);
}

/* ------------------------------------------------------------------------- */

void getPolySites(struct mutStorage *storage, long loc, long *Nsites, 
		  long **sites)
{
  long i, keep;
  if(storage != NULL){
    if(storage->locus >= loc)  /* recurse over storage tree */
      getPolySites(storage->Ltree, loc, Nsites, sites);
    if(storage->locus == loc){
      if((*Nsites) == 0){ /* first polymorphic site */
	assert((*sites) = malloc(sizeof(**sites)));
	(*sites)[0] = storage->event->site;
	(*Nsites) = 1;
      }
      else{
	keep = 1;
	for(i=0; i<*Nsites; i++) /* check that it doesn't already exist */
	  if(storage->event->site == (*sites)[i]){
	    keep = 0;
	    break;
	  }
	if(keep){
	  (*Nsites)++;
	  assert((*sites) = realloc((*sites), (*Nsites)*sizeof(**sites)));
	  (*sites)[(*Nsites)-1] = storage->event->site;
	}
      }
    }
    if(storage->locus <= loc)
      getPolySites(storage->Rtree, loc, Nsites, sites);
  }
}

/* ------------------------------------------------------------------------- */

void getHap(struct mutStorage *storage, long *Ncarry, int *carry, char **hap,
	    long locus, long chr, int pop, long *sites, long Nsites)
{
  long i, j, k, keep;
  if(storage != NULL){
    if(storage->Ltree != NULL)
      getHap(storage->Ltree, Ncarry, carry, hap, locus, chr, pop, sites,Nsites);
    if(storage->locus == locus){
      keep = 0; /* see if chr carries this mutation */
      for(i=0; i<storage->numCarry; i++){
	if(storage->pops[i] != pop)  continue;
	if(storage->chrs[i] == chr){
	  keep = 1;
	  if(storage->numCarry == 1){
	    (*Ncarry) = 1;
	    carry[0] = chr;
	  }
	  /* update haplotype */
	  for(j=0; j<Nsites; j++){
	    if(storage->event->site == sites[j]){
	      (*hap)[j] = storage->event->derNuc;
	      break;
	    }
	  }
	  break;
	}
      }
      if(keep){   /* ensure all carriers of haplotype carry mutation */
	if((*Ncarry) == 0){
	  for(j=0; j<storage->numCarry; j++)
	    if(storage->pops[j] == pop){
	      keep = 1;
	      for(k=0; k<j; k++)
		if(storage->chrs[j] == carry[k]){
		  keep = 0;
		  break;
		}
	      if(keep)
		carry[(*Ncarry)++] = storage->chrs[j];
	    }
	}
	else{
	  for(j=0; j<(*Ncarry); j++){
	    keep = 0;
	    for(k=0; k<storage->numCarry; k++)
	      if(storage->pops[k] == pop && storage->chrs[k] == carry[j]){
		keep = 1;
		break;
	      }
	    if(!keep){ /* chromosome does not carry mutation on haplotype */
	      (*Ncarry)--;
	      for(k=j; k<(*Ncarry); k++)
		carry[k] = carry[k+1];
	      j--;
	    }
	  }
	}
      }
      else{  /* ensure that individuals do not carry extra mutations */
	for(j=0; j<(*Ncarry); j++){
	  keep = 1;
	  for(k=0; k<storage->numCarry; k++)
	    if(storage->pops[k] == pop && storage->chrs[k] == carry[j]){
	      keep = 0;
	      break;
	    }
	  if(!keep){ /* chromosome carries extra mutations */
	    (*Ncarry)--;
	    for(k=j; k<(*Ncarry); k++)
	      carry[k] = carry[k+1];
	    j--;
	  }
	}
      }
    }
    if(storage->Rtree != NULL)
      getHap(storage->Rtree, Ncarry, carry, hap, locus, chr, pop, sites,Nsites);
  }
}

/* ------------------------------------------------------------------------- */

void copyStorageToBigHead(struct mutStorage *storage, int pop, long chr,
			  long ind, long locus)
{
  int k;
  long i;
  if(storage == NULL)
    return;
  if(storage->locus >= locus)
    copyStorageToBigHead(storage->Ltree, pop, chr, ind, locus);
  if(storage->locus <= locus)
    copyStorageToBigHead(storage->Rtree, pop, chr, ind, locus);
  if(storage->locus == locus){
    for(i=0; i<storage->numCarry; i++){
      if(storage->pops[i] == pop && storage->chrs[i] == ind){
	/* if mutation already exists in population, point to it! */
	k = 0;
	for(i=0; i<mutAr[locus]->numMuts; i++){
	  if(mutAr[locus]->muts[i]->event->free == '0'){ /* not free */
	    if(storage->event->site == mutAr[locus]->muts[i]->event->site &&
	       storage->event->gen == mutAr[locus]->muts[i]->event->gen && 
	       storage->event->ancNuc == mutAr[locus]->muts[i]->event->ancNuc &&
	       storage->event->derNuc == mutAr[locus]->muts[i]->event->derNuc){
	      /* mutation already exists! */
	      k = 1;
	      mutAr[locus]->mutIndex = i;
	      break;
	    }
	  }
	}
	addHistoryNode(storage->event, &popn[pop].BigHead[chr][locus]->Rtree,
		       NULL, pop, chr, locus, k);
	if(k == 0 && (storage->event->nonsyn == '0' ||
		      storage->event->nonsyn == '1')){ /* adding new mutation */
	  popn[pop].polySites[locus][storage->event->site]++;
	  storage->event->genDead[pop] = -LONG_MAX;
#ifdef VERBOSE_DEBUG
	  fprintf(errfile,"update PS[%ld]=%d\n",storage->event->site, 
		  popn[pop].polySites[locus][storage->event->site]);
	  if(popn[pop].polySites[locus][storage->event->site]>1){
	    //getchar();
	  }
#endif
	}
	break;
      }
    }
  }
}

/* ------------------------------------------------------------------------- */


void checkIDhaplotype(struct mutStorage *storage, int pop, long ref, long test, 
		      long locus, int *keep)
{
  long i, j;
  if(ref == test) return;
  else if(storage != NULL){
    if(storage->locus == locus){
      j = 0;
      for(i=0; i<storage->numCarry; i++){
	if(storage->pops[i] == pop &&
	   (storage->chrs[i] == test || storage->chrs[i] == ref))
	  j++;
      }
      if(j != 0 && j != 2){
	if(j != 1){
	  fprintf(errfile,"\nin checkIDhaplotype (%ld), j=%ld for pop%d, ref=%ld, test=%ld, locus=%ld, site=%lu!!  numCarry=%ld\n", INITIALSEED, j, pop, ref, test, locus,storage->event->site, storage->numCarry);
	  for(j=0; j<storage->numCarry; j++){
	    fprintf(errfile,"%d.%ld,",storage->pops[j],storage->chrs[j]);
	  }
	  fprintf(errfile,"\n");
	  fflush(errfile);
	  abort();
	}
	(*keep) = 0;
	return;
      }
    }
    if(storage->Ltree != NULL)
      checkIDhaplotype(storage->Ltree, pop, ref, test, locus, keep);
    if(storage->Rtree != NULL)
      checkIDhaplotype(storage->Rtree, pop, ref, test, locus, keep);
  }
}


/* ------------------------------------------------------------------------- */

void freeHistory(struct history *tree, const int pop, const long chr,
		 const long locus, const int step, const int inc)
{
  int D = (gpars.P == 4 ? 2 : 1);/*hap/diploid males have 1 X, tet have 2Xs*/
  long i, j;
  if(tree == NULL)
    return;
  freeHistory(tree->Ltree, pop, chr, locus, step, inc);
  freeHistory(tree->Rtree, pop, chr, locus, step, inc);

  if(tree->event == NULL)
    return;
#ifdef VERBOSE_DEBUG
  printf("freeHistory: site=%ld.%ld, free=%c, fixed=%d, checkGen=%ld(%ld), \
checkStep=%d(%d), checkPop=%d(%d), NC=%ld\n",
	 locus, tree->event->site, tree->event->free, tree->event->fixed[pop], tree->event->checkGen,
	 popn[pop].gen,tree->event->checkStep,step,tree->event->checkPop,pop, tree->event->numCarriers[pop]);
  fflush(stdout);
#endif
  if(tree->event->free == '0' && tree->event->fixed[pop] == '0' && 1
     /* (tree->event->checkGen != popn[pop].gen || tree->event->checkStep != step */
     /*  || tree->event->checkPop != pop) */){
    tree->event->checkGen = popn[pop].gen;
    tree->event->checkStep = step;
    tree->event->numCarriers[pop] = 0;
    for(i=0; i<tree->event->maxHaps[pop]; i++){
      if(tree->event->hapFreq[pop][i] != NULL){
	if((*tree->event->hapFreq[pop][i]) == 0){
	  tree->event->hapFreq[pop][i] = NULL;
	  if(tree->event->nextHap[pop] > i)
	    tree->event->nextHap[pop] = i;
	}
	else
	  tree->event->numCarriers[pop] += (*tree->event->hapFreq[pop][i]);
      }
    }
#ifdef VERBOSE_DEBUG
    if(tree->event->index == 4){
      fprintf(errfile,"(freeHistory) site=%ld; numCopies=%ld; polysites=%d\n",
	      tree->event->site, tree->event->numCarriers[pop], 
	      popn[pop].polySites[locus][tree->event->site]);
      fflush(errfile);
    }
#endif
    if(tree->event->numCarriers[pop] == 0){
      if(tree->event->nonsyn == '0' || tree->event->nonsyn == '1'){
	if(popn[pop].polySites[locus][tree->event->site]>0){
	  if(tree->event->genDead[pop] != popn[pop].gen){
	    popn[pop].polySites[locus][tree->event->site]--;
	    tree->event->genDead[pop] = popn[pop].gen;
#ifdef VERBOSE_DEBUG
	    fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!!\n",pop, locus,
		    tree->event->site,
		    popn[pop].polySites[locus][tree->event->site]);
	    fflush(errfile);
#endif
	  }
#ifdef VERBOSE_DEBUG
	  else{
	    fprintf(errfile,"ps[%d][%ld][%ld] -> %d  ALREADY!!!\n",pop, locus,
		    tree->event->site,
		    popn[pop].polySites[locus][tree->event->site]);
	    fflush(errfile);
	  }
#endif
	}
#ifdef VERBOSE_DEBUG
	else{
	  fprintf(errfile,"ps[%d][%ld][%ld] would -> <0 (%d)!!!\n",pop, locus,
		  tree->event->site,
		  popn[pop].polySites[locus][tree->event->site]);
	  fflush(errfile);
	}
#endif
      }
#ifdef UNIT_TEST
      {
	if(popn[pop].polySites[locus][tree->event->site] < 0){
	  fprintf(errfile,"ps[%d][%ld][%ld] = %d!!!\n",pop, locus,
		  tree->event->site,
		  popn[pop].polySites[locus][tree->event->site]);
	  fflush(errfile);
	  abort();
	}
#ifdef VERBOSE_DEBUG
	fprintf(errfile,"popn[%d].polySites[%ld][%ld] = %d\n",pop, locus,
		tree->event->site,popn[pop].polySites[locus][tree->event->site]);
	fflush(errfile);
#endif
	if(popn[pop].polySites[locus][tree->event->site] < 0){
	  fprintf(errfile,"ps[%d][%ld][%ld] = %d < 0!!\n",pop, locus, tree->event->site,
		  popn[pop].polySites[locus][tree->event->site]);
	  PrintErrorStats(pop, locus);
	  abort();
	}
	for(i=0; i<=popn[pop].maxchr; i++){
	  if(i == chr || popn[pop].numCopies[i][locus] == 0)
	    continue;
	  int foo = 0;
	  carryMut(popn[pop].BigHead[i][locus]->Rtree, tree->event->index,&foo);
	  if(foo > 0){
	    fprintf(errfile,"error! BigHead[%d][%ld][%ld] points to index %ld\n"
		    ,pop,i,locus,tree->event->index);
	    fflush(errfile);
	    j=0;
	    PrintHistTree(popn[pop].BigHead[i][locus]->Rtree, &j, pop);
	    fflush(errfile);
	    abort();
	  }
	}
      }
#endif
      j = 0;
      for(i=0; i<gpars.NPOP; i++){
	if(tree->event->numCarriers[i] != 0 || tree->event->fixed[i] == '1'){
	  j = 1;
	  break;
	}
      }
      if(j == 0){ /* free in all populations */
	tree->event->free = '1';
	if(tree->event->index < mutAr[locus]->mutIndex)
	  mutAr[locus]->mutIndex = tree->event->index;
      }
      tree->event->checkPop = pop;
      tree->event->checkGen = popn[pop].gen;
      tree->event->checkStep = -1;
    }
    else if((tree->event->axy == 'A' &&
	     tree->event->numCarriers[pop] != gpars.P*ppars[pop].N) ||
	    (tree->event->axy == 'X' &&        /* #female Xs + #male Xs */
	     tree->event->numCarriers[pop] !=
	     (gpars.P*ppars[pop].MALES+(ppars[pop].N-ppars[pop].MALES)*D)) ||
	    (tree->event->axy == 'Y' &&        /* just #male Ys */
	     tree->event->numCarriers[pop] !=
	     D*(ppars[pop].N-ppars[pop].MALES))){
      tree->event->checkPop = pop;
      tree->event->checkGen = popn[pop].gen;
      tree->event->checkStep = step;
    }
  }
  tree->event = NULL;
}

/* ------------------------------------------------------------------------- */

void pushTrash(struct history **tree)
{
  if(*tree != NULL){
    if(trash->Rtree != NULL)
      (*tree)->Parent = trash->Rtree;
    else
      (*tree)->Parent = NULL;
    trash->Rtree = (*tree);
  }
}

/* ------------------------------------------------------------------------- */

void createNewHistory(const long chr, const long locus,int pop, int IBD, int FZ)
{
  long i, tmpChr;
  /* find an unused chromosome */
  if(IBD==1 && FZ == 1 &&
     popn[pop].numCopies[popn[pop].parentLoc[chr][locus]][locus] <= 1)
    return; /* individual is already sole copy, nothing to do */
  else{
    for(tmpChr=popn[pop].nextchr[locus]; tmpChr<=popn[pop].maxchr; tmpChr++){
      if(popn[pop].numCopies[tmpChr][locus]==0 &&
	 popn[pop].BigHead[tmpChr][locus]->Rtree==NULL &&
	 tmpChr!=popn[pop].parentLoc[chr][locus] /* && */
/* 	 popn[pop].chrGen[tmpChr][locus]!=popn[pop].gen */){ /*first free element*/
	popn[pop].nextchr[locus] = tmpChr;
	break;
      }
    }
  }
  
#ifdef UNIT_TEST
  if(tmpChr > popn[pop].maxchr){
    fprintf(errfile,"new slot = %ld (old = %ld:%ld) ",
	    tmpChr,chr,popn[pop].parentLoc[chr][locus]);
    fprintf(errfile,"maxchr[p=%d;l=%ld] = %ld -> %ld ", pop, locus,
	    popn[pop].maxchr,popn[pop].maxchr+1);
    fprintf(errfile,"\n");
  }
  fflush(errfile);
#endif
  
  if(tmpChr > popn[pop].maxchr){  /* need to realloc several components */
    popn[pop].maxchr++;
    popn[pop].nextchr[locus] = 0;

    if(gpars.substMod == 5){
      assert(popn[pop].hpSUM = realloc(popn[pop].hpSUM,
				       (tmpChr+1)*sizeof(double**)));
      assert(popn[pop].hpSUM[tmpChr] = malloc(gpars.R*sizeof(double*)));
    }
    assert(popn[pop].numCopies = 
	   realloc(popn[pop].numCopies, (tmpChr+1)*sizeof(long*)));
    assert(popn[pop].chrGen = 
	   realloc(popn[pop].chrGen, (tmpChr+1)*sizeof(long*)));
    assert(popn[pop].numCopies[tmpChr] = malloc(gpars.R*sizeof(long)));
    assert(popn[pop].chrGen[tmpChr] = malloc(gpars.R*sizeof(long)));
    assert(popn[pop].indfit = 
	   realloc(popn[pop].indfit, (tmpChr+1)*sizeof(double*)));
    assert(popn[pop].indfit[tmpChr] = malloc(gpars.R*sizeof(double)));
    assert(popn[pop].BigHead = 
	   realloc(popn[pop].BigHead, (tmpChr+1)*sizeof(struct history**)));
    assert(popn[pop].BigHead[tmpChr] = malloc(gpars.R*sizeof(struct history*)));
    assert(popn[pop].extremeMuts = 
	   realloc(popn[pop].extremeMuts, (tmpChr+1)*sizeof(long**)));
    assert(popn[pop].extremeMuts[tmpChr] = malloc(gpars.R*sizeof(long*)));
    if(gpars.TRACKANC){
      assert(popn[pop].ancestry = 
	     realloc(popn[pop].ancestry, (tmpChr+1)*sizeof(int**)));
      assert(popn[pop].ancestry[tmpChr] = malloc(gpars.R*sizeof(int*)));
    }
    for(i=0; i<gpars.R; i++){
      if(gpars.substMod == 5){
	assert(popn[pop].hpSUM[tmpChr][i] = malloc(gpars.L[i]*sizeof(double)));
	popn[pop].hpSUM[tmpChr][i][0] = 0.0;
      }
      popn[pop].numCopies[tmpChr][i] = 0;
      popn[pop].chrGen[tmpChr][i] = popn[pop].gen;
      popn[pop].indfit[tmpChr][i] = (gpars.ADDITIVE ? 0.0 : 1.0);
      popn[pop].BigHead[tmpChr][i] = NULL;
      popn[pop].BigHead[tmpChr][i] = popTrash();
      assert(popn[pop].extremeMuts[tmpChr][i] = malloc(2*sizeof(long)));
      popn[pop].extremeMuts[tmpChr][i][0] = gpars.L[i];
      popn[pop].extremeMuts[tmpChr][i][1] = 0;
      if(gpars.TRACKANC){
	assert(popn[pop].ancestry[tmpChr][i] = malloc(gpars.L[i]*sizeof(int)));
      }
    }
  }
  
  if(IBD == 1){
#ifdef VERBOSE_DEBUG
    i=0;
    fprintf(errfile,"tree before copy: %ld, n=%ld\n",popn[pop].parentLoc[chr][locus],
	   popn[pop].numCopies[popn[pop].parentLoc[chr][locus]][locus]);
    PrintHistTree(popn[pop].BigHead[popn[pop].parentLoc
				    [chr][locus]][locus]->Rtree, &i, pop);
    fflush(errfile);
#endif
    if(popn[pop].BigHead[popn[pop].parentLoc[chr][locus]][locus]->Rtree!=NULL){
      copyHistory(&popn[pop].BigHead[tmpChr][locus],
		  &popn[pop].BigHead[popn[pop].parentLoc[chr][locus]][locus],
		  locus, pop, tmpChr, pop, 0, 0, 0);
    }
    if(gpars.substMod == 5){
      memcpy(popn[pop].hpSUM[tmpChr][locus], 
	     popn[pop].hpSUM[popn[pop].parentLoc[chr][locus]][locus],
	     gpars.L[locus]*sizeof(double));
    }
    popn[pop].indfit[tmpChr][locus] =
      popn[pop].indfit[popn[pop].parentLoc[chr][locus]][locus];
    popn[pop].extremeMuts[tmpChr][locus][0] =
      popn[pop].extremeMuts[popn[pop].parentLoc[chr][locus]][locus][0];
    popn[pop].extremeMuts[tmpChr][locus][1] =
      popn[pop].extremeMuts[popn[pop].parentLoc[chr][locus]][locus][1];
    if(gpars.TRACKANC){
      memcpy(popn[pop].ancestry[tmpChr][locus], 
	     popn[pop].ancestry[popn[pop].parentLoc[chr][locus]][locus],
	     gpars.L[locus]*sizeof(int));
    }
  }
  popn[pop].numCopies[tmpChr][locus] = 1;
  popn[pop].numCopies[popn[pop].parentLoc[chr][locus]][locus]--;
  if(popn[pop].numCopies[popn[pop].parentLoc[chr][locus]][locus]==0){
    popn[pop].chrGen[popn[pop].parentLoc[chr][locus]][locus] = popn[pop].gen;
    newDead(popn[pop].parentLoc[chr][locus], locus);
  }
#ifdef VERBOSE_DEBUG
  fprintf(errfile,"after copy: n[%ld][%ld]=%ld, n[%ld][%ld]=%ld\n",
	 popn[pop].parentLoc[chr][locus], locus,
	 popn[pop].numCopies[popn[pop].parentLoc[chr][locus]][locus],
	 tmpChr, locus, popn[pop].numCopies[tmpChr][locus]);
  fflush(errfile);
  i=0;
  PrintHistTree(popn[pop].BigHead[tmpChr][locus]->Rtree, &i, pop);
  fflush(errfile);
#endif /* VERBOSE_DEBUG */
  if(FZ==1 && popn[pop].numCopies[popn[pop].parentLoc[chr][locus]][locus]==0){
    freeHistory(popn[pop].BigHead[popn[pop].parentLoc[chr][locus]]
		[locus]->Rtree,pop,popn[pop].parentLoc[chr][locus],locus,0,1);
    if(popn[pop].nextchr[locus] > popn[pop].parentLoc[chr][locus])
      popn[pop].nextchr[locus] = popn[pop].parentLoc[chr][locus];
    pushTrash(&popn[pop].BigHead[popn[pop].parentLoc[chr][locus]]
	      [locus]->Rtree);
    popn[pop].BigHead[popn[pop].parentLoc[chr][locus]][locus]->Rtree = NULL;
  }
  popn[pop].parentLoc[chr][locus] = tmpChr;
}

/* ------------------------------------------------------------------------- */

void regeneratePOP()
{
  int *carry;
  long i, j, r, *sites=NULL, Ncarry, Nsites=0;
  char *newHap;
  int *done, keep;
  
  assert(done = malloc(gpars.P*ppars[0].N*sizeof(*done)));
  assert(carry = malloc(gpars.P*ppars[0].N*sizeof(*carry)));
  for(i=0; i<gpars.P*ppars[0].N; i++)
    done[i] = 0;
  for(r=0; r<gpars.R; r++){
    /* first find all polymorphic sites */
    Nsites = 0;
    getPolySites(storage, r, &Nsites, &sites);

    /* generate unique haplotypes */
    assert(newHap = malloc((Nsites+1)*sizeof(*newHap)));
    for(i=0; i<gpars.P*ppars[0].N; i++){
      if(done[i] > r)  continue;
      for(j=0; j<Nsites; j++)  newHap[j] = '4';/* non-nuc value */
      newHap[Nsites] = '\0'; /* string termination char */
      Ncarry = 0;
      getHap(storage, &Ncarry, carry, &newHap, r, i, 0, sites, Nsites);
      done[i]++;
      if(Ncarry > 0){
	createNewHistory(i, r, 0, 0, 1);
	copyStorageToBigHead(storage, 0, popn[0].parentLoc[i][r], i, r);
	popn[0].indfit[popn[0].parentLoc[i][r]][r] = (gpars.ADDITIVE ? 0.0:1.0);
	getFitTree(popn[0].BigHead[popn[0].parentLoc[i][r]][r]->Rtree,
		   &popn[0].indfit[popn[0].parentLoc[i][r]][r], 0, r);
	for(j=0; j<Nsites; j++){
	  if(newHap[j] != '4'){
	    if(popn[0].extremeMuts[popn[0].parentLoc[i][r]][r][0]>sites[j])
	      popn[0].extremeMuts[popn[0].parentLoc[i][r]][r][0] = sites[j];
	    if(popn[0].extremeMuts[popn[0].parentLoc[i][r]][r][1]<sites[j])
	      popn[0].extremeMuts[popn[0].parentLoc[i][r]][r][1] = sites[j];
	  }
	}
	
	for(j=0; j<Ncarry; j++){
	  keep = 1;
	  checkIDhaplotype(storage, 0, i, carry[j], r, &keep);
	  if(!keep) continue;
	  if(carry[j] == i)  continue;
#ifdef UNIT_TEST
	  if(done[carry[j]] > r){
	    fprintf(errfile,"sfs_code error(%ld):carry[%ld]=%d already done!\n",
		    INITIALSEED, j, carry[j]);
	    fprintf(errfile,"done[%d] = %d\n",carry[j], done[carry[j]]);
	    abort();
	  }
#endif
	  done[carry[j]]++;
	  popn[0].numCopies[0][r]--;
	  if(popn[0].numCopies[0][r]==0)
	    popn[0].chrGen[0][r] = popn[0].gen;
	  popn[0].parentLoc[carry[j]][r] = popn[0].parentLoc[i][r];
	  popn[0].numCopies[popn[0].parentLoc[i][r]][r]++;
	}
      }
    }
    for(i=0; i<gpars.P*ppars[0].N; i++)
      if(done[i] <= r){
	fprintf(errfile,"DOH (%ld):  %ld is not done\n",INITIALSEED,i);
	abort();
      }
    /* free instead of having to realloc */
    free(sites);
    sites = NULL;
    free(newHap);
  }
  free(done);
  free(carry);
}

/* ------------------------------------------------------------------------- */

void renewConSeq()
{
  char *nucs, codon[3];
  int k;
  long i, r;
  double *ST;            /* absolute mut rates, transition probs, 
			    stat freqs, rel. mut rates, hit probs */
  assert(ST = calloc(64,sizeof(*ST)));
  GenContextRate(ST,NULL);
  assert(nucs=malloc(4*sizeof(char)));
  for(r = 0; r < gpars.R; r++){
    for(i=0;i<gpars.L[r];i++){
      if(i<gpars.L[r]-2 && popn[0].polySites[r][i+2] > 0)
	continue;
      else if(i<gpars.L[r]-1 && popn[0].polySites[r][i+1] > 0)
	continue;
      else if(popn[0].polySites[r][i] > 0)
	continue;
      else if(i>0 && popn[0].polySites[r][i-1] > 0)
	continue;
      else if(i>1 && popn[0].polySites[r][i-2] > 0)
	continue;
      else{ /* no polymorphisms w/in 2bp, update site */
	if(i==0)
	  nucs[0] = popn[0].conSeq[r][gpars.L[r]-1];
	else
	  nucs[0] = popn[0].conSeq[r][i-1];
	nucs[1] = popn[0].conSeq[r][i];
	if(i==gpars.L[r]-1)
	  nucs[2] = popn[0].conSeq[r][0];
	else
	  nucs[2] = popn[0].conSeq[r][i+1];
	k = 0;
	if((nucs[0]=='0' && nucs[1]=='1') || (nucs[1]=='0' && nucs[2]=='1'))
	  k=1;  /* CpG */

	popn[0].conSeq[r][i] = MutantNuc(nucs, k, 0);
	if(gpars.ANNOTATE[r]=='C'){
	  if(i%3 == 0){
	    codon[0] = popn[0].conSeq[r][i];
	    codon[1] = popn[0].conSeq[r][i+1];
	    codon[2] = popn[0].conSeq[r][i+2];
	  }
	  else if(i%3 == 1){
	    codon[0] = popn[0].conSeq[r][i-1];
	    codon[1] = popn[0].conSeq[r][i];
	    codon[2] = popn[0].conSeq[r][i+1];
	  }
	  else{
	    codon[0] = popn[0].conSeq[r][i-2];
	    codon[1] = popn[0].conSeq[r][i-1];
	    codon[2] = popn[0].conSeq[r][i];
	  }
	  k = 16*(codon[0]-'0') + 4*(codon[1]-'0') + (codon[2]-'0');
	  if(k == 39 || k == 45 || k == 47) {
	    /* stop codon, draw another nucleotide at this site */
	    i--;
	    continue;
	  }
	}
      }
    }
  }
  free(ST);
  free(nucs);
}

/* ------------------------------------------------------------------------- */

struct history* Delete_max_hist(struct history *tree, struct history **data,
				const long locus) 
{
  struct history *temp;
  if(tree->Rtree == NULL){
    (*data)->event = tree->event;
    temp = tree->Ltree;
    if(temp != NULL)
      temp->Parent = tree->Parent;
    tree->Rtree = NULL;
    tree->Ltree = NULL;
    free(tree);
    tree = NULL;
    return temp;
  } 
  else{
    tree->Rtree = Delete_max_hist(tree->Rtree, data, locus);
    return tree;
  }
}  /* Delete_max_hist */

/* ------------------------------------------------------------------------- */

struct history* Delete_node_hist(struct history *tree, const long locus)
{
  struct history *temp, *data;
  if(tree == NULL)
    return NULL;
  if((tree->Ltree == NULL) && (tree->Rtree == NULL)){
    free(tree);
    return NULL;
  }
  else if(tree->Ltree == NULL){
    temp = tree->Rtree;
    if(temp != NULL)
      temp->Parent = tree->Parent;
    else{
      fprintf(errfile, "error in Delete_node_hist, null Rtree!\n");
      abort();
    }
    tree->Rtree = NULL;
    tree->Ltree = NULL;
    free(tree);
    return temp;
  }
  else if (tree->Rtree == NULL){
    temp = tree->Ltree;
    if(temp != NULL)
      temp->Parent = tree->Parent;
    else{
      fprintf(errfile, "error in Delete_node_hist, null Rtree!\n");
      abort();
    }
    tree->Rtree = NULL;
    tree->Ltree = NULL;
    free(tree);
    return temp;
  }
  else{ /* node has two children */
    data = popTrash();
    tree->Ltree = Delete_max_hist(tree->Ltree, &data, locus);
    tree->event = NULL;
    tree->event = data->event;
    data->event = NULL;
    data->Rtree = NULL;
    data->Ltree = NULL;
    free(data);
    return tree;
  }
}  /* Delete_node_hist */


/* ------------------------------------------------------------------------- */

struct history* freeHistoryNode(struct event *event, struct history *tree,
				const long locus)
{
  struct history *tmp;
  if(tree == NULL)
    return NULL;
  else if(event->index == tree->event->index){
    if(tree->Rtree == NULL && tree->Ltree == NULL){
      free(tree);
      return NULL;
    }
    else{
      tmp = popTrash();
      tmp->Rtree = tree->Parent;
      tree = Delete_node_hist(tree, locus);
      if(tree != NULL)
	tree->Parent = tmp->Rtree;
      tmp->Rtree = NULL;
      free(tmp);
    }
  }
  else if (event->site < tree->event->site ||
	   (event->site == tree->event->site && event->gen < tree->event->gen))
    tree->Ltree = freeHistoryNode(event, tree->Ltree, locus);
  else if (event->site > tree->event->site ||
	   (event->site == tree->event->site && event->gen > tree->event->gen))
    tree->Rtree = freeHistoryNode(event, tree->Rtree, locus);
  else{
    tree->Ltree = freeHistoryNode(event, tree->Ltree, locus);
    tree->Rtree = freeHistoryNode(event, tree->Rtree, locus);
  }
  return tree;
}  /* freeHistoryNode */

/* ------------------------------------------------------------------------- */

struct history* freeHistoryNodeRange(const long min, const long max,
				     struct history *tree, const long locus,
				     const int pop, const long chr)
{
  long i;
  struct history *tmp;
  if(tree == NULL)
    return NULL;
  if(min <= tree->event->site)
    tree->Ltree = freeHistoryNodeRange(min, max, tree->Ltree, locus, pop, chr);
  if(max >= tree->event->site)
    tree->Rtree = freeHistoryNodeRange(min, max, tree->Rtree, locus, pop, chr);
  if(tree->event->site >= min && tree->event->site <= max){
    /* remove reference to numCopies */
    for(i=0; i<tree->event->maxHaps[pop]; i++){
      if(tree->event->hapFreq[pop][i] == &popn[pop].numCopies[chr][locus]){
	tree->event->hapFreq[pop][i] = NULL;
	if(tree->event->nextHap[pop] > i)
	  tree->event->nextHap[pop] = i;
	tree->event->numCarriers[pop] -= popn[pop].numCopies[chr][locus];
	break;
      }
    }
    if(tree->Rtree == NULL && tree->Ltree == NULL){
      tree->event = NULL;
      free(tree);
      return NULL;
    }
    else{
      tmp = popTrash();
      tmp->Rtree = tree->Parent;
      tree = Delete_node_hist(tree, locus);
      if(tree != NULL)
	tree->Parent = tmp->Rtree;
      tmp->Rtree = NULL;
      free(tmp);
    }
  }
  return tree;
}  /* freeHistoryNodeRange */

/* ------------------------------------------------------------------------- */

void adjustCarriers(int pop)
{
  long i, j, k, r;

  for(r=0; r<gpars.R; r++){
    for(i=0; i<mutAr[r]->numMuts; i++){
      for(j=0; j<gpars.NPOP; j++){
	mutAr[r]->muts[i]->event->numCarriers[j] = 0;
	for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[j]; k++){
	  if(mutAr[r]->muts[i]->event->hapFreq[j][k] != NULL && 
	     *(mutAr[r]->muts[i]->event->hapFreq[j][k]) > 0){
	    mutAr[r]->muts[i]->event->numCarriers[j] +=
	      *(mutAr[r]->muts[i]->event->hapFreq[j][k]);
	  }
	}
      }
    }
  }
}

/* ------------------------------------------------------------------------- */
/* this function randomly chooses a haplotype to search for substitutions
   across all loci (a substitution is carried by every chromosome, randomly
   choosing one to look at ensures that intermediate frequency mutations
   are also updated).
   initially call this function with pop, chr=-1, locus, tree=NULL.
*/
void checkForFixed(const int pop, long chr, long locus, struct history *tree,
		   int inc, int BURNING)
{
  int D = (gpars.P == 4 ? 2 : 1);/*hap/diploid males have 1 X, tet have 2Xs*/
  long i;
  if(chr == -1){
    /* choose first segregating haplotype */
    for(chr=0; chr<popn[pop].maxchr; chr++){
      if(popn[pop].numCopies[chr][locus] == 0){
	if(chr < popn[pop].nextchr[locus]){
	  popn[pop].nextchr[locus] = chr;
	}
      }
      else{
	break;
      }
    }
    checkForFixed(pop, chr, locus, popn[pop].BigHead[chr][locus]->Rtree, inc, BURNING);
  }
  else{
    if(tree == NULL)
      return;
    /* recurse tree in post-order fashion to maintain balance */
    checkForFixed(pop, chr, locus, tree->Ltree, inc, BURNING);
    checkForFixed(pop, chr, locus, tree->Rtree, inc, BURNING);

    if(tree->event != NULL && tree->event->free == '0' && 
       tree->event->fixed[pop] == '0' &&
       (inc!=0 || tree->event->checkGen != popn[pop].gen ||
	tree->event->checkPop != pop)){
      tree->event->numCarriers[pop] = 0;
      for(i=0; i<tree->event->maxHaps[pop]; i++){
	if(tree->event->hapFreq[pop][i] != NULL){
	  if((*tree->event->hapFreq[pop][i]) == 0){
	    tree->event->hapFreq[pop][i] = NULL;
	    if(tree->event->nextHap[pop] > i)
	      tree->event->nextHap[pop] = i;
	  }
	  else
	    tree->event->numCarriers[pop] += (*tree->event->hapFreq[pop][i]);
	}
      }
    }
    if((tree->event->axy == 'A' &&
	tree->event->numCarriers[pop] == gpars.P*ppars[pop].N) ||
       (tree->event->axy == 'X' &&        /* #female Xs + #male Xs */
	tree->event->numCarriers[pop] ==
	(gpars.P*ppars[pop].MALES+(ppars[pop].N-ppars[pop].MALES)*D)) ||
       (tree->event->axy == 'Y' &&        /* just #male Ys */
	tree->event->numCarriers[pop] ==
	D*(ppars[pop].N-ppars[pop].MALES))){
#ifdef UNIT_TEST
      fprintf(errfile,"REMOVING SITE %ld POP %d LOCUS %ld freq=%ld/%ld (%d)!!!\n",
	     tree->event->site, pop, locus, tree->event->numCarriers[pop],
	     gpars.P*ppars[pop].N, 
	     popn[pop].polySites[locus][tree->event->site]);
#endif
      if(!BURNING || gpars.KEEPFIXED){
	tree->event->fixed[pop] = '1';
	if(tree->event->nonsyn == '0' || tree->event->nonsyn == '1'){
	  /* update consensus sequence only if substitution */
	  popn[pop].conSeq[locus][tree->event->site] = tree->event->derNuc;

	}
      }
      else{
	tree->event->free = '1';
      }
      tree->event->genFix[pop] = popn[pop].gen;
      
      if(tree->event->nonsyn == '0' || tree->event->nonsyn == '1'){
	/* update polySites */
	if(popn[pop].polySites[locus][tree->event->site]>0){
	  popn[pop].polySites[locus][tree->event->site]--;
#ifdef VERBOSE_DEBUG
	  fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!!\n",pop, locus,
		  tree->event->site,
		  popn[pop].polySites[locus][tree->event->site]);
	  fflush(errfile);
#endif
	}
#ifdef UNIT_TEST
#ifdef VERBOSE_DEBUG
	else{
	  fprintf(errfile,"ps[%d][%ld][%ld] would -> <0 (%d)!!!\n",pop, locus,
		  tree->event->site,
		  popn[pop].polySites[locus][tree->event->site]);
	  fflush(errfile);
	}
#endif
	if(popn[pop].polySites[locus][tree->event->site] < 0){
	  fprintf(errfile,"ps[%d][%ld][%ld] = %d < 0!\n", pop, locus, tree->event->site,
		  popn[pop].polySites[locus][tree->event->site]);
	  abort();
	}
#endif
      }

      /* now zero numCarriers */
      tree->event->numCarriers[pop]= 0;
      tree->event->nextHap[pop] = 0;
      for(i=0; i<tree->event->maxHaps[pop]; i++)
	tree->event->hapFreq[pop][i] = NULL;
      if(!BURNING || gpars.KEEPFIXED)
	cpyToFixed(&tree, &popn[pop].fixed[locus]->Rtree, pop);
      
      /* free from all haplotypes & update fitness */
      for(i=0; i<=popn[pop].maxchr; i++){
	if(fabs(tree->event->fit) > DBL_MIN){
	  if(gpars.ADDITIVE)
	    popn[pop].indfit[i][locus] -= tree->event->fit;
	  else
	    popn[pop].indfit[i][locus] /= (1.0 + tree->event->fit);
	}
	if(i != chr)
	  popn[pop].BigHead[i][locus]->Rtree =
	    freeHistoryNode(tree->event, popn[pop].BigHead[i][locus]->Rtree,
			    locus);
      }
      
#ifdef VERBOSE_DEBUG
      fprintf(errfile,"before free (%c:%d):\n", tree->event->fixed[pop],
	     popn[pop].polySites[locus][tree->event->site]);
      fflush(errfile);
      i=0;
      PrintHistTree(popn[pop].BigHead[chr][locus]->Rtree, &i, pop);
#endif
      /* now free from chr */
      popn[pop].BigHead[chr][locus]->Rtree =
	freeHistoryNode(tree->event, popn[pop].BigHead[chr][locus]->Rtree,
			locus);
#ifdef VERBOSE_DEBUG
      fprintf(errfile,"after free:\n");
      fflush(errfile);
      i=0;
      PrintHistTree(popn[pop].BigHead[chr][locus]->Rtree, &i, pop);
#endif      
    }
  }
}

/* ------------------------------------------------------------------------- */

void freeDead(const int p, const int step, int BURNING)
{
  long i, r;
  
  while(dead->next != NULL){
    i = dead->next->chr;
    r = dead->next->loc;
    if(popn[p].numCopies[i][r]==0 && popn[p].BigHead[i][r]->Rtree!=NULL){
#ifdef UNIT_TEST
      long j;
      for(j=0; j<ppars[p].N; j++){
	if(popn[p].parentLoc[j][r] == i){
	  fprintf(errfile,"popn[%d].parentLoc[%ld][%ld]=%ld, but numCopies == 0!!\n",p,j,r,i);
	  fflush(errfile);
	  abort();
	}
      }
#endif
#ifdef VERBOSE_DEBUG
      long foo=0;
      fprintf(errfile,"(freeDead) Freeing BigHead[%d][%ld][%ld]\n",p,i,r);
      PrintHistTree(popn[p].BigHead[i][r]->Rtree, &foo,p);
      fprintf(errfile,"\n");
      fflush(errfile);
#endif
      freeHistory(popn[p].BigHead[i][r]->Rtree, p, i, r, step, 1);
      if(i<popn[p].nextchr[r])
	popn[p].nextchr[r] = i;
      pushTrash(&popn[p].BigHead[i][r]->Rtree);
      popn[p].BigHead[i][r]->Rtree = NULL;
    }
    if(dead->next->next == NULL){
      free(dead->next);
      dead->next = NULL;
      break;
    }
    else{
      dead->next = dead->next->next;
      free(dead->next->prev);
      dead->next->prev = NULL;
    }
  }
  for(r=0; r<gpars.R; r++)
    checkForFixed(p, -1, r, NULL, 0, BURNING);
}

/* ------------------------------------------------------------------------- */

void freeTrash(struct history *tree)
{
  if(tree == NULL)
    return;
  else{
    if(tree->Parent != NULL){
      if(tree->Parent->Ltree == tree) /* we don't want to recursively return to this node */
	tree->Parent->Ltree = NULL;
      if(tree->Parent->Rtree == tree)
	tree->Parent->Rtree = NULL;
      freeTrash(tree->Parent);
    }
    if(tree->Rtree != NULL){
      tree->Rtree->Parent = NULL;
      freeTrash(tree->Rtree);
    }
    if(tree->Ltree != NULL){
      tree->Ltree->Parent = NULL;
      freeTrash(tree->Ltree);
    }
    tree->Rtree = NULL;
    tree->Ltree = NULL;
    tree->Parent = NULL;
    free(tree);
    tree = NULL;
  }
}

/* ------------------------------------------------------------------------- */

void cpyToFixed(struct history **tree, struct history **fixed, const int pop)
{
  if((*fixed) == NULL){ /* reached leaf */
    (*fixed) = popTrash();
    (*fixed)->event = (*tree)->event;
  }
  else if((*fixed)->event->site < (*tree)->event->site)
    cpyToFixed(tree, &(*fixed)->Rtree, pop);
  else if((*fixed)->event->site > (*tree)->event->site)
    cpyToFixed(tree, &(*fixed)->Ltree, pop);
  else if((*fixed)->event->gen < (*tree)->event->gen)
    cpyToFixed(tree, &(*fixed)->Rtree, pop);
  else if((*fixed)->event->gen > (*tree)->event->gen)
    cpyToFixed(tree, &(*fixed)->Ltree, pop);
  else if((*fixed)->event->genFix[pop] < (*tree)->event->genFix[pop])
    cpyToFixed(tree, &(*fixed)->Rtree, pop);
  else if((*fixed)->event->genFix[pop] > (*tree)->event->genFix[pop])
    cpyToFixed(tree, &(*fixed)->Rtree, pop);
  else{
    fprintf(errfile,"error in cpyToFixed, mutation already exists!\n");
    abort();
  }
}

/* ------------------------------------------------------------------------- */

void freeFixedHistory(struct history *tree, const int pop)
{
  int i;
  if(tree != NULL){
    freeFixedHistory(tree->Ltree, pop);
    freeFixedHistory(tree->Rtree, pop);
    if(tree->event != NULL){
      tree->event->free = '1';
      for(i=0; i<gpars.NPOP; i++)
	tree->event->fixed[i] = '0';
    }
    tree->event = NULL;
    tree->Rtree = NULL;
    tree->Ltree = NULL;
    tree->Parent = NULL;
    free(tree);
  }
}

/* ------------------------------------------------------------------------- */

void clearMutArray()
{
  long i, j, k, r;
  for(r=0; r<gpars.R; r++){
    for(i=0; i<mutAr[r]->numMuts; i++){
      mutAr[r]->muts[i]->event->free = '1';
      for(j=0; j<gpars.NPOP; j++){
	mutAr[r]->muts[i]->event->genFix[j] = 0;
	mutAr[r]->muts[i]->event->genDead[j] = -LONG_MAX;
	mutAr[r]->muts[i]->event->fixed[j] = '0';
	mutAr[r]->muts[i]->event->nextHap[j] = 0;
	mutAr[r]->muts[i]->event->numCarriers[j] = 0;
	for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[j]; k++)
	  mutAr[r]->muts[i]->event->hapFreq[j][k] = NULL;
      }
    }
    mutAr[r]->mutIndex = 0;
  }
}

/* ------------------------------------------------------------------------- */

void freeMutArray()
{
  long i, j, k, m;
  for(i=0; i<gpars.R; i++){
    for(j=0; j<mutAr[i]->numMuts; j++){
      for(k=0; k<gpars.NPOP; k++){
	for(m=0; m<mutAr[i]->muts[j]->event->maxHaps[k]; m++)
	  mutAr[i]->muts[j]->event->hapFreq[k][m] = NULL;
	free(mutAr[i]->muts[j]->event->hapFreq[k]);
      }
      free(mutAr[i]->muts[j]->event->genFix);
      free(mutAr[i]->muts[j]->event->genDead);
      free(mutAr[i]->muts[j]->event->fixed);
      free(mutAr[i]->muts[j]->event->nextHap);
      free(mutAr[i]->muts[j]->event->numCarriers);
      free(mutAr[i]->muts[j]->event->maxHaps);
      free(mutAr[i]->muts[j]->event->hapFreq);
      if(mutAr[i]->muts[j]->event->nucs != NULL)
	free(mutAr[i]->muts[j]->event->nucs);
      free(mutAr[i]->muts[j]->event);
      free(mutAr[i]->muts[j]);
    }
    free(mutAr[i]->muts);
    free(mutAr[i]);
  }
  free(mutAr);
  mutAr = NULL;
}

/* ------------------------------------------------------------------------- */

void genNewHitProbTot(double **hpTOT,double **hp,char **conSeq, 
		      long popSize, const long locus)
{
  int j;
  long i, k, r;
  char c0,c1,c2;
  
  for(k=0;k<popSize;k++)
    for(r=0;r<gpars.R;r++){
      hpTOT[k][r]=0.0;
      c0=conSeq[r][0];
      c1=conSeq[r][1];
      for(i=1;i<gpars.L[locus]-1;i++){
	c2=conSeq[r][i+1];
	hp[r][i]=0.0;
	for(j = 0;j < 4;j++)
	  if(j != c1-'0'){
	    hpTOT[k][r] += Q[16*(c0-'0')+4*(c1-'0')+(c2-'0')][j];
	    hp[r][i] += Q[16*(c0-'0')+4*(c1-'0')+(c2-'0')][j];
	  }
	c0=c1;
	c1=c2;
      }
    }
}

/* ------------------------------------------------------------------------- */

void exitNOW(const char *s)
{
  fprintf(errfile,"%s\nSEED=%ld",s,INITIALSEED);
  fflush(errfile);
  abort();
}

/* ------------------------------------------------------------------------- */

void updateMutFreqs()
{
  int n, pop, D = (gpars.P == 4 ? 2 : 1);
  long i, j, k, r;
  for(r=0; r<gpars.R; r++){
    for(i=0; i<mutAr[r]->numMuts; i++){
      if(mutAr[r]->muts[i]->event->free == '1')
	continue; /* nothing to do */
      n = 0; /* number of pops mutation is not segregating in */
      for(pop=0; pop<gpars.NPOP; pop++){
	if(!ppars[pop].ALIVE){
	  n++;
	  continue;
	} 
	else if(mutAr[r]->muts[i]->event->fixed[pop] == '1')
	  continue; /* nothing to count */
	
	mutAr[r]->muts[i]->event->numCarriers[pop] = 0;

	for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
	  if(mutAr[r]->muts[i]->event->hapFreq[pop][j] != NULL &&
	     (*mutAr[r]->muts[i]->event->hapFreq[pop][j]) > 0){
	    mutAr[r]->muts[i]->event->numCarriers[pop] += 
	      (*mutAr[r]->muts[i]->event->hapFreq[pop][j]);
	  }
	  else{ /* dead chromosome, remove it */
	      for(k=0; k<=popn[pop].maxchr; k++){
		if(mutAr[r]->muts[i]->event->hapFreq[pop][j] ==
		   &popn[pop].numCopies[k][r]){
#ifdef VERBOSE_DEBUG
		  fprintf(errfile,"(updateMutFreqs) Freeing BigHead[%d][%ld][%ld]\n",pop,k,r);
		  long m = 0;
		  PrintHistTree(popn[pop].BigHead[k][r]->Rtree, &m, pop);
		  fprintf(errfile,"\n");
		  fflush(errfile);
#endif
		  if(k<popn[pop].nextchr[r])
		    popn[pop].nextchr[r] = k;
		  pushTrash(&popn[pop].BigHead[k][r]->Rtree);
		  popn[pop].BigHead[k][r]->Rtree = NULL;
		  break;
		}
	      }
	    mutAr[r]->muts[i]->event->hapFreq[pop][j] = NULL;
	    if(j < mutAr[r]->muts[i]->event->nextHap[pop])
	      mutAr[r]->muts[i]->event->nextHap[pop] = j;
	  }
	}

	if(mutAr[r]->muts[i]->event->numCarriers[pop] == 0){/* mutation lost */
	  if(mutAr[r]->muts[i]->event->nonsyn == '0' || 
	     mutAr[r]->muts[i]->event->nonsyn == '1'){
#ifdef UNIT_TEST
	    fprintf(errfile,"REMOVING DEAD SITE %ld LOCUS %ld freq=%ld/%ld (%d)!!!\n",
		   mutAr[r]->muts[i]->event->site, r, 
		   mutAr[r]->muts[i]->event->numCarriers[pop],
		   gpars.P*ppars[pop].N, 
		   popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
#endif
	    if(popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]>0){
	      if(mutAr[r]->muts[i]->event->genDead[pop] != popn[pop].gen){
		popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]--;
		mutAr[r]->muts[i]->event->genDead[pop] = popn[pop].gen;
#ifdef VERBOSE_DEBUG
		fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!!\n",pop, r,
			mutAr[r]->muts[i]->event->site,
			popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
		fflush(errfile);
#endif
	      }
#ifdef VERBOSE_DEBUG
	      else{
		fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!! ALREADY\n",pop, r,
			mutAr[r]->muts[i]->event->site,
			popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
		fflush(errfile);
	      }
#endif
	    }
#ifdef VERBOSE_DEBUG
	    else{
	      fprintf(errfile,"ps[%d][%ld][%ld] would -> <0 (%d)!!!\n",pop, r,
		     mutAr[r]->muts[i]->event->site,
		      popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
	      fflush(errfile);
	    }
#endif

	  }
	  n++;
	  if(i<mutAr[r]->mutIndex)
	    mutAr[r]->mutIndex = i;
	}
	else if((mutAr[r]->muts[i]->event->axy == 'A' &&
		 mutAr[r]->muts[i]->event->numCarriers[pop] ==
		 gpars.P*ppars[pop].N) ||
		(mutAr[r]->muts[i]->event->axy == 'X' && /* #fem X + #male X */
		 mutAr[r]->muts[i]->event->numCarriers[pop] ==
		 (gpars.P*ppars[pop].MALES+(ppars[pop].N-ppars[pop].MALES)*D))
		|| (mutAr[r]->muts[i]->event->axy == 'Y' && /* just #male Ys */
		    mutAr[r]->muts[i]->event->numCarriers[pop] ==
		    D*(ppars[pop].N-ppars[pop].MALES))){

#ifdef UNIT_TEST
	  fprintf(errfile,"REMOVING FIXED SITE %ld LOCUS %ld freq=%ld/%ld (%d)!!!\n",
		 mutAr[r]->muts[i]->event->site, r, 
		 mutAr[r]->muts[i]->event->numCarriers[pop],
		 gpars.P*ppars[pop].N, 
		 popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
#endif
	  mutAr[r]->muts[i]->event->fixed[pop] = '1';
	  mutAr[r]->muts[i]->event->genFix[pop] = popn[pop].gen;
	  
	  if(mutAr[r]->muts[i]->event->nonsyn == '0' || 
	     mutAr[r]->muts[i]->event->nonsyn == '1'){
	    /* update consensus sequence only if substitution */
	    popn[pop].conSeq[r][mutAr[r]->muts[i]->event->site] =
	      mutAr[r]->muts[i]->event->derNuc;
	    /* update polySites */
	    if(popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]>0){
	      popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]--;
#ifdef VERBOSE_DEBUG
	      fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!!\n",pop, r,
		      mutAr[r]->muts[i]->event->site,
		      popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
	      fflush(errfile);
#endif
	    }
#ifdef UNIT_TEST
#ifdef VERBOSE_DEBUG
	    else{
	      fprintf(errfile,"ps[%d][%ld][%ld] would -> <0 (%d)!!!\n",pop, r,
		     mutAr[r]->muts[i]->event->site,
		      popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
	      fflush(errfile);
	    }
#endif
	    if(popn[pop].polySites[r][mutAr[r]->muts[i]->event->site] < 0){
	      fprintf(errfile,"ps[%d][%ld][%ld] = %d < 0!\n", pop, r,
		     mutAr[r]->muts[i]->event->site,
		     popn[pop].polySites[r]
		     [mutAr[r]->muts[i]->event->site]);
	      abort();
	    }
#endif
	  }   

	  /* now zero numCarriers */
	  mutAr[r]->muts[i]->event->numCarriers[pop] = 0;
	  mutAr[r]->muts[i]->event->nextHap[pop] = 0;
	  for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[pop]; k++)
	    mutAr[r]->muts[i]->event->hapFreq[pop][k] = NULL;
	  cpyToFixed(&mutAr[r]->muts[i], &popn[pop].fixed[r]->Rtree, pop);
	  
	  /* free from all haplotypes & update fitness */
	  for(k=0; k<=popn[pop].maxchr; k++){
	    if(fabs(mutAr[r]->muts[i]->event->fit) > DBL_MIN){
	      if(gpars.ADDITIVE)
		popn[pop].indfit[k][r] -= mutAr[r]->muts[i]->event->fit;
	      else
		popn[pop].indfit[k][r] /= (1.0 + mutAr[r]->muts[i]->event->fit);
	    }
	    popn[pop].BigHead[k][r]->Rtree =
	      freeHistoryNode(mutAr[r]->muts[i]->event,
			      popn[pop].BigHead[k][r]->Rtree, r);
	  }
	}
      }
      if(n == gpars.NPOP) /* free in all populations! */
	mutAr[r]->muts[i]->event->free = '1';
    }
  }
}

/* ------------------------------------------------------------------------- */

void copyFixed(struct history *tree, int pop, long chr, long locus)
{
  if(tree == NULL)
    return;
  copyFixed(tree->Rtree, pop, chr, locus);
  copyFixed(tree->Ltree, pop, chr, locus);
  
  if(tree->event->fixed[pop] == '0'){ /* not fixed, must copy */
    mutAr[locus]->mutIndex = tree->event->index;
    tree->event->genDead[pop] = -LONG_MAX;
    if(tree->event->numCarriers[pop] == 0 &&
       (tree->event->nonsyn == '0' || tree->event->nonsyn == '1')){
      popn[pop].polySites[locus][tree->event->site]++;
      tree->event->genDead[pop] = -LONG_MAX;
    }
#ifdef VERBOSE_DEBUG
    fprintf(errfile,"update PS[%ld]=%d\n",tree->event->site, 
	   popn[pop].polySites[locus][tree->event->site]);
    if(popn[pop].polySites[locus][tree->event->site]>1){
      //getchar();
    }
#endif    
/*     tree->event->numCarriers[pop]++; */
    addHistoryNode(tree->event, &popn[pop].BigHead[chr][locus]->Rtree, NULL,
		   pop, chr, locus, 1);
  }
}

/* ------------------------------------------------------------------------- */

void copyNchrs(long femalesIN, long malesIN, const int popTo, const int popFrom,
	       long femaleToStart, long maleToStart)
{
  int p, tmpIN;
  long j, t, rem=0;

#ifdef UNIT_TEST
  long foo;
#endif  
  
  /* copy females first */
  tmpIN = femalesIN;
  if(femaleToStart >= 0){
    rem = femaleToStart-gpars.P;
  }
  while(tmpIN > 0){
    while((t = ((long)(ran1(&gpars.seed)*ppars[popFrom].MALES))*gpars.P) == 
	  gpars.P*ppars[popFrom].MALES);
    
    /* now must replace individuals in popTo */
    if(femaleToStart >= 0)
      rem += gpars.P;
    else /* choose replaced individual randomly */
      rem = ((long)(ran1(&gpars.seed)*ppars[popTo].MALES))*gpars.P;

    if(rem == gpars.P*ppars[popTo].MALES)  rem -= 2;
#ifdef UNIT_TEST
    fprintf(errfile,"MIGRATION(F):  %d:%ld(%ld) -> %d:%ld(%ld)\n",popFrom,t,
	   popn[popFrom].parentLoc[t][0], popTo,rem,
	   popn[popTo].parentLoc[rem][0]);
    fflush(errfile);
#endif
    for(j=0; j<gpars.R; j++){
      for(p=0; p<gpars.P; p++){
#ifdef UNIT_TEST
	{
	  fprintf(errfile,"migrant = %ld:%ld -> ",rem+p,popn[popTo].parentLoc[rem+p][j]);
	  fflush(errfile);
	  foo=0;
	  PrintHistTree(popn[popFrom].BigHead[popn[popFrom].parentLoc
					    [t+p][j]][j]->Rtree, &foo, popFrom);
	  fflush(errfile);
	}
#endif
	createNewHistory(rem+p, j, popTo, 0, 1);
#ifdef UNIT_TEST
	{
	  fprintf(errfile,"%ld:%ld (%ld)\n",rem+p,popn[popTo].parentLoc[rem+p][j],
		 popn[popTo].numCopies[popn[popTo].parentLoc[rem+p][j]][j]);
	  if(popn[popTo].BigHead[popn[popTo].parentLoc[rem+p][j]][j]->Rtree !=
	     NULL){
	    fprintf(errfile,"error, didn't clear bighead in migrate!\n");
	    abort();
	  }
	}
#endif
	if(gpars.TRACKANC){
	  memcpy(popn[popTo].ancestry[popn[popTo].parentLoc[rem+p][j]][j],
		 popn[popFrom].ancestry[popn[popFrom].parentLoc[t+p][j]][j],
		 gpars.L[j]*sizeof(int));
	}
	if(popn[popFrom].BigHead[popn[popFrom].parentLoc[t+p][j]][j]->Rtree !=
	   NULL)
	  copyHistory(&popn[popTo].BigHead[popn[popTo].parentLoc[rem+p][j]][j],
		      &popn[popFrom].BigHead[popn[popFrom].parentLoc[t+p][j]][j]
		      , j, popTo, popn[popTo].parentLoc[rem+p][j], popFrom, 1, 
		      0, 0);
	copyFixed(popn[popFrom].fixed[j]->Rtree, popTo, 
		  popn[popTo].parentLoc[rem+p][j], j);
	checkStillFixed(popn[popTo].fixed[j]->Rtree, popTo, popFrom, j,
			popn[popTo].parentLoc[rem+p][j]);
#ifdef UNIT_TEST
	foo=0;
	PrintHistTree(popn[popTo].BigHead[popn[popTo].parentLoc[rem+p][j]]
		      [j]->Rtree,&foo,popTo);
	fflush(errfile);
#endif

	if(gpars.substMod == 5)
	  memcpy(popn[popTo].hpSUM[popn[popTo].parentLoc[rem+p][j]][j],
		 popn[popFrom].hpSUM[popn[popFrom].parentLoc[t+p][j]][j],
		 gpars.L[j]*sizeof(double));
	popn[popTo].indfit[popn[popTo].parentLoc[rem+p][j]][j] = 
	  popn[popFrom].indfit[popn[popFrom].parentLoc[t+p][j]][j];
	popn[popTo].extremeMuts[popn[popTo].parentLoc[rem+p][j]][j][0] =
	  popn[popFrom].extremeMuts[popn[popFrom].parentLoc[t+p][j]][j][0];
	popn[popTo].extremeMuts[popn[popTo].parentLoc[rem+p][j]][j][1] =
	  popn[popFrom].extremeMuts[popn[popFrom].parentLoc[t+p][j]][j][1];
      }
    }
    tmpIN--;
  }

  /* copy males */
  tmpIN = malesIN;
  if(maleToStart >= 0){
    rem = maleToStart - gpars.P;
  }
  while(tmpIN > 0){
    while((t = gpars.P*((long)(ppars[popFrom].N*
			       (1.0-ran1(&gpars.seed)*
				(1.0-ppars[popFrom].pFEMALES))))) >=
	  gpars.P*ppars[popFrom].N-1 || t < gpars.P*ppars[popFrom].MALES);
#ifdef UNIT_TEST
    fprintf(errfile,"MIGRATION(M):  %d:%ld -> ",popFrom,t);
    fflush(errfile);
#endif
    /* now must replace individuals in popTo */
    if(maleToStart >= 0)
      rem += gpars.P;
    else
      while((rem = (gpars.P*((long)(ppars[popTo].N*
				   (1.0-ran1(&gpars.seed)*
				    (1.0-ppars[popTo].pFEMALES)))))) >=
	    gpars.P*ppars[popTo].N-1 || rem < gpars.P*ppars[popTo].MALES){
#ifdef UNIT_TEST
	fprintf(errfile,"rejecting %d:%ld [%ld, %ld]\n",popTo,rem, gpars.P*ppars[popTo].N-1,
	       gpars.P*ppars[popTo].MALES);
	fflush(errfile);
#endif
      }
#ifdef UNIT_TEST
    fprintf(errfile,"%d:%ld\n",popTo,rem);
    fflush(errfile);
#endif
    for(j=0; j<gpars.R; j++){
      for(p=0; p<gpars.P; p++){
#ifdef UNIT_TEST
	fprintf(errfile,"migrant = %ld:%ld -> ",rem+p,popn[popTo].parentLoc[rem+p][j]);
	fflush(errfile);
	foo=0;
	PrintHistTree(popn[popFrom].BigHead[popn[popFrom].parentLoc
					  [t+p][j]][j]->Rtree, &foo, popFrom);
	fflush(errfile);
#endif
	createNewHistory(rem+p, j, popTo, 0, 1);
	if(gpars.TRACKANC){
	  memcpy(popn[popTo].ancestry[popn[popTo].parentLoc[rem+p][j]][j],
		 popn[popFrom].ancestry[popn[popFrom].parentLoc[t+p][j]][j],
		 gpars.L[j]*sizeof(int));
	}
	if(popn[popFrom].BigHead[popn[popFrom].parentLoc[t+p][j]][j]->Rtree
	   != NULL)
	  copyHistory(&popn[popTo].BigHead[popn[popTo].parentLoc[rem+p][j]][j],
		      &popn[popFrom].BigHead[popn[popFrom].parentLoc[t+p][j]][j]
		      , j, popTo, popn[popTo].parentLoc[rem+p][j], popFrom, 1,
		      0, 0);
	copyFixed(popn[popFrom].fixed[j]->Rtree, popTo, 
		  popn[popTo].parentLoc[rem+p][j], j);
	checkStillFixed(popn[popTo].fixed[j]->Rtree, popTo, popFrom, j,
			popn[popTo].parentLoc[rem+p][j]);
#ifdef UNIT_TEST
	fprintf(errfile,"%ld:%ld (%ld)\n",rem+p,popn[popTo].parentLoc[rem+p][j],
	       popn[popTo].numCopies[popn[popTo].parentLoc[rem+p][j]][j]);
	fprintf(errfile,"%ld:%ld",rem+p,popn[popTo].parentLoc[rem+p][j]);
	fflush(errfile);
	foo=0;
	PrintHistTree(popn[popTo].BigHead[popn[popTo].parentLoc
					  [rem+p][j]][j]->Rtree, &foo, popTo);
	fflush(errfile);
	fprintf(errfile,"\n");
	fprintf(errfile,"Fixed Tree:\n");
	foo=0;
	PrintHistTree(popn[popTo].fixed[j]->Rtree, &foo, popTo);
#endif
	if(gpars.substMod == 5)
	  memcpy(popn[popTo].hpSUM[popn[popTo].parentLoc[rem+p][j]][j],
		 popn[popFrom].hpSUM[popn[popFrom].parentLoc[t+p][j]][j],
		 gpars.L[j]*sizeof(double));
	popn[popTo].indfit[popn[popTo].parentLoc[rem+p][j]][j] = 
	  popn[popFrom].indfit[popn[popFrom].parentLoc[t+p][j]][j];
	popn[popTo].extremeMuts[popn[popTo].parentLoc[rem+p][j]][j][0] =
	  popn[popFrom].extremeMuts[popn[popFrom].parentLoc[t+p][j]][j][0];
	popn[popTo].extremeMuts[popn[popTo].parentLoc[rem+p][j]][j][1] =
	  popn[popFrom].extremeMuts[popn[popFrom].parentLoc[t+p][j]][j][1];
      }
    }
    tmpIN--;
  }
}

/* ------------------------------------------------------------------------- */

void checkStillFixed(struct history *tree, int popTo, int popFrom, long locus,
		     long chr)
{
  long i;
  int foo;
  if(tree == NULL)
    return;
  checkStillFixed(tree->Ltree, popTo, popFrom, locus, chr);
  checkStillFixed(tree->Rtree, popTo, popFrom, locus, chr);
  
  if(tree->event->fixed[popFrom] == '1') /* fixed in both, nothing to do */
    return;
  foo = 0;
  carryMut(popn[popTo].BigHead[chr][locus]->Rtree, tree->event->index, &foo);
  if(foo == 1){
    tree->event->numCarriers[popTo] = 0;
    for(i = 0; i<tree->event->maxHaps[popTo]; i++){
      if(tree->event->hapFreq[popTo][i] != NULL){
	if(tree->event->hapFreq[popTo][i]==&popn[popTo].numCopies[chr][locus] ||
	   (*tree->event->hapFreq[popTo][i]) == 0){
	  tree->event->hapFreq[popTo][i] = NULL;
	  if(tree->event->nextHap[popTo] > i)
	    tree->event->nextHap[popTo] = i;
	}
	else
	  tree->event->numCarriers[popTo] += (*tree->event->hapFreq[popTo][i]);
      }
    }
    if(tree->event->numCarriers[popTo] == 0 && (tree->event->nonsyn == '0' ||
						tree->event->nonsyn == '1')){
      if(popn[popTo].polySites[locus][tree->event->site]>0){
	if(tree->event->genDead[popTo] != popn[popTo].gen){
	  popn[popTo].polySites[locus][tree->event->site]--;
	  tree->event->genDead[popTo] = popn[popTo].gen;
#ifdef VERBOSE_DEBUG
	  fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!!\n",popTo, locus,
		  tree->event->site,
		  popn[popTo].polySites[locus][tree->event->site]);
	  fflush(errfile);
#endif
	}
#ifdef VERBOSE_DEBUG
	else{
	  fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!! ALREADY\n",popTo, locus,
		  tree->event->site,
		  popn[popTo].polySites[locus][tree->event->site]);
	  fflush(errfile);
	}
#endif
      }
#ifdef UNIT_TEST
#ifdef VERBOSE_DEBUG
      else{
	fprintf(errfile,"ps[%d][%ld][%ld] would -> <0 (%d)!!!\n",popTo, locus,
		tree->event->site,
		popn[popTo].polySites[locus][tree->event->site]);
	fflush(errfile);
      }
#endif
      if(popn[popTo].polySites[locus][tree->event->site] < 0){
	fprintf(errfile,"ps[%d][%ld][%ld] = %d!!!\n",popTo, locus,
		tree->event->site,
		popn[popTo].polySites[locus][tree->event->site]);
	fflush(errfile);
	abort();
      }
#endif
    }
    popn[popTo].BigHead[chr][locus]->Rtree =
      freeHistoryNode(tree->event, popn[popTo].BigHead[chr][locus]->Rtree,
		      locus);
  }
  else{ /* mutation no longer fixed in population */
    tree->event->fixed[popTo] = '0';
    popn[popTo].polySites[locus][tree->event->site]++;
    tree->event->genDead[popTo] = -LONG_MAX;
#ifdef VERBOSE_DEBUG
    fprintf(errfile,"update PS[%ld]=%d\n",tree->event->site, 
	   popn[popTo].polySites[locus][tree->event->site]);
    if(popn[popTo].polySites[locus][tree->event->site]>1){
      //getchar();
    }
#endif    
    tree->event->numCarriers[popTo] = 0;
    mutAr[locus]->mutIndex = tree->event->index;
    for(i=0; i<=popn[popTo].maxchr; i++){
      if(i == chr || popn[popTo].numCopies[i][locus] == 0)
	continue;
      foo=0;
      carryMut(popn[popTo].BigHead[i][locus], tree->event->index, &foo);
      addHistoryNode(tree->event, &popn[popTo].BigHead[i][locus]->Rtree, NULL,
		     popTo, i, locus, 1);
      
    }
    popn[popTo].fixed[locus]->Rtree = 
      freeHistoryNode(tree->event, popn[popTo].fixed[locus]->Rtree, locus);
  }
}

/* ------------------------------------------------------------------------- */

void getMinSS(struct history *head, int *minSS)
{
  if(head == NULL)
    *minSS = -1;
  else if(head->Ltree == NULL)
    *minSS = head->event->site;
  else
    getMinSS(head->Ltree, minSS);
}

/* ------------------------------------------------------------------------- */

void checkHistTreeOrder(struct history *head, long min, long max, int pop)
{
  if(head != NULL){
    if(head->event != NULL && head->event->free == '1'){
      long i=0;
      fprintf(errfile,"error in checkHistTreeOrder, head carries free mut!\n");
      PrintHistTree(head, &i, pop);
      abort();
    }
    if(head->Ltree != NULL){
      if(head->Ltree->event->site < min || 
	 head->Ltree->event->site > head->event->site){
	long i=0;
	fprintf(errfile,"error in checkHistTreeOrder (L:%ld,%ld):\n", min, max);
	PrintHistTree(head, &i, pop);
	abort();
      }
      checkHistTreeOrder(head->Ltree, min, head->event->site, pop);
    }
    if(head->Rtree != NULL){
      if(head->Rtree->event->site < head->event->site ||
	 head->Rtree->event->site > max){
	long i=0;
	fprintf(errfile,"error in checkHistTreeOrder (R:%ld,%ld):\n", min, max);
	PrintHistTree(head, &i, pop);
	abort();
      }
      checkHistTreeOrder(head->Rtree, head->event->site, max, pop);
    }
  }
}

/* ------------------------------------------------------------------------- */

void freeEvolEvents(struct EvolEvents *devents)
{
  if(devents == NULL)
    return;
  else{
    freeEvolEvents(devents->nextEvent);
    free(devents);
    devents = NULL;
  }
}

/* ------------------------------------------------------------------------- */

void addEvolEvent(struct EvolEvents tmp, struct EvolEvents **devents, 
		  struct EvolEvents **parent)
{
  int i;
  long j;
  if((*devents) != NULL && 
     (tmp.tau > (*devents)->tau+DBL_MIN || 
      (fabs(tmp.tau - (*devents)->tau) <= DBL_MIN && tmp.eventType!=9 &&
       ((*devents)->eventType==9 ||
	(tmp.eventType < 2 && (*devents)->eventType < 2 &&
	 tmp.popj > (*devents)->popj) ||
	tmp.eventType > (*devents)->eventType)))){
    addEvolEvent(tmp, &(*devents)->nextEvent, devents);
  }
  else{
    struct EvolEvents *new;
    assert(new = malloc(sizeof(*new)));
    new->popi = -1;
    new->popj = -1;
    new->nAncPops = -1;
    new->ancPops = NULL;
    new->maleFreqs = NULL;
    new->femaleFreqs = NULL;
    new->nu = -1;
    new->freq = -1;
    new->parIndex = -1;
    new->SS = NULL;
    new->newP.selDistType = NULL;
    new->newP.GAMMA = NULL;
    new->newP.ProbPos = NULL;
    new->newP.ProbDel = NULL;
    new->newP.ProbNeut = NULL;
    new->newP.ProbNegGamma = NULL;
    new->newP.ProbPosGamma = NULL;
    new->newP.alphaN = NULL;
    new->newP.lambdaN = NULL;
    new->newP.alphaP = NULL;
    new->newP.lambdaP = NULL;
    new->newP.normMean = NULL;
    new->newP.normVar = NULL;
    new->newP.f0 = NULL;

    new->newP.codomPop = NULL;
    new->newP.hDistType = NULL;
    new->newP.h = NULL;
    new->newP.hProbDom = NULL;
    new->newP.hProbGamma = NULL;
    new->newP.halpha = NULL;
    new->newP.hlambda = NULL;
    new->newP.hnormMean = NULL;
    new->newP.hnormVar = NULL;
    
    new->newG.mig_mat = NULL;
    new->eventType = tmp.eventType;
    new->tau = tmp.tau;
    new->locus = tmp.locus;
    new->site = tmp.site;
    if(tmp.eventType == 0){
      new->popi = tmp.popi;
      new->popj = tmp.popj;
    }
    else if(tmp.eventType == 1){
      new->popi = tmp.popi;
      new->popj = tmp.popj;
      new->nAncPops = tmp.nAncPops;
      if(new->nAncPops > 0){
	assert(new->ancPops = malloc(new->nAncPops*sizeof(int)));
	assert(new->maleFreqs = malloc(new->nAncPops*sizeof(float)));
	assert(new->femaleFreqs = malloc(new->nAncPops*sizeof(float)));
	for(i=0; i<new->nAncPops; i++){
	  new->ancPops[i] = tmp.ancPops[i];
	  new->maleFreqs[i] = tmp.maleFreqs[i];
	  new->femaleFreqs[i] = tmp.femaleFreqs[i];
	}
      }
      new->freq = tmp.freq;
      new->newP.Nt = tmp.newP.Nt;
    }
    else if(tmp.eventType == 3){
      new->popi = tmp.popi;
      new->nu = tmp.nu;
    }
    else if(tmp.eventType == 5){
      new->popi = tmp.popi;
    }
    else if(tmp.eventType == 6){
      new->popi = tmp.popi;
      new->parIndex = tmp.parIndex;
      if(new->parIndex == 0)
	new->newP.Nt = tmp.newP.Nt;
      else if(new->parIndex == 1)
	new->newP.popAlpha = tmp.newP.popAlpha;
      else if(new->parIndex == 2){
	new->newP.popAlpha = tmp.newP.popAlpha;
	new->newP.K = tmp.newP.K;
      }
      else if(new->parIndex == 3){
	new->newP.THETA = tmp.newP.THETA;
      }
      else if(new->parIndex == 4){
	new->newP.RHO = tmp.newP.RHO;
      }
      else if(new->parIndex == 5){
	new->newP.SELF = tmp.newP.SELF;
      }
      /* parIndex = 6 included in (7) below */
      else if(new->parIndex == 7){
	if(new->newP.selDistType == NULL)
	  assert(new->newP.selDistType = malloc(gpars.R*sizeof(int)));
	if(new->newP.codomPop == NULL)
	  assert(new->newP.codomPop = malloc(gpars.R*sizeof(int)));
	if(new->newP.hDistType == NULL)
	  assert(new->newP.hDistType = malloc(gpars.R*sizeof(int)));
	j=(tmp.locus == -1 ? 0 : tmp.locus);
	for(;(tmp.locus==-1&&j<gpars.R)||j==tmp.locus; j++){
	  new->newP.selDistType[j] = tmp.newP.selDistType[j];
	  switch(new->newP.selDistType[j]){
	  case 0:
	    break;
	  case 1:
	    if(new->newP.GAMMA == NULL){
	      assert(new->newP.GAMMA = malloc(gpars.R*sizeof(*new->newP.GAMMA)));
	      assert(new->newP.ProbPos =
		     malloc(gpars.R*sizeof(*new->newP.ProbPos)));
	      assert(new->newP.ProbDel =
		     malloc(gpars.R*sizeof(*new->newP.ProbDel)));
	      assert(new->newP.ProbNeut =
		     malloc(gpars.R*sizeof(*new->newP.ProbNeut)));
	    }
	    new->newP.GAMMA[j] = tmp.newP.GAMMA[j];
	    new->newP.ProbPos[j] = tmp.newP.ProbPos[j];
	    new->newP.ProbDel[j] = tmp.newP.ProbDel[j];
	    new->newP.ProbNeut[j] = tmp.newP.ProbNeut[j];
	    break;
	  case 2:
	    if(new->newP.ProbNegGamma == NULL){
	      assert(new->newP.ProbNegGamma = 
		     malloc(gpars.R*sizeof(*new->newP.ProbNegGamma)));
	      assert(new->newP.ProbPosGamma = 
		     malloc(gpars.R*sizeof(*new->newP.ProbPosGamma)));
	      assert(new->newP.alphaN = 
		     malloc(gpars.R*sizeof(*new->newP.alphaN)));
	      assert(new->newP.alphaP = 
		     malloc(gpars.R*sizeof(*new->newP.alphaP)));
	      assert(new->newP.lambdaN = 
		     malloc(gpars.R*sizeof(*new->newP.lambdaN)));
	      assert(new->newP.lambdaP = 
		     malloc(gpars.R*sizeof(*new->newP.lambdaP)));
	    }
	    new->newP.ProbNegGamma[j] = tmp.newP.ProbNegGamma[j];
	    new->newP.ProbPosGamma[j] = tmp.newP.ProbPosGamma[j];
	    new->newP.alphaN[j] = tmp.newP.alphaN[j];
	    new->newP.alphaP[j] = tmp.newP.alphaP[j];
	    new->newP.lambdaN[j] = tmp.newP.lambdaN[j];
	    new->newP.lambdaP[j] = tmp.newP.lambdaP[j];
	    break;
	  case 3:
	    if(new->newP.normMean == NULL){
	      assert(new->newP.normMean =
		     malloc(gpars.R*sizeof(*new->newP.normMean)));
	      assert(new->newP.normVar =
		     malloc(gpars.R*sizeof(*new->newP.normVar)));
	    }
	    new->newP.normMean[j] = tmp.newP.normMean[j];
	    new->newP.normVar[j] = tmp.newP.normVar[j];
	    break;
	  case 4:
	    break;
	  default:
	    fprintf(errfile,"sfs_code error (%ld): bad selDistType = %d\n",
		    INITIALSEED, new->newP.selDistType[j]);
	    abort();
	    break;
	  }

	  new->newP.codomPop[j] = tmp.newP.codomPop[j];
	  new->newP.hDistType[j] = tmp.newP.hDistType[j];
	  switch(new->newP.hDistType[j]){
	  case 0:
	    break;
	  case 1:
	    if(new->newP.h == NULL){
	      assert(new->newP.h =malloc(gpars.R*sizeof(*new->newP.h)));
	      assert(new->newP.hProbDom =
		     malloc(gpars.R*sizeof(*new->newP.hProbDom)));
	    }
	    new->newP.h[j] = tmp.newP.h[j];
	    new->newP.hProbDom[j] = tmp.newP.hProbDom[j];
	    break;
	  case 2:
	    if(new->newP.hProbGamma == NULL){
	      assert(new->newP.hProbGamma = 
		     malloc(gpars.R*sizeof(*new->newP.hProbGamma)));
	      assert(new->newP.halpha = 
		     malloc(gpars.R*sizeof(*new->newP.halpha)));
	      assert(new->newP.hlambda = 
		     malloc(gpars.R*sizeof(*new->newP.hlambda)));
	    }
	    new->newP.hProbGamma[j] = tmp.newP.hProbGamma[j];
	    new->newP.halpha[j] = tmp.newP.halpha[j];
	    new->newP.hlambda[j] = tmp.newP.hlambda[j];
	    break;
	  case 3:
	    if(new->newP.hnormMean == NULL){
	      assert(new->newP.hnormMean =
		     malloc(gpars.R*sizeof(*new->newP.hnormMean)));
	      assert(new->newP.hnormVar =
		     malloc(gpars.R*sizeof(*new->newP.hnormVar)));
	    }
	    new->newP.hnormMean[j] = tmp.newP.hnormMean[j];
	    new->newP.hnormVar[j] = tmp.newP.hnormVar[j];
	    break;
	  default:
	    fprintf(errfile,"sfs_code error (%ld): bad hDistType = %d\n",
		    INITIALSEED, new->newP.hDistType[j]);
	    abort();
	    break;
	  }
	}
      }
      else if(new->parIndex == 8){
	new->newP.KAPPA = tmp.newP.KAPPA;
      }
      else if(new->parIndex == 9){
	new->newP.PSI = tmp.newP.PSI;
      }
      else if(new->parIndex == 10){
	assert(new->newP.f0 = malloc(gpars.R*sizeof(*new->newP.f0)));
	if(tmp.locus == -1){
	  for(j=0; j<gpars.R; j++)
	    new->newP.f0[j] = tmp.newP.f0[j];
	}
	else
	  new->newP.f0[tmp.locus] = tmp.newP.f0[tmp.locus];
      }
      else if(new->parIndex == 11){
	new->newP.pMaleMig = tmp.newP.pMaleMig;
      }
      else if(new->parIndex == 12){
	new->newP.GenEffect = tmp.newP.GenEffect;
      }
      else if(new->parIndex == 13){
	new->newP.INSRATE = tmp.newP.INSRATE;
	new->newP.DELRATE = tmp.newP.DELRATE;
	new->newP.INDELlength = tmp.newP.INDELlength;
      }
      else if(new->parIndex == 14){
	new->newP.longINSRATE = tmp.newP.longINSRATE;
	new->newP.longDELRATE = tmp.newP.longDELRATE;
	new->newP.longINDELlength = tmp.newP.longINDELlength;
      }
      else if(new->parIndex == 15){
	new->newP.INVRATE = tmp.newP.INVRATE;
	new->newP.INVlength = tmp.newP.INVlength;
      }
      else if(new->parIndex == 16){
	new->newP.BGC = tmp.newP.BGC;
	new->newP.fGC = tmp.newP.fGC;
	new->newP.GCtract = tmp.newP.GCtract;
      }
      else if(new->parIndex == 17){
	new->newP.pMaleRec = tmp.newP.pMaleRec;
      }
      else if(new->parIndex == 18){
	for(i=0; i<4; i++)
	  new->newP.baseFreq[i] = tmp.newP.baseFreq[i];
      }
      else if(new->parIndex == 19){ /* sampling time */
	new->nu = tmp.nu;
	new->freq = tmp.freq;
	assert(new->SS = malloc(gpars.NPOP*sizeof(*new->SS)));
	for(i=0; i<gpars.NPOP; i++){
	  new->SS[i] = tmp.SS[i];
	}
      }
      else if(new->parIndex == 20){ /* switch selection coefficient */
	j=(tmp.locus == -1 ? 0 : tmp.locus);
	new->freq = tmp.freq;
	new->nu = tmp.nu;
	if(new->newP.selDistType == NULL)
	  assert(new->newP.selDistType = malloc(gpars.R*sizeof(int)));
	if(new->newP.codomPop == NULL)
	  assert(new->newP.codomPop = malloc(gpars.R*sizeof(int)));
	if(new->newP.hDistType == NULL)
	  assert(new->newP.hDistType = malloc(gpars.R*sizeof(int)));
	for(;(tmp.locus==-1&&j<gpars.R)||j==tmp.locus; j++){
	  new->newP.selDistType[j] = tmp.newP.selDistType[j];
	  new->newP.hDistType[j] = tmp.newP.hDistType[j];
	  new->newP.codomPop[j] = tmp.newP.codomPop[j];
	  switch(new->newP.selDistType[j]){
	  case 0:
	    break;
	  case 1:
	    if(new->newP.GAMMA == NULL){
	      assert(new->newP.GAMMA =malloc(gpars.R*sizeof(*new->newP.GAMMA)));
	      assert(new->newP.ProbPos =
		     malloc(gpars.R*sizeof(*new->newP.ProbPos)));
	      assert(new->newP.ProbDel =
		     malloc(gpars.R*sizeof(*new->newP.ProbDel)));
	      assert(new->newP.ProbNeut =
		     malloc(gpars.R*sizeof(*new->newP.ProbNeut)));
	    }
	    new->newP.GAMMA[j] = tmp.newP.GAMMA[j];
	    new->newP.ProbPos[j] = tmp.newP.ProbPos[j];
	    new->newP.ProbDel[j] = tmp.newP.ProbDel[j];
	    new->newP.ProbNeut[j] = tmp.newP.ProbNeut[j];
	    break;
	  case 2:
	    if(new->newP.ProbNegGamma == NULL){
	      assert(new->newP.ProbNegGamma = 
		     malloc(gpars.R*sizeof(*new->newP.ProbNegGamma)));
	      assert(new->newP.ProbPosGamma = 
		     malloc(gpars.R*sizeof(*new->newP.ProbPosGamma)));
	      assert(new->newP.alphaN = 
		     malloc(gpars.R*sizeof(*new->newP.alphaN)));
	      assert(new->newP.alphaP = 
		     malloc(gpars.R*sizeof(*new->newP.alphaP)));
	      assert(new->newP.lambdaN = 
		     malloc(gpars.R*sizeof(*new->newP.lambdaN)));
	      assert(new->newP.lambdaP = 
		     malloc(gpars.R*sizeof(*new->newP.lambdaP)));
	    }
	    new->newP.ProbNegGamma[j] = tmp.newP.ProbNegGamma[j];
	    new->newP.ProbPosGamma[j] = tmp.newP.ProbPosGamma[j];
	    new->newP.alphaN[j] = tmp.newP.alphaN[j];
	    new->newP.alphaP[j] = tmp.newP.alphaP[j];
	    new->newP.lambdaN[j] = tmp.newP.lambdaN[j];
	    new->newP.lambdaP[j] = tmp.newP.lambdaP[j];
	    break;
	  case 3:
	    if(new->newP.normMean == NULL){
	      assert(new->newP.normMean =
		     malloc(gpars.R*sizeof(*new->newP.normMean)));
	      assert(new->newP.normVar =
		     malloc(gpars.R*sizeof(*new->newP.normVar)));
	    }
	    new->newP.normMean[j] = tmp.newP.normMean[j];
	    new->newP.normVar[j] = tmp.newP.normVar[j];
	    break;
	  case 4:
	    break;
	  default:
	    fprintf(errfile,"sfs_code error (%ld): bad selDistType = %d\n",
		    INITIALSEED, new->newP.selDistType[j]);
	    abort();
	    break;
	  }

	  new->newP.codomPop[j] = tmp.newP.codomPop[j];
	  new->newP.hDistType[j] = tmp.newP.hDistType[j];
	  switch(new->newP.hDistType[j]){
	  case 0:
	    break;
	  case 1:
	    if(new->newP.h == NULL){
	      assert(new->newP.h =malloc(gpars.R*sizeof(*new->newP.h)));
	      assert(new->newP.hProbDom =
		     malloc(gpars.R*sizeof(*new->newP.hProbDom)));
	    }
	    new->newP.h[j] = tmp.newP.h[j];
	    new->newP.hProbDom[j] = tmp.newP.hProbDom[j];
	    break;
	  case 2:
	    if(new->newP.hProbGamma == NULL){
	      assert(new->newP.hProbGamma = 
		     malloc(gpars.R*sizeof(*new->newP.hProbGamma)));
	      assert(new->newP.halpha = 
		     malloc(gpars.R*sizeof(*new->newP.halpha)));
	      assert(new->newP.hlambda = 
		     malloc(gpars.R*sizeof(*new->newP.hlambda)));
	    }
	    new->newP.hProbGamma[j] = tmp.newP.hProbGamma[j];
	    new->newP.halpha[j] = tmp.newP.halpha[j];
	    new->newP.hlambda[j] = tmp.newP.hlambda[j];
	    break;
	  case 3:
	    if(new->newP.hnormMean == NULL){
	      assert(new->newP.hnormMean =
		     malloc(gpars.R*sizeof(*new->newP.hnormMean)));
	      assert(new->newP.hnormVar =
		     malloc(gpars.R*sizeof(*new->newP.hnormVar)));
	    }
	    new->newP.hnormMean[j] = tmp.newP.hnormMean[j];
	    new->newP.hnormVar[j] = tmp.newP.hnormVar[j];
	    break;
	  default:
	    fprintf(errfile,"sfs_code error (%ld): bad hDistType = %d\n",
		    INITIALSEED, new->newP.hDistType[j]);
	    abort();
	    break;
	  }
	}
      }
      else{
	fprintf(errfile,"unfortunately, parIndex %d has not been implemented: 1\n",
	       new->parIndex);
	abort();
      }
    }
    else if(tmp.eventType == 7){
      new->parIndex = tmp.parIndex;
      if(tmp.parIndex == 0){
	new->popi = tmp.popi;
	new->popj = tmp.popj;
	assert(new->newG.mig_mat =
	       malloc(gpars.NPOP*sizeof(*new->newG.mig_mat)));
	for(i=0; i<gpars.NPOP; i++){
	  assert(new->newG.mig_mat[i] =
		 malloc(gpars.NPOP*sizeof(**new->newG.mig_mat)));
	  for(j=0; j<gpars.NPOP; j++)
	    new->newG.mig_mat[i][j] = tmp.newG.mig_mat[i][j];
	}
      }
      else if(tmp.parIndex == 1){
	new->popi = tmp.popi;
	new->locus = tmp.locus;
	new->site = tmp.site;
      }
    }
    else if(tmp.eventType == 8){ /* admixture */
      new->popi = tmp.popi;
      new->newP.Nt = tmp.newP.Nt;
      new->nAncPops = tmp.nAncPops;
      assert(new->ancPops = malloc(new->nAncPops*sizeof(int)));
      assert(new->maleFreqs = malloc(new->nAncPops*sizeof(float)));
      assert(new->femaleFreqs = malloc(new->nAncPops*sizeof(float)));
      for(i=0; i<new->nAncPops; i++){
	  new->ancPops[i] = tmp.ancPops[i];
	  new->maleFreqs[i] = tmp.maleFreqs[i];
	  new->femaleFreqs[i] = tmp.femaleFreqs[i];
      }
    }
    else if(tmp.eventType == 9){
      new->popi = tmp.popi;
      new->locus = tmp.locus;
      new->site = tmp.site;
      new->gamma = tmp.gamma;
      new->h = tmp.h;
    }
    else{
      fprintf(errfile,"sfs_code error:  unfortunately, event type %d has not \
been implemented yet...\n",tmp.eventType);
      abort();
    }
    if((*devents) == NULL && parent == NULL){
      new->nextEvent = NULL;
      (*devents) = new;
    }
    else if((*devents) == NULL){
      new->nextEvent = NULL;
      (*devents) = new;
      (*parent)->nextEvent = (*devents);
    }
    else{
      struct EvolEvents *tmp2;
      assert(tmp2 = malloc(sizeof(*tmp2)));
      tmp2->popi = -1;
      tmp2->popj = -1;
      tmp2->nAncPops = -1;
      tmp2->ancPops = NULL;
      tmp2->maleFreqs = NULL;
      tmp2->femaleFreqs = NULL;
      tmp2->nu = -1;
      tmp2->freq = -1;
      tmp2->parIndex = -1;
      tmp2->SS = NULL;
      tmp2->newG.mig_mat = NULL;
      tmp2->newP.selDistType = NULL;
      tmp2->newP.hDistType = NULL;
      tmp2->newP.codomPop = NULL;
      tmp2->newP.GAMMA = NULL;
      tmp2->newP.ProbPos = NULL;
      tmp2->newP.ProbDel = NULL;
      tmp2->newP.ProbNeut = NULL;
      tmp2->newP.ProbNegGamma = NULL;
      tmp2->newP.ProbPosGamma = NULL;
      tmp2->newP.alphaN = NULL;
      tmp2->newP.lambdaN = NULL;
      tmp2->newP.alphaP = NULL;
      tmp2->newP.lambdaP = NULL;
      tmp2->newP.normMean = NULL;
      tmp2->newP.normVar = NULL;
      tmp2->newP.f0 = NULL;
      tmp2->newP.hDistType = NULL;
      tmp2->newP.h = NULL;
      tmp2->newP.hProbDom = NULL;
      tmp2->newP.hProbGamma = NULL;
      tmp2->newP.halpha = NULL;
      tmp2->newP.hlambda = NULL;
      tmp2->newP.hnormMean = NULL;
      tmp2->newP.hnormVar = NULL;
      tmp2->eventType = (*devents)->eventType;
      tmp2->tau = (*devents)->tau;
      tmp2->locus = (*devents)->locus;
      tmp2->site = (*devents)->site;
      if((*devents)->eventType == 0){
	tmp2->popi = (*devents)->popi;
	tmp2->popj = (*devents)->popj;
      }
      else if((*devents)->eventType == 1){
	tmp2->popi = (*devents)->popi;
	tmp2->popj = (*devents)->popj;
	tmp2->freq = (*devents)->freq;
	tmp2->newP.Nt = (*devents)->newP.Nt;
      }
      else if((*devents)->eventType == 7){
	tmp2->parIndex = (*devents)->parIndex;
	tmp2->popi = (*devents)->popi;
	if((*devents)->parIndex == 0){
	  tmp2->popj = (*devents)->popj;
	  assert(tmp2->newG.mig_mat = malloc(gpars.NPOP*sizeof(double *)));
	  for(i=0; i<gpars.NPOP; i++){
	    assert(tmp2->newG.mig_mat[i] =
		   malloc(gpars.NPOP*sizeof(double)));
	    memcpy(tmp2->newG.mig_mat[i], (*devents)->newG.mig_mat[i],
		   gpars.NPOP*sizeof(double));
	  }
	}
	else if((*devents)->parIndex == 1){
	  tmp2->locus = (*devents)->locus;
	  tmp2->site = (*devents)->site;
	}
      }
      else if((*devents)->eventType == 8){
	tmp2->popi = (*devents)->popi;
	tmp2->newP.Nt = (*devents)->newP.Nt;
	tmp2->nAncPops = (*devents)->nAncPops;
	assert(tmp2->ancPops = malloc(tmp2->nAncPops*sizeof(int)));
	assert(tmp2->maleFreqs = malloc(tmp2->nAncPops*sizeof(float)));
	assert(tmp2->femaleFreqs = malloc(tmp2->nAncPops*sizeof(float)));
	for(i=0; i<tmp2->nAncPops; i++){
	  tmp2->ancPops[i] = (*devents)->ancPops[i];
	  tmp2->maleFreqs[i] = (*devents)->maleFreqs[i];
	  tmp2->femaleFreqs[i] = (*devents)->femaleFreqs[i];
	}
      }
      else if((*devents)->eventType == 9){
	tmp2->popi = (*devents)->popi;
	tmp2->locus = (*devents)->locus;
	tmp2->site = (*devents)->site;
	tmp2->gamma = (*devents)->gamma;
      }
      else if((*devents)->eventType == 3){
	tmp2->popi = (*devents)->popi;
	tmp2->nu = (*devents)->nu;
      }
      else if((*devents)->eventType == 5){
	tmp2->popi = (*devents)->popi;
      }
      else if((*devents)->eventType == 6){
	tmp2->popi = (*devents)->popi;
	tmp2->parIndex = (*devents)->parIndex;
	if(tmp2->parIndex == 0)
	  tmp2->newP.Nt = (*devents)->newP.Nt;
	else if(tmp2->parIndex == 1)
	  tmp2->newP.popAlpha = (*devents)->newP.popAlpha;
	else if(tmp2->parIndex == 2){
	  tmp2->newP.popAlpha = (*devents)->newP.popAlpha;
	  tmp2->newP.K = (*devents)->newP.K;
	}
	else if(tmp2->parIndex == 3){
	  tmp2->newP.THETA = (*devents)->newP.THETA;
	}
	else if(tmp2->parIndex == 4){
	  tmp2->newP.RHO = (*devents)->newP.RHO;
	}
	else if(tmp2->parIndex == 5){
	  tmp2->newP.SELF = (*devents)->newP.SELF;
	}
	/* parIndex = 6 included in (7) below */
	else if(tmp2->parIndex == 7){
	  j=((*devents)->locus == -1 ? 0 : (*devents)->locus);
	  if(tmp2->newP.selDistType == NULL)
	    assert(tmp2->newP.selDistType = malloc(gpars.R*sizeof(int)));
	  if(tmp2->newP.codomPop == NULL)
	    assert(tmp2->newP.codomPop = malloc(gpars.R*sizeof(int)));
	  if(tmp2->newP.hDistType == NULL)
	    assert(tmp2->newP.hDistType = malloc(gpars.R*sizeof(int)));
	  for(;((*devents)->locus==-1&&j<gpars.R)||j==(*devents)->locus; j++){
	    tmp2->newP.selDistType[j] = (*devents)->newP.selDistType[j];
	    tmp2->newP.hDistType[j] = (*devents)->newP.hDistType[j];
	    tmp2->newP.codomPop[j] = (*devents)->newP.codomPop[j];
	    switch(tmp2->newP.selDistType[j]){
	    case 0:
	      break;
	    case 1:
	      if(tmp2->newP.GAMMA == NULL){
		assert(tmp2->newP.GAMMA =
		       malloc(gpars.R*sizeof(*tmp2->newP.GAMMA)));
		assert(tmp2->newP.ProbPos =
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbPos)));
		assert(tmp2->newP.ProbDel =
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbDel)));
		assert(tmp2->newP.ProbNeut =
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbNeut)));
	      }
	      tmp2->newP.GAMMA[j] = (*devents)->newP.GAMMA[j];
	      tmp2->newP.ProbPos[j] = (*devents)->newP.ProbPos[j];
	      tmp2->newP.ProbDel[j] = (*devents)->newP.ProbDel[j];
	      tmp2->newP.ProbNeut[j] = (*devents)->newP.ProbNeut[j];
	      break;
	    case 2:
	      if(tmp2->newP.ProbNegGamma == NULL){
		assert(tmp2->newP.ProbNegGamma = 
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbNegGamma)));
		assert(tmp2->newP.ProbPosGamma = 
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbPosGamma)));
		assert(tmp2->newP.alphaN = 
		       malloc(gpars.R*sizeof(*tmp2->newP.alphaN)));
		assert(tmp2->newP.alphaP = 
		       malloc(gpars.R*sizeof(*tmp2->newP.alphaP)));
		assert(tmp2->newP.lambdaN = 
		       malloc(gpars.R*sizeof(*tmp2->newP.lambdaN)));
		assert(tmp2->newP.lambdaP = 
		       malloc(gpars.R*sizeof(*tmp2->newP.lambdaP)));
	      }
	      tmp2->newP.ProbNegGamma[j] = (*devents)->newP.ProbNegGamma[j];
	      tmp2->newP.ProbPosGamma[j] = (*devents)->newP.ProbPosGamma[j];
	      tmp2->newP.alphaN[j] = (*devents)->newP.alphaN[j];
	      tmp2->newP.alphaP[j] = (*devents)->newP.alphaP[j];
	      tmp2->newP.lambdaN[j] = (*devents)->newP.lambdaN[j];
	      tmp2->newP.lambdaP[j] = (*devents)->newP.lambdaP[j];
	      break;
	    case 3:
	      if(tmp2->newP.normMean == NULL){
		assert(tmp2->newP.normMean =
		       malloc(gpars.R*sizeof(*tmp2->newP.normMean)));
		assert(tmp2->newP.normVar =
		       malloc(gpars.R*sizeof(*tmp2->newP.normVar)));
	      }
	      tmp2->newP.normMean[j] = (*devents)->newP.normMean[j];
	      tmp2->newP.normVar[j] = (*devents)->newP.normVar[j];
	      break;
	    case 4:
	      break;
	    default:
	      fprintf(errfile,"sfs_code error (%ld): bad selDistType = %d\n",
		      INITIALSEED, tmp2->newP.selDistType[j]);
	      abort();
	      break;
	    }
	  
	    tmp2->newP.codomPop[j] = (*devents)->newP.codomPop[j];
	    tmp2->newP.hDistType[j] = (*devents)->newP.hDistType[j];
	    switch(tmp2->newP.hDistType[j]){
	    case 0:
	      break;
	    case 1:
	      if(tmp2->newP.h == NULL){
		assert(tmp2->newP.h =
		       malloc(gpars.R*sizeof(*tmp2->newP.h)));
		assert(tmp2->newP.hProbDom =
		       malloc(gpars.R*sizeof(*tmp2->newP.hProbDom)));
	      }
	      tmp2->newP.h[j] = (*devents)->newP.h[j];
	      tmp2->newP.hProbDom[j] = (*devents)->newP.hProbDom[j];
	      break;
	    case 2:
	      if(tmp2->newP.hProbGamma == NULL){
		assert(tmp2->newP.hProbGamma = 
		       malloc(gpars.R*sizeof(*tmp2->newP.hProbGamma)));
		assert(tmp2->newP.halpha = 
		       malloc(gpars.R*sizeof(*tmp2->newP.halpha)));
		assert(tmp2->newP.hlambda = 
		       malloc(gpars.R*sizeof(*tmp2->newP.hlambda)));
	      }
	      tmp2->newP.hProbGamma[j] = (*devents)->newP.hProbGamma[j];
	      tmp2->newP.halpha[j] = (*devents)->newP.halpha[j];
	      tmp2->newP.hlambda[j] = (*devents)->newP.hlambda[j];
	      break;
	    case 3:
	      if(tmp2->newP.hnormMean == NULL){
		assert(tmp2->newP.hnormMean =
		       malloc(gpars.R*sizeof(*tmp2->newP.hnormMean)));
		assert(tmp2->newP.hnormVar =
		       malloc(gpars.R*sizeof(*tmp2->newP.hnormVar)));
	      }
	      tmp2->newP.hnormMean[j] = (*devents)->newP.hnormMean[j];
	      tmp2->newP.hnormVar[j] = (*devents)->newP.hnormVar[j];
	      break;
	    default:
	      fprintf(errfile,"sfs_code error (%ld): bad hDistType = %d\n",
		      INITIALSEED, tmp2->newP.hDistType[j]);
	      abort();
	      break;
	    }
	  }
	}
	else if(tmp2->parIndex == 8){
	  tmp2->newP.KAPPA = (*devents)->newP.KAPPA;
	}
	else if(tmp2->parIndex == 9){
	  tmp2->newP.PSI = (*devents)->newP.PSI;
	}
	else if(tmp2->parIndex == 10){
	  assert(tmp2->newP.f0 = malloc(gpars.R*sizeof(*tmp2->newP.f0)));
	  if((*devents)->locus == -1){
	    for(j=0; j<gpars.R; j++)
	      tmp2->newP.f0[j] = (*devents)->newP.f0[j];
	  }
	  else
	    tmp2->newP.f0[(*devents)->locus] =
	      (*devents)->newP.f0[(*devents)->locus];
	}
	else if(tmp2->parIndex == 11){
	  tmp2->newP.pMaleMig = (*devents)->newP.pMaleMig;
	}
	else if(tmp2->parIndex == 12){
	  tmp2->newP.GenEffect = (*devents)->newP.GenEffect;
	}
	else if(tmp2->parIndex == 13){
	  tmp2->newP.INSRATE = (*devents)->newP.INSRATE;
	  tmp2->newP.DELRATE = (*devents)->newP.DELRATE;
	  tmp2->newP.INDELlength = (*devents)->newP.INDELlength;
	}
	else if(tmp2->parIndex == 14){
	  tmp2->newP.longINSRATE = (*devents)->newP.longINSRATE;
	  tmp2->newP.longDELRATE = (*devents)->newP.longDELRATE;
	  tmp2->newP.longINDELlength = (*devents)->newP.longINDELlength;
	}
	else if(tmp2->parIndex == 15){
	  tmp2->newP.INVRATE = (*devents)->newP.INVRATE;
	  tmp2->newP.INVlength = (*devents)->newP.INVlength;
	}
	else if(tmp2->parIndex == 16){
	  tmp2->newP.BGC = (*devents)->newP.BGC;
	  tmp2->newP.fGC = (*devents)->newP.fGC;
	  tmp2->newP.GCtract = (*devents)->newP.GCtract;
	}
	else if(tmp2->parIndex == 17){
	  tmp2->newP.pMaleRec = (*devents)->newP.pMaleRec;
	}
	else if(tmp2->parIndex == 18){
	  for(i=0; i<4; i++)
	    tmp2->newP.baseFreq[i] = (*devents)->newP.baseFreq[i];
	}
	else if(tmp2->parIndex == 19){ /* sampling time */
	  tmp2->nu = (*devents)->nu;
	  tmp2->freq = (*devents)->freq;
	  assert(tmp2->SS = malloc(gpars.NPOP*sizeof(int)));
	  for(i=0; i<gpars.NPOP; i++){
	    tmp2->SS[i] = (*devents)->SS[i];
	  }
	}
	else if(tmp2->parIndex == 20){ /* switch selection coefficient */
	  j=((*devents)->locus == -1 ? 0 : (*devents)->locus);
	   tmp2->freq = (*devents)->freq;
	   tmp2->nu = (*devents)->nu;
	   if(tmp2->newP.selDistType == NULL)
	    assert(tmp2->newP.selDistType = malloc(gpars.R*sizeof(int)));
	   if(tmp2->newP.codomPop == NULL)
	    assert(tmp2->newP.codomPop = malloc(gpars.R*sizeof(int)));
	   if(tmp2->newP.hDistType == NULL)
	    assert(tmp2->newP.hDistType = malloc(gpars.R*sizeof(int)));
	  for(;((*devents)->locus==-1&&j<gpars.R)||j==(*devents)->locus; j++){
	    tmp2->newP.selDistType[j] = (*devents)->newP.selDistType[j];
	    switch(tmp2->newP.selDistType[j]){
	    case 0:
	      break;
	    case 1:
	      if(tmp2->newP.GAMMA == NULL){
		assert(tmp2->newP.GAMMA =
		       malloc(gpars.R*sizeof(*tmp2->newP.GAMMA)));
		assert(tmp2->newP.ProbPos =
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbPos)));
		assert(tmp2->newP.ProbDel =
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbDel)));
		assert(tmp2->newP.ProbNeut =
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbNeut)));
	      }
	      tmp2->newP.GAMMA[j] = (*devents)->newP.GAMMA[j];
	      tmp2->newP.ProbPos[j] = (*devents)->newP.ProbPos[j];
	      tmp2->newP.ProbDel[j] = (*devents)->newP.ProbDel[j];
	      tmp2->newP.ProbNeut[j] = (*devents)->newP.ProbNeut[j];
	      break;
	    case 2:
	      if(tmp2->newP.ProbNegGamma == NULL){
		assert(tmp2->newP.ProbNegGamma = 
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbNegGamma)));
		assert(tmp2->newP.ProbPosGamma = 
		       malloc(gpars.R*sizeof(*tmp2->newP.ProbPosGamma)));
		assert(tmp2->newP.alphaN = 
		       malloc(gpars.R*sizeof(*tmp2->newP.alphaN)));
		assert(tmp2->newP.alphaP = 
		       malloc(gpars.R*sizeof(*tmp2->newP.alphaP)));
		assert(tmp2->newP.lambdaN = 
		       malloc(gpars.R*sizeof(*tmp2->newP.lambdaN)));
		assert(tmp2->newP.lambdaP = 
		       malloc(gpars.R*sizeof(*tmp2->newP.lambdaP)));
	      }
	      tmp2->newP.ProbNegGamma[j] = (*devents)->newP.ProbNegGamma[j];
	      tmp2->newP.ProbPosGamma[j] = (*devents)->newP.ProbPosGamma[j];
	      tmp2->newP.alphaN[j] = (*devents)->newP.alphaN[j];
	      tmp2->newP.alphaP[j] = (*devents)->newP.alphaP[j];
	      tmp2->newP.lambdaN[j] = (*devents)->newP.lambdaN[j];
	      tmp2->newP.lambdaP[j] = (*devents)->newP.lambdaP[j];
	      break;
	    case 3:
	      if(tmp2->newP.normMean == NULL){
		assert(tmp2->newP.normMean =
		       malloc(gpars.R*sizeof(*tmp2->newP.normMean)));
		assert(tmp2->newP.normVar =
		       malloc(gpars.R*sizeof(*tmp2->newP.normVar)));
	      }
	      tmp2->newP.normMean[j] = (*devents)->newP.normMean[j];
	      tmp2->newP.normVar[j] = (*devents)->newP.normVar[j];
	      break;
	    case 4:
	      break;
	    default:
	      fprintf(errfile,"sfs_code error (%ld): bad selDistType = %d\n",
		      INITIALSEED, tmp2->newP.selDistType[j]);
	      abort();
	      break;
	    }
	  
	    tmp2->newP.codomPop[j] = (*devents)->newP.codomPop[j];
	    tmp2->newP.hDistType[j] = (*devents)->newP.hDistType[j];
	    switch(tmp2->newP.hDistType[j]){
	    case 0:
	      break;
	    case 1:
	      if(tmp2->newP.h == NULL){
		assert(tmp2->newP.h =
		       malloc(gpars.R*sizeof(*tmp2->newP.h)));
		assert(tmp2->newP.hProbDom =
		       malloc(gpars.R*sizeof(*tmp2->newP.hProbDom)));
	      }
	      tmp2->newP.h[j] = (*devents)->newP.h[j];
	      tmp2->newP.hProbDom[j] = (*devents)->newP.hProbDom[j];
	      break;
	    case 2:
	      if(tmp2->newP.hProbGamma == NULL){
		assert(tmp2->newP.hProbGamma = 
		       malloc(gpars.R*sizeof(*tmp2->newP.hProbGamma)));
		assert(tmp2->newP.halpha = 
		       malloc(gpars.R*sizeof(*tmp2->newP.halpha)));
		assert(tmp2->newP.hlambda = 
		       malloc(gpars.R*sizeof(*tmp2->newP.hlambda)));
	      }
	      tmp2->newP.hProbGamma[j] = (*devents)->newP.hProbGamma[j];
	      tmp2->newP.halpha[j] = (*devents)->newP.halpha[j];
	      tmp2->newP.hlambda[j] = (*devents)->newP.hlambda[j];
	      break;
	    case 3:
	      if(tmp2->newP.hnormMean == NULL){
		assert(tmp2->newP.hnormMean =
		       malloc(gpars.R*sizeof(*tmp2->newP.hnormMean)));
		assert(tmp2->newP.hnormVar =
		       malloc(gpars.R*sizeof(*tmp2->newP.hnormVar)));
	      }
	      tmp2->newP.hnormMean[j] = (*devents)->newP.hnormMean[j];
	      tmp2->newP.hnormVar[j] = (*devents)->newP.hnormVar[j];
	      break;
	    default:
	      fprintf(errfile,"sfs_code error (%ld): bad hDistType = %d\n",
		      INITIALSEED, tmp2->newP.hDistType[j]);
	      abort();
	      break;
	    }
	  }
	}
	else{
	  fprintf(errfile,"unfortunately, parIndex %d has not been implemented: 2\n",
		 tmp2->parIndex);
	  abort();
	}
      }
      
      tmp2->nextEvent = (*devents)->nextEvent;
      if((*devents)->newG.mig_mat != NULL){
	for(i=0; i<gpars.NPOP; i++)
	  free((*devents)->newG.mig_mat[i]);
	free((*devents)->newG.mig_mat);
	(*devents)->newG.mig_mat = NULL;
      }
      if((*devents)->newP.selDistType != NULL){
	if((*devents)->newP.GAMMA != NULL){
	  free((*devents)->newP.GAMMA);
	  (*devents)->newP.GAMMA = NULL;
	  free((*devents)->newP.ProbPos);
	  free((*devents)->newP.ProbDel);
	  free((*devents)->newP.ProbNeut);
	}
	else if((*devents)->newP.ProbNegGamma != NULL){
	  free((*devents)->newP.ProbNegGamma);
	  (*devents)->newP.ProbNegGamma = NULL;
	  free((*devents)->newP.ProbPosGamma);
	  (*devents)->newP.ProbPosGamma = NULL;
	  free((*devents)->newP.alphaP);
	  free((*devents)->newP.alphaN);
	  free((*devents)->newP.lambdaP);
	  free((*devents)->newP.lambdaN);
	}
	else if((*devents)->newP.normMean != NULL){
	  free((*devents)->newP.normMean);
	  (*devents)->newP.normMean = NULL;
	  free((*devents)->newP.normVar);
	}
	free((*devents)->newP.selDistType);
	(*devents)->newP.selDistType = NULL;
      }

      if((*devents)->newP.hDistType != NULL){
	if((*devents)->newP.h != NULL){
	  free((*devents)->newP.h);
	  (*devents)->newP.h = NULL;
	  free((*devents)->newP.hProbDom);
	}
	else if((*devents)->newP.hProbGamma != NULL){
	  free((*devents)->newP.hProbGamma);
	  (*devents)->newP.hProbGamma = NULL;
	  free((*devents)->newP.halpha);
	  free((*devents)->newP.hlambda);
	}
	else if((*devents)->newP.hnormMean != NULL){
	  free((*devents)->newP.hnormMean);
	  (*devents)->newP.hnormMean = NULL;
	  free((*devents)->newP.hnormVar);
	}
	free((*devents)->newP.codomPop);
	(*devents)->newP.codomPop = NULL;
	free((*devents)->newP.hDistType);
	(*devents)->newP.hDistType = NULL;
      }
      
      if((*devents)->ancPops != NULL){
	free((*devents)->ancPops);
	free((*devents)->maleFreqs);
	free((*devents)->femaleFreqs);
      }
      free(*devents);
      new->nextEvent = tmp2;
      (*devents) = new;
      if(parent != NULL) (*parent)->nextEvent = (*devents);
    }
    return;
  }
}

/* ------------------------------------------------------------------------- */

void checkEvolEventTimes(struct EvolEvents *list, double time, long PNanc)
{
  if(list == NULL)
    return;
  /*  if(list->tau - time < FLT_MIN)
      list->tau = time + 1.0/PNanc;*/
  if(list->nextEvent != NULL)
    checkEvolEventTimes(list->nextEvent, list->tau, PNanc);
  /* just to make sure different events occur in different generations */
}

/* ------------------------------------------------------------------------- */

void copyEvolEvents(struct EvolEvents **to, struct EvolEvents *from)
{
  if(from != NULL){
    addEvolEvent(*from, to, NULL);
  }
  else
    return;
}

/* ------------------------------------------------------------------------- */

void printEvolEvents(struct EvolEvents *tmp)
{
  if(tmp != NULL){
    printf("event type = %d\n",tmp->eventType);
    if(tmp->eventType == 6 || tmp->eventType == 7){
      printf("parIndex = %d\n", tmp->parIndex);
    }
    printf("\ttau = %f\n",tmp->tau);
    if(tmp->eventType == 6 && tmp->parIndex == 19){
      int i;
      printf("\tpopi = %d\n",tmp->popi);
      printf("\ttau_d = %f\n",tmp->freq);
      printf("\ttau_1 = %f\n",tmp->nu);
      printf("\tSS=");
      for(i=0; i<gpars.NPOP; i++){
	printf("%d:%ld; ", i, tmp->SS[i]);
      }
      printf("\n");
    }
    else{
      printf("\tpopi = %d\n",tmp->popi);
      printf("\tpopj = %d\n",tmp->popj);
      printf("\tnu = %f\n",tmp->nu);
      printf("\tsite = %ld\n",tmp->site);
      printf("\tfreq = %f\n",tmp->freq);
    }
    printEvolEvents(tmp->nextEvent);
  }
  else{
    printf("\n");
  }
}

/* ------------------------------------------------------------------------- */

void argcheck(int arg, int opt, int argc, char *argv[])
{
  if((arg >= argc ) || (argv[arg][0] == '-')){
    fprintf(errfile,"not enough arguments after %s (entry %d)\n", argv[opt],
	    arg) ;
    fprintf(errfile,"For usage type: sfs_code -h<return>\n");
    exit(0);
  }
}

/* ------------------------------------------------------------------------- */

void parseCommandSFSCODE(int argc, char *argv[], int it)
{
  int arg, ind, p1, p2, *popInit, toBeTracked=0, errPrinted=0;
  int pop, popSET, locSET, aBeforeL=-1;
  long i, j, k, loc;
  char tmpStr[100], args[100];
  struct EvolEvents tmpEvolEvent;
  tmpEvolEvent.locus = -1;
  tmpEvolEvent.site = -1;
  tmpEvolEvent.nAncPops = 0;
  tmpEvolEvent.ancPops = NULL;
  tmpEvolEvent.maleFreqs = NULL;
  tmpEvolEvent.femaleFreqs = NULL;
  tmpEvolEvent.SS = NULL;
  tmpEvolEvent.newP.selDistType = NULL;
  tmpEvolEvent.newP.GAMMA = NULL;
  tmpEvolEvent.newP.ProbPos = NULL;
  tmpEvolEvent.newP.ProbDel = NULL;
  tmpEvolEvent.newP.ProbNeut = NULL;
  tmpEvolEvent.newP.ProbNegGamma = NULL;
  tmpEvolEvent.newP.ProbPosGamma = NULL;
  tmpEvolEvent.newP.alphaN = NULL;
  tmpEvolEvent.newP.lambdaN = NULL;
  tmpEvolEvent.newP.alphaP = NULL;
  tmpEvolEvent.newP.lambdaP = NULL;
  tmpEvolEvent.newP.normMean = NULL;
  tmpEvolEvent.newP.normVar = NULL;
  tmpEvolEvent.newP.f0 = NULL;
  tmpEvolEvent.newP.codomPop = NULL;
  tmpEvolEvent.newP.hDistType = NULL;
  tmpEvolEvent.newP.h = NULL;
  tmpEvolEvent.newP.hProbDom = NULL;
  tmpEvolEvent.newP.hProbGamma = NULL;
  tmpEvolEvent.newP.halpha = NULL;
  tmpEvolEvent.newP.hlambda = NULL;
  tmpEvolEvent.newP.hnormMean = NULL;
  tmpEvolEvent.newP.hnormVar = NULL;
  tmpEvolEvent.newG.mig_mat = NULL;
  tmpEvolEvent.nextEvent = NULL;

  popInit = malloc(gpars.NPOP*sizeof(int));
  for(i=0; i<gpars.NPOP; i++)
    popInit[i] = 0;

  arg = 3;
  while(arg < argc){
    if(argv[arg][0] == '-'){
#ifdef WINDOWS      
      strcpy_s(args, 100, &argv[arg][1]);
#else
      strcpy(args, &argv[arg][1]);
#endif
      if(strcmp(args,"-help") == 0 || strcmp(args,"h") == 0){
	helpMenu();
	exit(0);
      }
      else if(strcmp(args,"-noSeq") == 0 || strcmp(args,"A") == 0){
	gpars.PRINTSEQ = 0;
	arg++;
      }
      else if(strcmp(args,"-printGen") == 0){
	gpars.PRINTGEN = 1;
	arg++;
      }
      else if(strcmp(args,"-printLocus") == 0){
	argcheck(arg+1, arg, argc, argv); /* check for 2 arguments */
	if(atoi(argv[arg+1]) >= gpars.R){
	  fprintf(errfile,"Error specifying locus to print.  Note that the %ld loci in this simulation are numbered from 0,...,%ld\n",gpars.R, gpars.R-1);
	  exit(1);
	}
	gpars.PRINTLOC[atoi(argv[arg+1])] = '1';
	arg+=2;
      }
      else if(strcmp(args, "-VCF") == 0){
	gpars.OUTVCF = 1;
	arg++;
      }
      else if(strcmp(args,"-outfile") == 0 || strcmp(args,"o") == 0){
	argcheck(arg+1, arg, argc, argv); /* check for 2 arguments */
	arg++;
	if(strcmp(argv[arg], "a") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* check for 2 arguments */
	  arg++;
#ifdef WINDOWS
	  strcpy_s(tmpStr, 100, argv[arg]); /* filename */
#else
	  strcpy(tmpStr, argv[arg]); /* filename */
#endif
	  if(outfile == stdout){
#ifdef WINDOWS
	    fopen_s(&outfile, tmpStr, "a"); /* w/a */
#else
	    if((outfile = fopen(tmpStr,"a"))==NULL){ /* w/a */
	      fprintf(errfile,"cannot open file %s\n",tmpStr);
	      abort();
	    }
#endif
	  }
	}
	else{
#ifdef WINDOWS
	  strcpy_s(tmpStr, 100, argv[arg]); /* filename */
#else
	  strcpy(tmpStr, argv[arg]); /* filename */
#endif
	  if(outfile == stdout)
#ifdef WINDOWS
	    fopen_s(&outfile, tmpStr, "w"); /* w/a */
#else
	  if((outfile = fopen(tmpStr,"w"))==NULL){ /* w/a */
	    fprintf(errfile,"cannot open file %s\n",tmpStr);
	    abort();
	  }
#endif
	}
	arg++;
      }
      else if(strcmp(args,"-errfile") == 0 || strcmp(args,"e") == 0){
	argcheck(arg+1, arg, argc, argv); /* check for 2 arguments */
	arg++;
	if(strcmp(argv[arg], "a") == 0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
#ifdef WINDOWS
	  strcpy_s(tmpStr, 100, argv[arg]); /* filename */
#else
	  strcpy(tmpStr, argv[arg]); /* filename */
#endif
	  if(errfile == stderr)
#ifdef WINDOWS
	    fopen_s(&errfile,tmpStr,"a");
#else
	  if((errfile = fopen(tmpStr,"a"))==NULL){ /* w/a */
	    fprintf(errfile,"cannot open file %s\n",tmpStr);
	    abort();
	  }
#endif
	}
	else{
#ifdef WINDOWS
	  strcpy_s(tmpStr, 100, argv[arg]);
#else
	  strcpy(tmpStr, argv[arg]);
#endif
	  if(errfile == stderr){
#ifdef WINDOWS
	    fopen_s(&errfile,tmpStr,"w");
#else
	    if((errfile = fopen(tmpStr,"w"))==NULL){ /* w/a */
	      fprintf(errfile,"cannot open file %s\n",tmpStr);
	      abort();
	    }
#endif
	  }
	}
	arg++;
      }
      else if(strcmp(args,"-popFreq") == 0 || strcmp(args, "F") == 0){
#ifdef WINDOWS
	fprintf(errfile,"option --popFreq (-F) not implemented on cluster, \
sorry!!\n");
	abort();
#else
	arg++;
	argcheck(arg, arg-1, argc, argv);
	/* get output file name */
	if(strcmp(argv[arg], "a") == 0){
	  arg++;
	  argcheck(arg, arg-2, argc, argv);
	  if(freqFile == NULL){
	    if((freqFile = fopen(argv[arg],"a"))==NULL){
	      fprintf(errfile,"cannot open file %s\n",argv[arg]);
	      abort();
	    }
	  }
	}
	else{
	  if(freqFile == NULL){
	    if((freqFile = fopen(argv[arg],"w"))==NULL){
	      fprintf(errfile,"cannot open file %s\n",argv[arg]);
	      abort();
	    }
	  }
	}
	arg++;
#endif
      }
      else if(strcmp(args,"-ReportBurnFixed") == 0 || strcmp(args,"R") == 0){
	arg++;
	gpars.KEEPFIXED = 1;
      }
      else if(strcmp(args,"-additive") == 0 || strcmp(args,"Z") == 0){
	gpars.ADDITIVE = 1;
	arg++;
      }
      else if(strcmp(args,"-selDistType") == 0 || strcmp(args,"W") == 0){
	int ODD = -1;
	argcheck(arg+1, arg, argc, argv);
	arg++;
	pop = 0;
	for(i=0; i<gpars.NPOP; i++){
	  if(ppars[i].stayNeut == 1){
	    pop++;
	  }
	  else{
	    break;
	  }
	}
	if(pop >= gpars.NPOP){
	  fprintf(errfile,"error:  requesting selection when all populations are neutral\n");
	  abort();
	}
	popSET = 0;
	locSET = 0;
	if(strcmp(argv[arg], "P") == 0){
	  pop = atoi(argv[arg+1]);
	  popSET = 1;
	  arg += 2;
	  if(ppars[pop].stayNeut == 1){
	    fprintf(errfile,"error:  requesting selection on neutral population %d\n", pop);
	    abort();
	  }
	}
	loc = 0;
	if(strcmp(argv[arg], "L") == 0){
	  if(strcmp(argv[arg+1], "O") == 0){
	    ODD = 1;
	    loc = 1;
	  }
	  else if(strcmp(argv[arg+1], "E") == 0)
	    ODD = 0;
	  else{
	    loc = atol(argv[arg+1]);
	    locSET = 1;
	  }
	  arg += 2;
	}
	if(strcmp(argv[arg], "H") == 0){
	  gpars.allCodom = 0;
	  arg++;
	  ind = atoi(argv[arg]);
	  arg++;
	  if(ind < 0 || ind > 3){
	    fprintf(errfile,"sfs_code error:  hDistType must be 0,1,2, or 3");
	    fprintf(errfile, "\n");
	    abort();
	  }
	  ppars[pop].hDistType[loc] = ind;
	  if(!locSET || ODD >= 0){
	    for(i=loc+1; i<gpars.R; i++){
	      if(ODD == -1 || i%2 == ODD)
		ppars[pop].hDistType[i] = ppars[pop].hDistType[loc];
	    }
	  }
	  if(ind == 0);  /* codominance */
	  else if(ind == 1){    /* point mass */
	    double val = atof(argv[arg]);
	    arg++;
	    if(val < 0){
	      fprintf(errfile,"sfs_code error:  dominance coefficient h < 0... read %lf\n",val);
	      abort();
	    }
	    ppars[pop].codomPop[loc] = 0;
	    ppars[pop].h[loc] = val;
	    ppars[pop].hProbDom[loc] = atof(argv[arg]);
	    arg++;
	    if(!locSET){
	      for(i = 1; i<gpars.R; i++){
		ppars[pop].codomPop[i] = 0;
		ppars[pop].h[i] = ppars[pop].h[0];
		ppars[pop].hProbDom[i] = ppars[pop].hProbDom[0];
	      }
	    }
	  }
	  else if(ind == 2){ /* Gamma distribution */
	    ppars[pop].codomPop[loc] = 0;
	    ppars[pop].hProbGamma[loc] = atof(argv[arg]);
	    arg++;
	    ppars[pop].halpha[loc] = atof(argv[arg]);
	    arg++;
	    ppars[pop].hlambda[loc] = atof(argv[arg]);
	    arg++;
	    if(!locSET){
	      for(i = 1; i<gpars.R; i++){
		ppars[pop].codomPop[i] = 0;
		ppars[pop].hProbGamma[i] = ppars[pop].hProbGamma[0];
		ppars[pop].halpha[i] = ppars[pop].halpha[0];
		ppars[pop].hlambda[i] = ppars[pop].hlambda[0];
	      }
	    }
	  }
	  else if(ind == 3){
	    ppars[pop].codomPop[loc] = 0;
	    ppars[pop].hnormMean[loc] = atof(argv[arg]);
	    arg++;
	    ppars[pop].hnormVar[loc] = atof(argv[arg]);
	    arg++;
	    if(!locSET){
	      for(i = 1; i<gpars.R; i++){
		ppars[pop].codomPop[i] = 0;
		ppars[pop].hnormMean[i] = ppars[pop].hnormMean[0];
		ppars[pop].hnormVar[i] = ppars[pop].hnormVar[0];
	      }
	    }
	  }
	}
	if(!popSET){
	  for(i=0; i<gpars.NPOP; i++){
	    if(i == pop || ppars[i].stayNeut == 1)
	      continue;
	    ppars[i].hDistType[loc] = ppars[pop].hDistType[loc];
	    if(!locSET){
	      for(j=1; j<gpars.R; j++)
		ppars[i].hDistType[j] = ppars[pop].hDistType[0];
	    }
	    if(ind == 1){
	      if(locSET){
		ppars[i].codomPop[loc] = 0;
		ppars[i].h[loc] = ppars[pop].h[loc];
		ppars[i].hProbDom[loc] = ppars[pop].hProbDom[loc];
	      }
	      else{
		for(j=0; j<gpars.R; j++){
		  ppars[i].codomPop[j] = 0;
		  ppars[i].h[j] = ppars[pop].h[0];
		  ppars[i].hProbDom[j] = ppars[pop].hProbDom[0];
		}
	      }
	    }
	    else if(ind == 2){
	      if(locSET){
		ppars[i].codomPop[loc] = 0;
		ppars[i].hProbGamma[loc] = ppars[pop].hProbGamma[loc];
		ppars[i].halpha[loc] = ppars[pop].halpha[loc];
		ppars[i].hlambda[loc] = ppars[pop].hlambda[loc];
	      }
	      else{
		for(j=0; j<gpars.R; j++){
		  ppars[i].codomPop[j] = 0;
		  ppars[i].hProbGamma[j] = ppars[pop].hProbGamma[loc];
		  ppars[i].halpha[j] = ppars[pop].halpha[loc];
		  ppars[i].hlambda[j] = ppars[pop].hlambda[loc];
		}
	      }
	    }
	    else if(ind == 3){
	      if(locSET){
		ppars[i].codomPop[loc] = 0;
		ppars[i].hnormMean[loc] = ppars[pop].hnormMean[loc];
		ppars[i].hnormVar[loc] = ppars[pop].hnormVar[loc];
	      }
	      else{
		for(j=0; j<gpars.R; j++){
		  ppars[i].codomPop[j] = 0;
		  ppars[i].hnormMean[j] = ppars[pop].hnormMean[0];
		  ppars[i].hnormVar[j] = ppars[pop].hnormVar[0];
		}
	      }
	    }
	  }
	}
	if(arg >= argc){
	  fprintf(errfile,"sfs_code error:  insufficient arguments to -W\n");
	  abort();
	}
	ind = atoi(argv[arg]);
	if(ind < 0 || ind > 4){
	  fprintf(errfile,"sfs_code error:  selDistType must be 0,1,2,3, or 4");
	  fprintf(errfile,"\n");
	  abort();
	}
	ppars[pop].selDistType[loc] = ind;
	if(!locSET || ODD >= 0){
	  for(i=loc+1; i<gpars.R; i++){
	    if(ODD == -1 || i%2 == ODD)
	      ppars[pop].selDistType[i] = ppars[pop].selDistType[loc];
	  }
	}
	if(ind == 0)  arg++;  /* neutrality */
	else if(ind == 1){    /* 3-point mass */
	  double val = atof(argv[arg+1]);
	  ppars[pop].neutpop = 0;
	  arg++;
	  if(val < 0)  val *= -1.0;
	  ppars[pop].GAMMA[loc] = val;
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  ppars[pop].ProbPos[loc] = atof(argv[arg]);
	  argcheck(arg+1, arg-4, argc, argv);
	  arg++;
	  ppars[pop].ProbDel[loc] = atof(argv[arg]);
	  ppars[pop].ProbNeut[loc] =
	    1-ppars[pop].ProbPos[loc]-ppars[pop].ProbDel[loc];
	  if(ppars[pop].ProbPos[loc] < 0 || ppars[pop].ProbPos[loc] > 1 ||
	     ppars[pop].ProbDel[loc] < 0 || ppars[pop].ProbDel[loc] > 1 ||
	     ppars[pop].ProbPos[loc] + ppars[pop].ProbDel[loc] > 1){
	    fprintf(errfile,"sfs_code error:  For option --selDistType <pop> \
<GAMMA> <p_pos> <p_neg>, p_pos and p_neg in (0,1) with p_pos + p_neg <= 1.\n");
	    abort();
	  }
	  if(!locSET){
	    for(i = 1; i<gpars.R; i++){
	      ppars[pop].GAMMA[i] = ppars[pop].GAMMA[0];
	      ppars[pop].ProbPos[i] = ppars[pop].ProbPos[0];
	      ppars[pop].ProbDel[i] = ppars[pop].ProbDel[0];
	      ppars[pop].ProbNeut[i] = ppars[pop].ProbNeut[0];
	    }
	  }
	  arg++;
	}
	else if(ind == 2){ /* mixture of Gamma distributions */
	  ppars[pop].neutpop = 0;
	  arg++;
	  ppars[pop].ProbPosGamma[loc] = atof(argv[arg]);
	  if(ppars[pop].ProbPosGamma[loc] < 0){
	    ppars[pop].DoubleNegGamma = 1;
	    ppars[pop].ProbPosGamma[loc] = -ppars[pop].ProbPosGamma[loc];
	  }
	  ppars[pop].ProbNegGamma[loc] = 1.0-ppars[pop].ProbPosGamma[loc];
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  ppars[pop].alphaP[loc] = atof(argv[arg]);
	  argcheck(arg+1, arg-4, argc, argv);
	  arg++;
	  ppars[pop].lambdaP[loc] = atof(argv[arg]);
	  argcheck(arg+1, arg-5, argc, argv);
	  arg++;
	  ppars[pop].alphaN[loc] = atof(argv[arg]);
	  argcheck(arg+1, arg-6, argc, argv);
	  arg++;
	  ppars[pop].lambdaN[loc] = atof(argv[arg]);
	  arg++;
	  if(arg < argc && (isdigit(argv[arg][0]) || argv[arg][0] == '.')){
	    ppars[pop].ProbNegGamma[loc] = atof(argv[arg]);
	    if(ppars[pop].ProbNegGamma[loc]+ppars[pop].ProbPosGamma[loc] < 0 ||
	       ppars[pop].ProbNegGamma[loc]+ppars[pop].ProbPosGamma[loc] > 1){
	      fprintf(errfile, "Error, p_pos + p_neg must be in [0,1]; read %f \
%f at locus %ld\n", ppars[pop].ProbPosGamma[loc], ppars[pop].ProbNegGamma[loc], 
		      loc);
	      abort();
	    }
	    arg++;
	  }
	  if(!locSET){
	    for(i = 1; i<gpars.R; i++){
	      ppars[pop].ProbNegGamma[i] = ppars[pop].ProbNegGamma[0];
	      ppars[pop].ProbPosGamma[i] = ppars[pop].ProbPosGamma[0];
	      ppars[pop].alphaP[i] = ppars[pop].alphaP[0];
	      ppars[pop].lambdaP[i] = ppars[pop].lambdaP[0];
	      ppars[pop].alphaN[i] = ppars[pop].alphaN[0];
	      ppars[pop].lambdaN[i] = ppars[pop].lambdaN[0];
	    }
	  }
	}
	else if(ind == 3){
	  ppars[pop].neutpop = 0;
	  arg++;
	  ppars[pop].normMean[loc] = atof(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  ppars[pop].normVar[loc] = atof(argv[arg]);
	  arg++;
	  if(!locSET){
	    for(i = 1; i<gpars.R; i++){
	      ppars[pop].normMean[i] = ppars[pop].normMean[0];
	      ppars[pop].normVar[i] = ppars[pop].normVar[0];
	    }
	  }
	}
	else if(ind == 4){
	  ppars[pop].neutpop = 0;
	  arg++;
	}
	else{
	  fprintf(errfile,"sfs_code error:  --selDistType <pop> <type> must\
have <type>==0,1,2,3,4\n");
	  abort();
	}
	if(!popSET){
	  for(i=0; i<gpars.NPOP; i++){
	    if(i == pop || ppars[i].stayNeut == 1)
	      continue;
	    ppars[i].selDistType[loc] = ppars[pop].selDistType[loc];
	    if(!locSET){
	      for(j=1; j<gpars.R; j++)
		ppars[i].selDistType[j] = ppars[pop].selDistType[0];
	    }
	    if(ind == 1){
	      ppars[i].neutpop = 0;
	      if(locSET){
		ppars[i].GAMMA[loc] = ppars[pop].GAMMA[loc];
		ppars[i].ProbPos[loc] = ppars[pop].ProbPos[loc];
		ppars[i].ProbDel[loc] = ppars[pop].ProbDel[loc];
		ppars[i].ProbNeut[loc] = ppars[pop].ProbNeut[loc];
	      }
	      else{
		for(j=0; j<gpars.R; j++){
		  ppars[i].GAMMA[j] = ppars[pop].GAMMA[0];
		  ppars[i].ProbPos[j] = ppars[pop].ProbPos[0];
		  ppars[i].ProbDel[j] = ppars[pop].ProbDel[0];
		  ppars[i].ProbNeut[j] = ppars[pop].ProbNeut[0];
		}
	      }
	    }
	    else if(ind == 2){
	      ppars[i].neutpop = 0;
	      if(locSET){
		ppars[i].ProbNegGamma[loc] = ppars[pop].ProbNegGamma[loc];
		ppars[i].ProbPosGamma[loc] = ppars[pop].ProbPosGamma[loc];
		ppars[i].alphaN[loc] = ppars[pop].alphaN[loc];
		ppars[i].lambdaN[loc] = ppars[pop].lambdaN[loc];
		ppars[i].alphaP[loc] = ppars[pop].alphaP[loc];
		ppars[i].lambdaP[loc] = ppars[pop].lambdaP[loc];
	      }
	      else{
		for(j=0; j<gpars.R; j++){
		  ppars[i].ProbNegGamma[j] = ppars[pop].ProbNegGamma[loc];
		  ppars[i].ProbPosGamma[j] = ppars[pop].ProbPosGamma[loc];
		  ppars[i].alphaN[j] = ppars[pop].alphaN[loc];
		  ppars[i].lambdaN[j] = ppars[pop].lambdaN[loc];
		  ppars[i].alphaP[j] = ppars[pop].alphaP[loc];
		  ppars[i].lambdaP[j] = ppars[pop].lambdaP[loc];
		}
	      }
	    }
	    else if(ind == 3){
	      ppars[i].neutpop = 0;
	      if(locSET){
		ppars[i].normMean[loc] = ppars[pop].normMean[loc];
		ppars[i].normVar[loc] = ppars[pop].normVar[loc];
	      }
	      else{
		for(j=0; j<gpars.R; j++){
		  ppars[i].normMean[j] = ppars[pop].normMean[0];
		  ppars[i].normVar[j] = ppars[pop].normVar[0];
		}
	      }
	    }
	    else if(ind == 4){
	      ppars[i].neutpop = 0;
	    }
	  }
	}
      }
      else if(strcmp(args,"TW") == 0){ /* make change selDistType */
	ind = arg; /* save it for use in argcheck() */
	argcheck(arg+1, ind, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 7;
	tmpEvolEvent.tau = atof(argv[arg]);
	tmpEvolEvent.nu = -1; /* dummy variable for tracking option */
	argcheck(arg+1, ind, argc, argv);

	arg++;
	pop = 0;
	while(ppars[pop].stayNeut == 1 && pop<gpars.NPOP){
	  pop++;
	}
	if(pop >= gpars.NPOP){
	  fprintf(errfile,"error:  requesting selection -TW when all populations are neutral\n");
	  abort();
	}
	
	popSET = 0;
	locSET = 0;
	if(strcmp(argv[arg], "P") == 0){
	  pop = tmpEvolEvent.popi = atoi(argv[arg+1]);
	  if(ppars[pop].stayNeut == 1){
	    fprintf(errfile,"error: requesting selection -TW on neutral population %d\n", pop);
	    abort();
	  }
	  popSET = 1;
	  arg += 2;
	  if(ppars[pop].stayNeut == 1){
	    fprintf(errfile,"error:  requesting selection on neutral population %d\n", pop);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.popi = -1;
	}
	loc = 0;
	if(strcmp(argv[arg], "L") == 0){
	  loc = tmpEvolEvent.locus = atol(argv[arg+1]);
	  locSET = 1;
	  arg += 2;
	}
	else
	  tmpEvolEvent.locus = -1;
	
	if(tmpEvolEvent.newP.codomPop == NULL){
	  assert(tmpEvolEvent.newP.codomPop =
		 malloc(gpars.R*sizeof(*tmpEvolEvent.newP.codomPop)));
	  assert(tmpEvolEvent.newP.hDistType =
		 malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hDistType)));
	}
	tmpEvolEvent.newP.codomPop[loc] = 1;
	tmpEvolEvent.newP.hDistType[loc] = 0;

	if(strcmp(argv[arg], "H") == 0){
	  arg++;
	  if(tmpEvolEvent.newP.hDistType == NULL){
	    assert(tmpEvolEvent.newP.hDistType =
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hDistType)));
	  }
	  ind = tmpEvolEvent.newP.hDistType[loc] = atoi(argv[arg]);
	  arg++;
	  if(!locSET){
	    for(i=1; i<gpars.R; i++){
	      tmpEvolEvent.newP.hDistType[i] = tmpEvolEvent.newP.hDistType[0];
	      tmpEvolEvent.newP.codomPop[i] = tmpEvolEvent.newP.codomPop[0];
	    }
	  }
	  if(ind < 0 || ind >= 4){
	    fprintf(errfile,"sfs_code error:  hDistType must be 0,1,2, or 3");
	    fprintf(errfile,"\n");
	    abort();
	  }
	
	  if(ind == 0);  /* codominance */
	  else if(ind == 1){    /* point mass */
	    double val = atof(argv[arg]);
	    arg++;
	    if(tmpEvolEvent.newP.h == NULL){
	      assert(tmpEvolEvent.newP.h = 
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.h)));
	      assert(tmpEvolEvent.newP.hProbDom = 
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hProbDom)));
	    }
	    if(val < 0){
	      fprintf(errfile,"sfs_code error:  h value must be non-negative.  Read %lf\n",val);
	      abort();
	    }
	    tmpEvolEvent.newP.codomPop[loc] = 0;
	    tmpEvolEvent.newP.h[loc] = val;
	    tmpEvolEvent.newP.hProbDom[loc] = atof(argv[arg]);
	    arg++;
	    if(!locSET){
	      for(i = 1; i<gpars.R; i++){
		tmpEvolEvent.newP.codomPop[i] = 0;
		tmpEvolEvent.newP.h[i] = tmpEvolEvent.newP.h[0];
		tmpEvolEvent.newP.hProbDom[i] = tmpEvolEvent.newP.hProbDom[0];
	      }
	    }
	  }
	  else if(ind == 2){ /* Gamma distributions */
	    if(tmpEvolEvent.newP.hProbGamma == NULL){
	      assert(tmpEvolEvent.newP.hProbGamma = 
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hProbGamma)));
	      assert(tmpEvolEvent.newP.halpha = 
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.halpha)));
	      assert(tmpEvolEvent.newP.hlambda = 
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hlambda)));
	    }
	    
	    tmpEvolEvent.newP.codomPop[loc] = 0;
	    tmpEvolEvent.newP.hProbGamma[loc] = atof(argv[arg]);
	    arg++;
	    tmpEvolEvent.newP.halpha[loc] = atof(argv[arg]);
	    arg++;
	    tmpEvolEvent.newP.hlambda[loc] = atof(argv[arg]);
	    arg++;
	    if(!locSET){
	      for(i=1; i<gpars.R; i++){
		tmpEvolEvent.newP.codomPop[i] = 0;
		tmpEvolEvent.newP.hProbGamma[i] = 
		  tmpEvolEvent.newP.hProbGamma[0];
		tmpEvolEvent.newP.halpha[i] = tmpEvolEvent.newP.halpha[0];
		tmpEvolEvent.newP.hlambda[i] = tmpEvolEvent.newP.hlambda[0];
	      }
	    }
	  }
	  else if(ind == 3){
	    if(tmpEvolEvent.newP.hnormMean == NULL){
	      assert(tmpEvolEvent.newP.hnormMean =
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hnormMean)));
	      assert(tmpEvolEvent.newP.hnormVar =
		     malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hnormVar)));
	    }
	    tmpEvolEvent.newP.codomPop[loc] = 0;
	    tmpEvolEvent.newP.hnormMean[loc] = atof(argv[arg]);
	    arg++;
	    tmpEvolEvent.newP.hnormVar[loc] = atof(argv[arg]);
	    arg++;
	    if(!locSET){
	      for(i=1; i<gpars.R; i++){
		tmpEvolEvent.newP.codomPop[i] = 0;
		tmpEvolEvent.newP.hnormMean[i] = tmpEvolEvent.newP.hnormMean[0];
		tmpEvolEvent.newP.hnormVar[i] = tmpEvolEvent.newP.hnormVar[0];
	      }
	    }
	  }
	  else{
	    fprintf(errfile,"sfs_code error:  --hDistType <pop> <type> must\
have <type>==0,1,2,3\n");
	    abort();
	  }
	}
	
	if(strcmp(argv[arg], "F") == 0){
	  tmpEvolEvent.parIndex = 20;
	  tmpEvolEvent.freq = atof(argv[arg+1]);
	  if(tmpEvolEvent.popi == -1){
	    tmpEvolEvent.popi = 0;
	  }
	  arg += 2;
	  if(arg<argc && strcmp(argv[arg], "w") == 0){
	    gpars.stayNeut = 1;
	    arg++;
	  }
	  if(arg<argc && strcmp(argv[arg], "T") == 0){ /* track locus */
	    tmpEvolEvent.nu = 1; /* dummy to turn on tracking  */
	    if(toBeTracked != 0){
	      fprintf(errfile,"Can only use --trackTrajectory once.  This \
includes the track allele option in -TW\n");
	      abort();
	    }
	    toBeTracked = 1;
	    arg++;
	    if(arg<argc && strcmp(argv[arg],"R") == 0){
	      argcheck(arg+1, ind, argc, argv);
	      arg++;
	      gpars.trajMinFreq = atof(argv[arg]);
	      argcheck(arg+1, ind, argc, argv);
	      arg++;
	      gpars.trajMaxFreq = atof(argv[arg]);
	      arg++;
	      if(arg<argc && strcmp(argv[arg],"S") == 0){
		gpars.stopHitFreq = 1;
		arg++;
	      }
	    }
	    if(arg<argc && strcmp(argv[arg],"A") == 0){
	      gpars.autoRestart = 1;
	      arg++;
	      if(arg<argc && strcmp(argv[arg], "G") == 0){
		gpars.gensAfterFix = atol(argv[arg+1]);
		arg+=2;
	      }
	    }
	    if(arg<argc && strcmp(argv[arg],"M") == 0){
	      argcheck(arg+1, ind, argc, argv);
	      arg++;
	      if(INITIALSEED == 0)
		gpars.trajMaxReps = atol(argv[arg]);
	      arg++;
	    }
	    if(arg<argc && strcmp(argv[arg],"F") == 0){
	      argcheck(arg+1,ind,argc,argv);
	      arg++;
	      if(strcmp(argv[arg],"a") == 0){
		arg++;
		if(!(trajFile=fopen(argv[arg],"a"))){
		  fprintf(errfile,"cannot append to %s\n",argv[arg]);
		  abort();
		}
	      }
	      else{
		if(trajFile == NULL){
		  if(!(trajFile=fopen(argv[arg],"w"))){
		    fprintf(errfile,"cannot write to %s\n",argv[arg]);
		    abort();
		  }
		}
	      }
	      arg++;
	    }
	  }
	}

	if(tmpEvolEvent.newP.selDistType == NULL){
	  assert(tmpEvolEvent.newP.selDistType =
		 malloc(gpars.R*sizeof(*tmpEvolEvent.newP.selDistType)));
	}
	ind = tmpEvolEvent.newP.selDistType[loc] = atoi(argv[arg]);
	arg++;
	
	if(!locSET){
	  for(i=1; i<gpars.R; i++)
	    tmpEvolEvent.newP.selDistType[i] = tmpEvolEvent.newP.selDistType[0];
	}
	if(ind < 0 || ind > 4){
	  fprintf(errfile,"sfs_code error:  selDistType must be 0,1,2,3, or 4");
	  fprintf(errfile,"\n");
	  abort();
	}
	
	if(ind == 0);  /* neutrality */
	else if(ind == 1){    /* 3-point mass */
	  double val = atof(argv[arg]);
	  //arg++;
	  if(tmpEvolEvent.newP.GAMMA == NULL){
	    assert(tmpEvolEvent.newP.GAMMA = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.GAMMA)));
	    assert(tmpEvolEvent.newP.ProbPos = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.ProbPos)));
	    assert(tmpEvolEvent.newP.ProbDel = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.ProbDel)));
	    assert(tmpEvolEvent.newP.ProbNeut = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.ProbNeut)));
	  }
	  if(val < 0)  val *= -1.0;
	  tmpEvolEvent.newP.GAMMA[loc] = val;
	  tmpEvolEvent.newP.ProbPos[loc] = atof(argv[arg]);
	  arg++;
	  tmpEvolEvent.newP.ProbDel[loc] = atof(argv[arg]);
	  arg++;
	  
	  tmpEvolEvent.newP.ProbNeut[loc] =
	    1-tmpEvolEvent.newP.ProbPos[loc]-tmpEvolEvent.newP.ProbDel[loc];
	  if(tmpEvolEvent.newP.ProbPos[loc] > 1 ||
	     tmpEvolEvent.newP.ProbDel[loc] > 1 ||
	     (tmpEvolEvent.newP.ProbPos[loc] + 
	      tmpEvolEvent.newP.ProbDel[loc]) > 1){
	    fprintf(errfile,"sfs_code error:  For option --selDistType <pop> \
<GAMMA> <p_pos> <p_neg>, p_pos and p_neg in (0,1) with p_pos + p_neg <= 1.\n");
	    abort();
	  }
	  if(!locSET){
	    for(i = 1; i<gpars.R; i++){
	      tmpEvolEvent.newP.GAMMA[i] = tmpEvolEvent.newP.GAMMA[0];
	      tmpEvolEvent.newP.ProbPos[i] = tmpEvolEvent.newP.ProbPos[0];
	      tmpEvolEvent.newP.ProbDel[i] = tmpEvolEvent.newP.ProbDel[0];
	      tmpEvolEvent.newP.ProbNeut[i] = tmpEvolEvent.newP.ProbNeut[0];
	    }
	  }
	}
	else if(ind == 2){ /* mixture of Gamma distributions */
	  //arg++;
	  if(tmpEvolEvent.newP.ProbNegGamma == NULL){
	    assert(tmpEvolEvent.newP.ProbNegGamma = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.ProbNegGamma)));
	    assert(tmpEvolEvent.newP.ProbPosGamma = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.ProbPosGamma)));
	    assert(tmpEvolEvent.newP.alphaP = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.alphaP)));
	    assert(tmpEvolEvent.newP.alphaN = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.alphaN)));
	    assert(tmpEvolEvent.newP.lambdaP = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.lambdaP)));
	    assert(tmpEvolEvent.newP.lambdaN = 
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.lambdaN)));
	  }

	  tmpEvolEvent.newP.ProbPosGamma[loc] = atof(argv[arg]);
	  if(tmpEvolEvent.newP.ProbPosGamma[loc] < 0){
	    ppars[pop].DoubleNegGamma = 1;
	    tmpEvolEvent.newP.ProbPosGamma[loc] =
	      -tmpEvolEvent.newP.ProbPosGamma[loc];
	  }
	  tmpEvolEvent.newP.ProbNegGamma[loc] = 
	    1.0-tmpEvolEvent.newP.ProbPosGamma[loc];
	  arg++;
	  tmpEvolEvent.newP.alphaP[loc] = atof(argv[arg]);
	  arg++;
	  tmpEvolEvent.newP.lambdaP[loc] = atof(argv[arg]);
	  arg++;
	  tmpEvolEvent.newP.alphaN[loc] = atof(argv[arg]);
	  arg++;
	  tmpEvolEvent.newP.lambdaN[loc] = atof(argv[arg]);
	  arg++;

	  if(arg<argc && (isdigit(argv[arg][0]) || argv[arg][0] == '.')){
	    tmpEvolEvent.newP.ProbNegGamma[loc] = atof(argv[arg]);
	    if(tmpEvolEvent.newP.ProbNegGamma[loc] +
	       tmpEvolEvent.newP.ProbPosGamma[loc] < 0 ||
	       tmpEvolEvent.newP.ProbNegGamma[loc] +
	       tmpEvolEvent.newP.ProbPosGamma[loc] > 1){
	      fprintf(errfile, "Error, p_pos + p_neg must be in [0,1]; read %f \
%f at locus %ld\n", tmpEvolEvent.newP.ProbPosGamma[loc], 
		      tmpEvolEvent.newP.ProbNegGamma[loc], loc);
	      abort();
	    }
	    arg++;
	  }
	  if(!locSET){
	    for(i=1; i<gpars.R; i++){
	      tmpEvolEvent.newP.ProbNegGamma[i] = 
		tmpEvolEvent.newP.ProbNegGamma[0];
	      tmpEvolEvent.newP.ProbPosGamma[i] = 
		tmpEvolEvent.newP.ProbPosGamma[0];
	      tmpEvolEvent.newP.alphaP[i] = tmpEvolEvent.newP.alphaP[0];
	      tmpEvolEvent.newP.alphaN[i] = tmpEvolEvent.newP.alphaN[0];
	      tmpEvolEvent.newP.lambdaP[i] = tmpEvolEvent.newP.lambdaP[0];
	      tmpEvolEvent.newP.lambdaN[i] = tmpEvolEvent.newP.lambdaN[0];
	    }
	  }
	}
	else if(ind == 3){
	  if(tmpEvolEvent.newP.normMean == NULL){
	    assert(tmpEvolEvent.newP.normMean =
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.normMean)));
	    assert(tmpEvolEvent.newP.normVar =
		   malloc(gpars.R*sizeof(*tmpEvolEvent.newP.normVar)));
	  }
	  //arg++;
	  tmpEvolEvent.newP.normMean[loc] = atof(argv[arg]);
	  arg++;
	  tmpEvolEvent.newP.normVar[loc] = atof(argv[arg]);
	  arg++;
	  if(!locSET){
	    for(i=1; i<gpars.R; i++){
	      tmpEvolEvent.newP.normMean[i] = tmpEvolEvent.newP.normMean[0];
	      tmpEvolEvent.newP.normVar[i] = tmpEvolEvent.newP.normVar[0];
	    }
	  }
	}
	else if(ind == 4)  arg++;
	else{
	  fprintf(errfile,"sfs_code error:  --selDistType <pop> <type> must\
have <type>==0,1,2,3,4\n");
	  abort();
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);

	if(tmpEvolEvent.newP.GAMMA != NULL){
	  free(tmpEvolEvent.newP.GAMMA);
	  tmpEvolEvent.newP.GAMMA = NULL;
	  free(tmpEvolEvent.newP.ProbPos);
	  free(tmpEvolEvent.newP.ProbDel);
	  free(tmpEvolEvent.newP.ProbNeut);
	}
	if(tmpEvolEvent.newP.ProbNegGamma != NULL){
	  free(tmpEvolEvent.newP.ProbNegGamma);
	  tmpEvolEvent.newP.ProbNegGamma = NULL;
	  free(tmpEvolEvent.newP.ProbPosGamma);
	  tmpEvolEvent.newP.ProbPosGamma = NULL;
	  free(tmpEvolEvent.newP.alphaP);
	  free(tmpEvolEvent.newP.alphaN);
	  free(tmpEvolEvent.newP.lambdaP);
	  free(tmpEvolEvent.newP.lambdaN);
	}
	if(tmpEvolEvent.newP.normMean != NULL){
	  free(tmpEvolEvent.newP.normMean);
	  tmpEvolEvent.newP.normMean = NULL;
	  free(tmpEvolEvent.newP.normVar);
	}
	if(tmpEvolEvent.newP.selDistType != NULL){
	  free(tmpEvolEvent.newP.selDistType);
	  tmpEvolEvent.newP.selDistType = NULL;
	}
	if(tmpEvolEvent.newP.h != NULL){
	  free(tmpEvolEvent.newP.h);
	  tmpEvolEvent.newP.h = NULL;
	  free(tmpEvolEvent.newP.hProbDom);
	}
	if(tmpEvolEvent.newP.hProbGamma != NULL){
	  free(tmpEvolEvent.newP.hProbGamma);
	  tmpEvolEvent.newP.hProbGamma = NULL;
	  free(tmpEvolEvent.newP.halpha);
	  free(tmpEvolEvent.newP.hlambda);
	}
	if(tmpEvolEvent.newP.hnormMean != NULL){
	  free(tmpEvolEvent.newP.hnormMean);
	  tmpEvolEvent.newP.hnormMean = NULL;
	  free(tmpEvolEvent.newP.hnormVar);
	}
	if(tmpEvolEvent.newP.hDistType != NULL){
	  free(tmpEvolEvent.newP.hDistType);
	  tmpEvolEvent.newP.hDistType = NULL;
	}
	if(tmpEvolEvent.newP.codomPop != NULL){
	  free(tmpEvolEvent.newP.codomPop);
	  tmpEvolEvent.newP.codomPop = NULL;
	}
      }
      else if(strcmp(args,"-neutPop") == 0 || strcmp(args,"w") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	pop = atoi(argv[arg++]);

	ppars[pop].neutpop = 1;
	ppars[pop].stayNeut = 1;
	for(i=0; i<gpars.R; i++){
	  ppars[pop].codomPop[i] = 1;
	  ppars[pop].selDistType[i] = 0;
	  ppars[pop].f0[i] = 1.0;
	  ppars[pop].hDistType[i] = 0;
	}
      }
      else if(strcmp(args,"Tw") == 0){ /* make pop. neutral */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 7;
	tmpEvolEvent.tau = atof(argv[arg++]);
	argcheck(arg, arg-2, argc, argv);
 	tmpEvolEvent.popi = atoi(argv[arg++]);
	if(tmpEvolEvent.newP.selDistType == NULL){
	  assert(tmpEvolEvent.newP.selDistType = 
		 malloc(gpars.R*sizeof(*tmpEvolEvent.newP.selDistType)));
	}
	if(tmpEvolEvent.newP.hDistType == NULL){
	  assert(tmpEvolEvent.newP.hDistType = 
		 malloc(gpars.R*sizeof(*tmpEvolEvent.newP.hDistType)));
	}
	if(tmpEvolEvent.newP.codomPop == NULL){
	  assert(tmpEvolEvent.newP.codomPop = 
		 malloc(gpars.R*sizeof(*tmpEvolEvent.newP.codomPop)));
	}
	for(i=0; i<gpars.R; i++){
	  tmpEvolEvent.newP.selDistType[i] = 0;
	  tmpEvolEvent.newP.hDistType[i] = 0;
	  tmpEvolEvent.newP.codomPop[i] = 1;
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
	free(tmpEvolEvent.newP.selDistType);
	tmpEvolEvent.newP.selDistType = NULL;
	free(tmpEvolEvent.newP.hDistType);
	tmpEvolEvent.newP.hDistType = NULL;
	free(tmpEvolEvent.newP.codomPop);
	tmpEvolEvent.newP.codomPop = NULL;
      }
      else if(strcmp(args,"-constraint") == 0 || strcmp(args,"c") == 0){
	argcheck(arg+1, arg, argc, argv); /* constrained evolution */
	arg++;
	pop = 0;
	while(ppars[pop].stayNeut == 1 && pop<gpars.NPOP){
	  pop++;
	}
	if(pop >= gpars.NPOP){
	  fprintf(errfile,"error:  requesting constraint when all populations are neutral\n");
	  abort();
	}
	
	popSET = 0;
	locSET = 0;
	if(strcmp(argv[arg], "P") == 0){
	  pop = atoi(argv[arg+1]);
	  if(ppars[pop].stayNeut == 1){
	    fprintf(errfile,"error: requesting constraint on neutral population %d\n", pop);
	    abort();
	  }
	  popSET = 1;
	  arg += 2;
	  if(pop >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",pop,gpars.NPOP-1);
	    abort();
	  }
	}
	loc = 0;
	if(strcmp(argv[arg], "L") == 0){
	  loc = atoi(argv[arg+1]);
	  locSET = 1;
	  arg += 2;
	  if(loc >= gpars.R){
	    fprintf(errfile,"sfs_code error: cannot set parameters for locus %ld \
(max = %ld)\n",loc,gpars.R-1);
	    abort();
	  }
 	}
	ppars[pop].f0[loc] = atof(argv[arg]);
	if(ppars[pop].f0[loc] < 1-FLT_MIN)
	  ppars[pop].neutpop = 0;
	arg++;
	if(!locSET){
	  for(i=1; i<gpars.R; i++)
	    ppars[pop].f0[i] = ppars[pop].f0[0];
	}
	if(!popSET){
	  for(i=0; i<gpars.NPOP; i++){
	    if(i == pop || ppars[i].stayNeut == 1)
	      continue;
	    if(locSET){
	      ppars[i].f0[loc] = ppars[pop].f0[loc];
	      if(ppars[i].f0[loc] < 1-FLT_MIN)
		ppars[i].neutpop = 0;
	    }
	    else{
	      for(j=0; j<gpars.R; j++){
		ppars[i].f0[j] = ppars[pop].f0[0];
		if(ppars[i].f0[j] < 1-FLT_MIN)
		  ppars[i].neutpop = 0;
	      }
	    }
	  }
	}
      }
      else if(strcmp(args,"Tc") == 0){
	if(tmpEvolEvent.newP.f0 == NULL){
	  assert(tmpEvolEvent.newP.f0 =
		 malloc(gpars.R*sizeof(tmpEvolEvent.newP.f0)));
	}
	argcheck(arg+1, arg, argc, argv); /* constrained evolution */
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 10;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	pop = 0;
	popSET = 0;
	locSET = 0;
	if(strcmp(argv[arg], "P") == 0){
	  pop = tmpEvolEvent.popi = atoi(argv[arg+1]);
	  popSET = 1;
	  arg += 2;
	  if(pop >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",pop,gpars.NPOP-1);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.popi = -1;
	}
	loc = 0;
	if(strcmp(argv[arg], "L") == 0){
	  loc = tmpEvolEvent.locus = atol(argv[arg+1]);
	  locSET = 1;
	  arg += 2;
	  if(loc >= gpars.R){
	    fprintf(errfile,"sfs_code error: cannot set parameters for locus %ld \
(max = %ld)\n",loc,gpars.R-1);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.locus = -1;
 	}
	
	tmpEvolEvent.newP.f0[loc] = atof(argv[arg++]);
	if(!locSET){
	  for(i=1; i<gpars.R; i++){
	    tmpEvolEvent.newP.f0[i] = tmpEvolEvent.newP.f0[0];
	  }
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
	free(tmpEvolEvent.newP.f0);
	tmpEvolEvent.newP.f0 = NULL;
      }
      else if(strcmp(args,"-fittest") == 0){
	arg++;
	gpars.FITTEST = 1;
      }
      else if(strcmp(args,"-BURN") == 0 || strcmp(args,"B") == 0){
	argcheck(arg+1, arg, argc, argv); 
	arg++;
	gpars.BURN = atof(argv[arg++]);
      }
      else if(strcmp(args,"-BURN2") == 0 || strcmp(args,"b") == 0){
	argcheck(arg+1, arg, argc, argv); 
	arg++;
	gpars.BURN2 = atof(argv[arg++]);
      }
      else if(strcmp(args,"-length") == 0 || strcmp(args,"L") == 0){
	if(gpars.SKIPSITES != NULL){
	  fprintf(errfile,"error, use --length option before --mutation\n");
	  abort();
	}
	if(gpars.trajLoc != -1){
	  fprintf(errfile,"error, use --length option before --trackTrajectory\n");
	  abort();
	}
	argcheck(arg+1, arg, argc, argv);
	arg++;
	gpars.R = atol(argv[arg++]);
	if(gpars.R > 1){
	  if(aBeforeL != -1){
	    fprintf(errfile,"sfs_code error: You must specify --length (-L) before --annotate (-a).\n");
	    abort();
	  }
	  assert(gpars.L = realloc(gpars.L,gpars.R*sizeof(*gpars.L)));
	  assert(gpars.LINK = realloc(gpars.LINK,
				      (gpars.R-1)*sizeof(*gpars.LINK)));
	  assert(gpars.PRINTLOC = realloc(gpars.PRINTLOC,
					  (gpars.R)*sizeof(*gpars.PRINTLOC)));
	  assert(gpars.ANNOTATE =
		 realloc(gpars.ANNOTATE, gpars.R*sizeof(*gpars.ANNOTATE)));
	  assert(gpars.SEX = realloc(gpars.SEX, gpars.R*sizeof(*gpars.SEX)));

	  assert(mutAr = realloc(mutAr, gpars.R*sizeof(struct mutArray*)));
	  
	  for(i=0; i<gpars.NPOP; i++){
	    assert(ppars[i].h = realloc(ppars[i].h, gpars.R*sizeof(*ppars[i].h)));
	    assert(ppars[i].hProbDom = 
		   realloc(ppars[i].hProbDom,gpars.R*sizeof(*ppars[i].hProbDom)));
	    assert(ppars[i].codomPop =
		   realloc(ppars[i].codomPop,
			   gpars.R*sizeof(*ppars[i].codomPop)));
	    assert(ppars[i].hDistType =
		   realloc(ppars[i].hDistType,
			   gpars.R*sizeof(*ppars[i].hDistType)));
	    assert(ppars[i].hProbGamma =
		   realloc(ppars[i].hProbGamma,
			   gpars.R*sizeof(*ppars[i].hProbGamma)));
	    assert(ppars[i].halpha =
		   realloc(ppars[i].halpha,
			   gpars.R*sizeof(*ppars[i].halpha)));
	    assert(ppars[i].hlambda =
		   realloc(ppars[i].hlambda,
			   gpars.R*sizeof(*ppars[i].hlambda)));
	    assert(ppars[i].hnormMean =
		   realloc(ppars[i].hnormMean,
			   gpars.R*sizeof(*ppars[i].hnormMean)));
	    assert(ppars[i].hnormVar =
		   realloc(ppars[i].hnormVar,
			   gpars.R*sizeof(*ppars[i].hnormVar)));
	    
	    assert(ppars[i].GAMMA = 
		   realloc(ppars[i].GAMMA, gpars.R*sizeof(*ppars[i].GAMMA)));
	    assert(ppars[i].ProbPos = 
		   realloc(ppars[i].ProbPos,gpars.R*sizeof(*ppars[i].ProbPos)));
	    assert(ppars[i].ProbDel = 
		   realloc(ppars[i].ProbDel,gpars.R*sizeof(*ppars[i].ProbDel)));
	    assert(ppars[i].ProbNeut =
		   realloc(ppars[i].ProbNeut,
			   gpars.R*sizeof(*ppars[i].ProbNeut)));
	    assert(ppars[i].selDistType =
		   realloc(ppars[i].selDistType,
			   gpars.R*sizeof(*ppars[i].selDistType)));
	    assert(ppars[i].ProbNegGamma =
		   realloc(ppars[i].ProbNegGamma,
			   gpars.R*sizeof(*ppars[i].ProbNegGamma)));
	    assert(ppars[i].ProbPosGamma =
		   realloc(ppars[i].ProbPosGamma,
			   gpars.R*sizeof(*ppars[i].ProbPosGamma)));
	    assert(ppars[i].alphaN =
		   realloc(ppars[i].alphaN,
			   gpars.R*sizeof(*ppars[i].alphaN)));
	    assert(ppars[i].lambdaN =
		   realloc(ppars[i].lambdaN,
			   gpars.R*sizeof(*ppars[i].lambdaN)));
	    assert(ppars[i].alphaP =
		   realloc(ppars[i].alphaP,
			   gpars.R*sizeof(*ppars[i].alphaP)));
	    assert(ppars[i].lambdaP =
		   realloc(ppars[i].lambdaP,
			   gpars.R*sizeof(*ppars[i].lambdaP)));
	    assert(ppars[i].normMean =
		   realloc(ppars[i].normMean,
			   gpars.R*sizeof(*ppars[i].normMean)));
	    assert(ppars[i].normVar =
		   realloc(ppars[i].normVar,
			   gpars.R*sizeof(*ppars[i].normVar)));
	    assert(ppars[i].f0 =
		   realloc(ppars[i].f0,
			   gpars.R*sizeof(*ppars[i].f0)));
	    
	    for(j=1; j<gpars.R; j++){
	      ppars[i].codomPop[j] = ppars[i].codomPop[0];
	      ppars[i].hDistType[j] = ppars[i].hDistType[0];
	      ppars[i].h[j] = ppars[i].h[0];
	      ppars[i].hProbDom[j] = ppars[i].hProbDom[0];
	      ppars[i].hProbGamma[j] = ppars[i].hProbGamma[0];
	      ppars[i].halpha[j] = ppars[i].halpha[0];
	      ppars[i].hlambda[j] = ppars[i].hlambda[0];
	      ppars[i].hnormMean[j] = ppars[i].hnormMean[0];
	      ppars[i].hnormVar[j] = ppars[i].hnormVar[0];

	      ppars[i].selDistType[j] = ppars[i].selDistType[0];
	      ppars[i].GAMMA[j] = ppars[i].GAMMA[0];
	      ppars[i].ProbNeut[j] = ppars[i].ProbNeut[0];
	      ppars[i].ProbDel[j] = ppars[i].ProbDel[0];
	      ppars[i].ProbPos[j] = ppars[i].ProbPos[0];
	      ppars[i].ProbNegGamma[j] = ppars[i].ProbNegGamma[0];
	      ppars[i].ProbPosGamma[j] = ppars[i].ProbPosGamma[0];
	      ppars[i].alphaN[j] = ppars[i].alphaN[0];
	      ppars[i].lambdaN[j] = ppars[i].lambdaN[0];
	      ppars[i].alphaP[j] = ppars[i].alphaP[0];
	      ppars[i].lambdaP[j] = ppars[i].lambdaP[0];
	      ppars[i].normMean[j] = ppars[i].normMean[0];
	      ppars[i].normVar[j] = ppars[i].normVar[0];
	      ppars[i].f0[j] = ppars[i].f0[0];
	      if(i == 0){
		assert(mutAr[j] = malloc(sizeof(**mutAr)));
		/* create first mutation as default */
		mutAr[j]->numMuts = 1;
		mutAr[j]->mutIndex = 0;
		assert(mutAr[j]->muts = malloc(sizeof(struct history*)));
		mutAr[j]->muts[0] = popTrash();
		assert(mutAr[j]->muts[0]->event=malloc(sizeof(struct event)));
		mutAr[j]->muts[0]->event->site = 0;
		assert(mutAr[j]->muts[0]->event->genFix =
		       malloc(gpars.NPOP*sizeof(long)));
		assert(mutAr[j]->muts[0]->event->genDead =
		       malloc(gpars.NPOP*sizeof(long)));
		mutAr[j]->muts[0]->event->nSites = 0;
		mutAr[j]->muts[0]->event->nucs = NULL;
		assert(mutAr[j]->muts[0]->event->fixed = 
		       malloc(gpars.NPOP*sizeof(char)));
		mutAr[j]->muts[0]->event->free = '1'; /* free by default */
		mutAr[j]->muts[0]->event->index = 0;
		mutAr[j]->muts[0]->event->checkGen = -1;
		mutAr[j]->muts[0]->event->checkStep = -1;
		mutAr[j]->muts[0]->event->checkPop = -1;
		assert(mutAr[j]->muts[0]->event->nextHap =
		       malloc(gpars.NPOP*sizeof(long)));
		assert(mutAr[j]->muts[0]->event->numCarriers = 
		       malloc(gpars.NPOP*sizeof(long)));
		assert(mutAr[j]->muts[0]->event->maxHaps =
		       malloc(gpars.NPOP*sizeof(long)));
		assert(mutAr[j]->muts[0]->event->hapFreq=
		       malloc(gpars.NPOP*sizeof(long**)));
	      }
	      mutAr[j]->muts[0]->event->genFix[i] = 0;
	      mutAr[j]->muts[0]->event->genDead[i] = -LONG_MAX;
	      mutAr[j]->muts[0]->event->fixed[i] = '0'; /* free by default */
	      mutAr[j]->muts[0]->event->nextHap[i] = 0;
	      mutAr[j]->muts[0]->event->numCarriers[i] = 0;
	      mutAr[j]->muts[0]->event->maxHaps[i] = 1;
	      assert(mutAr[j]->muts[0]->event->hapFreq[i] =
		     malloc(sizeof(long*))); 
	      mutAr[j]->muts[0]->event->hapFreq[i][0] = NULL;
	    }
	  }
	}
	argcheck(arg, arg-2, argc, argv);
	gpars.L[0] = atol(argv[arg++]);
	gpars.LINK[0] = 0.0;
	gpars.PRINTLOC[0] = '0';
	if(arg >= argc || argv[arg][0] == '-')  ind = 0;
	else  ind = 1;
	for(i=1; i<gpars.R; i++){
	  if(ind == 0)
	    gpars.L[i] = gpars.L[0];
	  else if(ind < i)
	    gpars.L[i] = gpars.L[i%ind];
	  else
	    if(strcmp(argv[arg], "R") == 0){ /* start repeating */
	      gpars.L[i] = gpars.L[0];
	      arg++;
	    }
	    else{
	      if(strcmp(argv[arg], "-1") != 0)
		argcheck(arg, arg-i-2, argc, argv);
	      gpars.L[i] = atol(argv[arg++]);
	      ind++;
	    }
	  if(i<gpars.R-1)  gpars.LINK[i] = 0.0;
	  gpars.ANNOTATE[i] = 'C';
	  gpars.SEX[i] = '0';
	  gpars.PRINTLOC[i] = '0';
	  if(gpars.L[i]*DBL_MIN > 1){
	    fprintf(errfile,"Due to random number generator, can only simulate \
loci with length up to %e, locus %ld is of length %ld\n",floor(1/FLT_MIN),
		    i,gpars.L[i]);
	    exit(1);
	  }
	}
	if(gpars.L[0] >= 1/DBL_MIN){
	  fprintf(errfile,"Due to random number generator, can only simulate \
loci with length up to %e, locus %d is of length %ld\n",floor(1/FLT_MIN),
		  0,gpars.L[0]);
	  fprintf(errfile,"Split sequence into multiple (shorter) loci.\n");
	  abort();
	}
      }
      else if(strcmp(args,"-linkage") == 0 || strcmp(args,"l") == 0){
	if(gpars.R == 1){
	  fprintf(errfile,"sfs_code error:  option --linkage can only be used \
if number of loci > 1.  Please use option -L before --linkage\n");
	  abort();
	}
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg],"F") == 0){ /* get linkage from file */
#ifdef WINDOWS
	  fprintf(errfile,"error cannot specify linkage file on web cluster\n");
	  abort();
#else
	  FILE *linkfile;
	  if(gpars.R == 1){
	    fprintf(errfile,"please specify loci using -L before linkage\n");
	    abort();
	  }
	  arg++;
	  if((linkfile = fopen(argv[arg], "r"))==NULL){
	    fprintf(errfile,"cannot open file %s\n",argv[arg]);
	    abort();
	  }
	  arg++;
	  gpars.LinkDistType = 'p';
	  for(i=0; i<gpars.R-1; i++){
	    assert(fscanf(linkfile,"%lf",&gpars.LINK[i]) == 1);
	  }
	  fclose(linkfile);
	  continue;
#endif
	}
	gpars.LinkDistType = argv[arg][0];
	if(gpars.LinkDistType != 'p' && gpars.LinkDistType != 'g'){
	  fprintf(errfile,"sfs_code error in --linkage option: must use either \
\'p\' or \'g\', not %c\n",gpars.LinkDistType);
	  abort();
	}
	argcheck(arg, arg-2, argc, argv);
	arg++;
	gpars.LINK[0] = atof(argv[arg++]);
	if(arg >= argc || (argv[arg][0] == '-' && strcmp(argv[arg], "-1") != 0))
	  ind = 0;
	else
	  ind = 1;
	for(i=1; i<gpars.R-1; i++){
	  if(ind == 0)
	    gpars.LINK[i] = gpars.LINK[0];
	  else if(ind < i)
	    gpars.LINK[i] = gpars.LINK[i%ind];
	  else{
	    if(strcmp(argv[arg], "R") == 0){ /* start repeating */
	      gpars.LINK[i] = gpars.LINK[0];
	      arg++;
	      continue;
	    }
	    if(strcmp(argv[arg], "-1") != 0)
	      argcheck(arg, arg-i-2, argc, argv);
	    gpars.LINK[i] = atof(argv[arg++]);
	    ind++;
	  }
	  if(gpars.LinkDistType == 'g' && gpars.recMap != NULL){
	    fprintf(errfile,"SFS_CODE error:  if using a recombination map with\
 multiple loci, you must use physical distances between loci.  Use --linkage (-\
l) p <distances...>.");
	    abort();
	  }
	}
	if(arg < argc && argv[arg][0] != '-'){
	  fprintf(errfile,"sfs_code error:  Too many arguments to option\n\
--linkage.  Either 1 or R-1=%ld arguments are required.",gpars.R-1);
	  abort();
	}
      }
      else if(strcmp(args,"-mutation") == 0){
	if(gpars.SKIPSITES == NULL){
	  assert(gpars.SKIPSITES = malloc(gpars.R*sizeof(int*)));
	  for(j=0; j<gpars.R; j++){
	    assert(gpars.SKIPSITES[j] = malloc(gpars.L[j]*sizeof(int)));
	    for(i=0; i<gpars.L[j]; i++){
	      gpars.SKIPSITES[j][i] = 0;
	    }
	  }
	}
	ind = arg; /* save it for use in argcheck() */
	tmpEvolEvent.tau = -1;
	tmpEvolEvent.popi = -1;
	tmpEvolEvent.locus = -1;
	tmpEvolEvent.site = -1;
	tmpEvolEvent.gamma = 0.0;
	tmpEvolEvent.h = 1.0;
	argcheck(arg+1, ind, argc, argv);
	arg++;
	tmpEvolEvent.tau = atof(argv[arg]);
	arg++;
	if(arg<argc && strcmp(argv[arg],"P") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  if(tmpEvolEvent.popi >= gpars.NPOP){
	    fprintf(errfile,"error in --mutation, population too large, read %d\n", tmpEvolEvent.popi);
	    abort();
	  }
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"L") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.locus = atol(argv[arg]);
	  if(tmpEvolEvent.locus >= gpars.R){
	    fprintf(errfile,"error in --mutation, locus too large.  Specify number of loci (--length) before using this option\n");
	    abort();
	  }
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"S") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.site = atol(argv[arg]);
	  if(tmpEvolEvent.site >=
	     gpars.L[(tmpEvolEvent.locus>0?tmpEvolEvent.locus:0)]){
	    fprintf(errfile,"error in --mutation, site too large.  Specify locus lengths (--length) before using this option\n");
	    abort();
	  }
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"H") == 0){
	  arg++;
	  tmpEvolEvent.h = atof(argv[arg]);
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"G") == 0){
	  arg++;
	  tmpEvolEvent.gamma = atof(argv[arg]);
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg], "F") == 0){
#ifdef WINDOWS
	  fprintf(errfile,"error cannot specify trajectory file on web \
cluster\n");
	  abort();
#else
	  char c;
	  FILE *trajINfile;
	  int isInt=1;
	  arg++;
	  if((trajINfile=fopen(argv[arg],"r"))==NULL){
	    fprintf(errfile,"cannot open file %s\n",argv[--arg]);
	    abort();
	  }
	  while((c=fgetc(trajINfile))!=EOF){
	    if(c == '\n')
	      gpars.nTrajLines++;
	  }
	  fclose(trajINfile);
	  
	  if((trajINfile=fopen(argv[arg++],"r"))==NULL){
	    fprintf(errfile,"cannot open file %s\n",argv[--arg]);
	    abort();
	  }
	  assert(gpars.CONDTRAJ = malloc(gpars.nTrajLines*sizeof(*gpars.CONDTRAJ)));
	  for(i=0; i<gpars.nTrajLines; i++){
	    assert(gpars.CONDTRAJ[i] = malloc(2*sizeof(**gpars.CONDTRAJ)));
	  }
	  for(i=0; i<gpars.nTrajLines; i++){
	    if(fscanf(trajINfile, "%lf", &gpars.CONDTRAJ[i][0]) != 1){
	      fprintf(errfile,"did not see numerical value for time on line \
%ld\n", i);
	      abort();
	    }
	    if(gpars.CONDTRAJ[i][0] - (double)(int)gpars.CONDTRAJ[i][0] == 0.0){
	      isInt = 0;
	    }
	    if(fscanf(trajINfile, "%lf", &gpars.CONDTRAJ[i][1]) != 1){
	      fprintf(errfile,"did not see numerical value for frequency on line %ld for population %ld\n", i, j-1);
	      abort();
	    }
	    while((c=fgetc(trajINfile)) != EOF && (c != '\n'));//red through line
	  }
	  fclose(trajINfile);
	  if(isInt == 0){ /* assumes BACKWARDS time scaled by 4N */
	    double **TMP = NULL;
	    assert(TMP = malloc(gpars.nTrajLines*sizeof(*TMP)));
	    for(i=0; i<gpars.nTrajLines; i++){
	      assert(TMP[i] = malloc((gpars.NPOP+1)*sizeof(**TMP)));
	      TMP[i][0] = (gpars.CONDTRAJ[gpars.nTrajLines-1][0]-gpars.CONDTRAJ[gpars.nTrajLines-i-1][0])*4.0*ppars[0].N;
	      TMP[i][1] = gpars.CONDTRAJ[gpars.nTrajLines-i-1][1];
	    }
	    for(i=0; i<gpars.nTrajLines; i++){
	      gpars.CONDTRAJ[i][0] = floor(TMP[i][0]);
	      if(i>0){
		if(gpars.CONDTRAJ[i][0] <= gpars.CONDTRAJ[i-1][0]){
		  if(errPrinted==0){
		    fprintf(errfile, "error in trajectory:  generation times\
 not increasing, try increasing population size?\nconflicting lines:\n\
(%ld):\t%.0lf\n(%ld):\t%.0lf\n",i-1,floor(gpars.CONDTRAJ[i-1][0]), i, 
			    floor(gpars.CONDTRAJ[i][0]));
		    abort();
		    errPrinted++;
		  }
		}
	      }
	      gpars.CONDTRAJ[i][1] = TMP[i][1];
	      free(TMP[i]);
	    }
	    free(TMP);
	  }
	  
#ifdef VERBOSE_DEBUG
	  printf("isInt = %d\n",isInt);
	  for(i=0; i<gpars.nTrajLines; i++){
	    printf("%ld (%.0lf):\t%lf\n", i, floor(gpars.CONDTRAJ[i][0]),
		   gpars.CONDTRAJ[i][1]);
	  }
	  getchar();
#endif

#endif
	
	}
	if(tmpEvolEvent.locus == -1){
	  tmpEvolEvent.locus = (long)(ran1(&gpars.seed)*gpars.R);
	}
	if(tmpEvolEvent.site == -1){
	  tmpEvolEvent.site =
	    (long)(ran1(&gpars.seed)*gpars.L[tmpEvolEvent.locus]);
	}
	gpars.SKIPSITES[tmpEvolEvent.locus][tmpEvolEvent.site] = 1;
	tmpEvolEvent.eventType = 9;
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-trackAncestry") == 0){
	gpars.TRACKANC = 1;
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg],"a")==0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  if((ancFile=fopen(argv[arg],"a")) == 0){
	    fprintf(errfile,"cannot append to ancestry file: %s\n", argv[arg]);
	    abort();
	  }
	}
	else{
	  if((ancFile=fopen(argv[arg],"w")) == 0){
	    fprintf(errfile,"cannot write to ancestry file: %s\n", argv[arg]);
	    abort();
	  }
	}
	arg++;
      }
      else if(strcmp(args,"-trackTrajectory") == 0){
	if(toBeTracked != 0){
	  fprintf(errfile,"Can only use --trackTrajectory once.  This includes \
the track allele option in -TW\n");
	  abort();
	}
	toBeTracked = 1;
	ind = arg; /* save it for use in argcheck() */
	argcheck(arg+1, ind, argc, argv);
	arg++;
	tmpEvolEvent.tau = 0;
	tmpEvolEvent.popi = -1;
	tmpEvolEvent.locus = -1;
	tmpEvolEvent.site = -1;
	if(arg<argc && strcmp(argv[arg],"T") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.tau = atof(argv[arg]);
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"P") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"L") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.locus = atol(argv[arg]);
	  if(tmpEvolEvent.locus >= gpars.R){
	    fprintf(errfile,"error in --trackTrajectory: invalid locus %ld; must be >=0 and <%ld\n", tmpEvolEvent.locus, gpars.R);
	    abort();
	  }
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"S") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  tmpEvolEvent.site = atol(argv[arg]);
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"R") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  gpars.trajMinFreq = atof(argv[arg]);
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  gpars.trajMaxFreq = atof(argv[arg]);
#ifdef UNIT_TEST
	  printf("min-max freq = %f-%f\n",gpars.trajMinFreq, gpars.trajMaxFreq);
#endif
	  arg++;
	  if(arg<argc && strcmp(argv[arg],"S") == 0){
	    gpars.stopHitFreq = 1;
	    arg++;
	  }
	}
	if(arg<argc && strcmp(argv[arg],"A") == 0){
	  gpars.autoRestart = 1;
	  arg++;
	  if(arg<argc && strcmp(argv[arg], "G") == 0){
	    gpars.gensAfterFix = atol(argv[arg+1]);
	    arg+=2;
	  }
	}
	if(arg<argc && strcmp(argv[arg],"M") == 0){
	  argcheck(arg+1, ind, argc, argv);
	  arg++;
	  if(INITIALSEED == 0)
	    gpars.trajMaxReps = atol(argv[arg]);
	  arg++;
	}
	if(arg<argc && strcmp(argv[arg],"F") == 0){
	  argcheck(arg+1,ind,argc,argv);
	  arg++;
	  if(strcmp(argv[arg],"a") == 0){
	    arg++;
	    if(!(trajFile=fopen(argv[arg],"a"))){
	      fprintf(errfile,"cannot append to %s\n",argv[arg]);
	      abort();
	    }
	  }
	  else{
	    if(trajFile == NULL){
	      if(!(trajFile=fopen(argv[arg],"w"))){
		fprintf(errfile,"cannot write to %s\n",argv[arg]);
		abort();
	      }
	    }
	  }
	  arg++;
	}
	if(tmpEvolEvent.tau < 0){
	  /* track from beginning */
	  gpars.TRACKTRAJ = 1;
	  gpars.trajPop = tmpEvolEvent.popi;
	  gpars.trajLoc = tmpEvolEvent.locus;
	  gpars.trajSite = tmpEvolEvent.site;
	}
	else{
	  tmpEvolEvent.eventType = 7;
	  tmpEvolEvent.parIndex = 1;
	  if(tmpEvolEvent.tau > gpars.lastEvent){
	    gpars.lastEvent = tmpEvolEvent.tau;
	  }
	  addEvolEvent(tmpEvolEvent, &devents, NULL);
	}
      }
      else if(strcmp(args,"-annotate") == 0 || strcmp(args,"a") == 0){
	if(aBeforeL == -1)
	  aBeforeL = 0; /* just a dummy checker */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "F") == 0){ /* read from file */
#ifdef WINDOWS
	  fprintf(errfile,"error cannot specify annotation file on web \
cluster\n");
	  abort();
#else
	  char c, *line;
	  FILE *annfile;
	  long maxLen = 1000;
	  arg++;
	  if((annfile=fopen(argv[arg++],"r"))==NULL){
	    fprintf(errfile,"cannot open file %s\n",argv[--arg]);
	    abort();
	  }
	  assert(line = malloc(maxLen*sizeof(*line)));
	  if(fscanf(annfile, "%ld", &gpars.R) != 1){
	    fprintf(errfile,"first line of annotation file must have the \
number of loci\n");
	    abort();
	  }
	  while((c=fgetc(annfile)) != EOF && (c == ';' || isspace(c)));
	  ungetc(c,annfile);
#ifdef VERBOSE_DEBUG
	  fprintf(errfile," R = %ld\n",gpars.R);
#endif
	  if(gpars.R > 1){
	    assert(gpars.L = realloc(gpars.L,gpars.R*sizeof(*gpars.L)));
	    assert(gpars.LINK = realloc(gpars.LINK,
					(gpars.R-1)*sizeof(*gpars.LINK)));
	    assert(gpars.ANNOTATE =
		   realloc(gpars.ANNOTATE, gpars.R*sizeof(*gpars.ANNOTATE)));
	    assert(gpars.SEX = realloc(gpars.SEX,gpars.R*sizeof(*gpars.SEX)));
	    assert(mutAr = realloc(mutAr, gpars.R*sizeof(struct mutArray*)));
	    
	    for(i=0; i<gpars.NPOP; i++){
	      assert(ppars[i].h = 
		     realloc(ppars[i].h,gpars.R*sizeof(*ppars[i].h)));
	      assert(ppars[i].hProbDom = 
		     realloc(ppars[i].hProbDom,
			     gpars.R*sizeof(*ppars[i].hProbDom)));
	      assert(ppars[i].codomPop =
		     realloc(ppars[i].codomPop,
			     gpars.R*sizeof(*ppars[i].codomPop)));
	      assert(ppars[i].hDistType =
		     realloc(ppars[i].hDistType,
			     gpars.R*sizeof(*ppars[i].hDistType)));
	      assert(ppars[i].hProbGamma =
		     realloc(ppars[i].hProbGamma,
			     gpars.R*sizeof(*ppars[i].hProbGamma)));
	      assert(ppars[i].halpha =
		     realloc(ppars[i].halpha,
			     gpars.R*sizeof(*ppars[i].halpha)));
	      assert(ppars[i].hlambda =
		     realloc(ppars[i].hlambda,
			     gpars.R*sizeof(*ppars[i].hlambda)));
	      assert(ppars[i].hnormMean =
		     realloc(ppars[i].hnormMean,
			     gpars.R*sizeof(*ppars[i].hnormMean)));
	      assert(ppars[i].hnormVar =
		     realloc(ppars[i].hnormVar,
			     gpars.R*sizeof(*ppars[i].hnormVar)));

	      assert(ppars[i].GAMMA = 
		     realloc(ppars[i].GAMMA,gpars.R*sizeof(*ppars[i].GAMMA)));
	      assert(ppars[i].ProbPos = 
		     realloc(ppars[i].ProbPos,
			     gpars.R*sizeof(*ppars[i].ProbPos)));
	      assert(ppars[i].ProbDel = 
		     realloc(ppars[i].ProbDel,
			     gpars.R*sizeof(*ppars[i].ProbDel)));
	      assert(ppars[i].ProbNeut =
		     realloc(ppars[i].ProbNeut,
			     gpars.R*sizeof(*ppars[i].ProbNeut)));
	      assert(ppars[i].selDistType =
		     realloc(ppars[i].selDistType,
			     gpars.R*sizeof(*ppars[i].selDistType)));
	      assert(ppars[i].ProbNegGamma =
		     realloc(ppars[i].ProbNegGamma,
			     gpars.R*sizeof(*ppars[i].ProbNegGamma)));
	      assert(ppars[i].ProbPosGamma =
		     realloc(ppars[i].ProbPosGamma,
			     gpars.R*sizeof(*ppars[i].ProbPosGamma)));
	      assert(ppars[i].alphaN =
		     realloc(ppars[i].alphaN,
			     gpars.R*sizeof(*ppars[i].alphaN)));
	      assert(ppars[i].lambdaN =
		     realloc(ppars[i].lambdaN,
			     gpars.R*sizeof(*ppars[i].lambdaN)));
	      assert(ppars[i].alphaP =
		     realloc(ppars[i].alphaP,
			     gpars.R*sizeof(*ppars[i].alphaP)));
	      assert(ppars[i].lambdaP =
		     realloc(ppars[i].lambdaP,
			     gpars.R*sizeof(*ppars[i].lambdaP)));
	      assert(ppars[i].normMean =
		     realloc(ppars[i].normMean,
			     gpars.R*sizeof(*ppars[i].normMean)));
	      assert(ppars[i].normVar =
		     realloc(ppars[i].normVar,
			     gpars.R*sizeof(*ppars[i].normVar)));
	      assert(ppars[i].f0 =
		     realloc(ppars[i].f0,
			     gpars.R*sizeof(*ppars[i].f0)));
	    
	      for(j=1; j<gpars.R; j++){
		ppars[i].hDistType[j] = ppars[i].hDistType[0];
		ppars[i].h[j] = ppars[i].h[0];
		ppars[i].hProbDom[j] = ppars[i].hProbDom[0];
		ppars[i].hProbGamma[j] = ppars[i].hProbGamma[0];
		ppars[i].halpha[j] = ppars[i].halpha[0];
		ppars[i].hlambda[j] = ppars[i].hlambda[0];
		ppars[i].hnormMean[j] = ppars[i].hnormMean[0];
		ppars[i].hnormVar[j] = ppars[i].hnormVar[0];

		ppars[i].selDistType[j] = ppars[i].selDistType[0];
		ppars[i].GAMMA[j] = ppars[i].GAMMA[0];
		ppars[i].ProbNeut[j] = ppars[i].ProbNeut[0];
		ppars[i].ProbDel[j] = ppars[i].ProbDel[0];
		ppars[i].ProbPos[j] = ppars[i].ProbPos[0];
		ppars[i].ProbNegGamma[j] = ppars[i].ProbNegGamma[0];
		ppars[i].ProbPosGamma[j] = ppars[i].ProbPosGamma[0];
		ppars[i].alphaN[j] = ppars[i].alphaN[0];
		ppars[i].lambdaN[j] = ppars[i].lambdaN[0];
		ppars[i].alphaP[j] = ppars[i].alphaP[0];
		ppars[i].lambdaP[j] = ppars[i].lambdaP[0];
		ppars[i].normMean[j] = ppars[i].normMean[0];
		ppars[i].normVar[j] = ppars[i].normVar[0];
		ppars[i].f0[j] = ppars[i].f0[0];
		if(i == 0){
		  /* create first mutation as default */
		  assert(mutAr[j] = malloc(sizeof(struct mutArray)));
		  mutAr[j]->numMuts = 1;
		  mutAr[j]->mutIndex = 0;
		  assert(mutAr[j]->muts = malloc(sizeof(struct history*)));
		  mutAr[j]->muts[0] = popTrash();
		  assert(mutAr[j]->muts[0]->event = 
			 malloc(sizeof(struct event)));
		  mutAr[j]->muts[0]->event->site = 0;
		  assert(mutAr[j]->muts[0]->event->genFix = 
			 malloc(gpars.NPOP*sizeof(long)));
		  assert(mutAr[j]->muts[0]->event->genDead = 
			 malloc(gpars.NPOP*sizeof(long)));
		  mutAr[j]->muts[0]->event->free = '1'; /* free by default */
		  mutAr[j]->muts[0]->event->index = 0;
		  mutAr[j]->muts[0]->event->checkGen = -1; 
		  mutAr[j]->muts[0]->event->checkStep = -1; 
		  mutAr[j]->muts[0]->event->checkPop = -1;
		  mutAr[j]->muts[0]->event->nSites = 0;
		  mutAr[j]->muts[0]->event->nucs = NULL;
		  assert(mutAr[j]->muts[0]->event->fixed = 
			 malloc(gpars.NPOP*sizeof(char)));
		  assert(mutAr[j]->muts[0]->event->nextHap =
			 malloc(gpars.NPOP*sizeof(long)));
		  assert(mutAr[j]->muts[0]->event->numCarriers =
			 malloc(gpars.NPOP*sizeof(long)));
		  assert(mutAr[j]->muts[0]->event->maxHaps =
			 malloc(gpars.NPOP*sizeof(long)));
		  assert(mutAr[j]->muts[0]->event->hapFreq =
			 malloc(gpars.NPOP*sizeof(long**)));
		  gpars.SEX[j] = '0';
		  if(j < gpars.R-1)
		    gpars.LINK[j] = 0.0;
		}
		mutAr[j]->muts[0]->event->genFix[i] = 0;
		mutAr[j]->muts[0]->event->genDead[i] = -LONG_MAX;
		mutAr[j]->muts[0]->event->fixed[i] = '0'; /*free by default*/
		mutAr[j]->muts[0]->event->nextHap[i] = 0;
		mutAr[j]->muts[0]->event->numCarriers[i] = 0;
		mutAr[j]->muts[0]->event->maxHaps[i] = 1;
		assert(mutAr[j]->muts[0]->event->hapFreq[i] =
		       malloc(sizeof(long*))); 
		mutAr[j]->muts[0]->event->hapFreq[i][0] = NULL;
	      }
	    }
	  }
	  
	  popn[0].conSeq = NULL;
	  for(i=0; i<gpars.R; i++){
	    j = 0;
	    while((c=fgetc(annfile))!=EOF && 
		  (c=='-' || isdigit(c) || c=='.' || isspace(c))){/*length*/
	      line[j] = c;
	      j++;
	    }
	    line[j] = '\0';
	    if(c==';'){
	      if(i==0 || i==gpars.R){
		fprintf(errfile,"error, cannot specify inter-locus length \
before first locus or after last locus\n");
		abort();
	      }
	      if(sscanf(line, "%lf", &gpars.LINK[i-1]) != 1){
		fprintf(errfile,"did not read inter-locus distance %ld (%s)\n",
			i-1, line);
		abort();
	      }
	      while((c=fgetc(annfile))!=EOF && isspace(c));
	      if(c!=EOF)
		ungetc(c, annfile);
	      i--;
	      continue;
	    }
	    if(sscanf(line, "%ld", &gpars.L[i]) != 1){
	      fprintf(errfile,"did not read length of locus %ld (%s)\n",i,line);
	      abort();
	    }
	    if(gpars.L[i] > maxLen)
	      maxLen = gpars.L[i];
	    
	    /* next get sequence */
	    while((c=fgetc(annfile)) != EOF && isspace(c)); /*skip whitespace*/
	    if(c != ','){ /* sequence is specified! */
	      ungetc(c, annfile);
	      if(popn[0].conSeq == NULL)
		assert(popn[0].conSeq = malloc(gpars.R*sizeof(char*)));
	      assert(popn[0].conSeq[i] = malloc((gpars.L[i]+1)*sizeof(char)));
	      j = 0;
	      while((c=fgetc(annfile)) != EOF && c != ','){
		if(isspace(c))  continue;
		popn[0].conSeq[i][j] = FromCGTA(c);
#ifdef VERBOSE_DEBUG
		fprintf(errfile,"%c",popn[0].conSeq[i][j]);
#endif
		j++;
	      }
	      if(j != gpars.L[i]){
		fprintf(errfile,"error in locus %ld; expecting sequence \
length %ld, read %ld basepairs\n",i, gpars.L[i], j);
		abort();
	      }
	      popn[0].conSeq[i][j] = '\0';
	    }

	    /* now get annotation */
	    while((c=fgetc(annfile)) != EOF && isspace(c)); /*skip whitespace*/
	    if(c != 'N' && c != 'C'){
	      fprintf(errfile, "annotation improperly specified in locus %ld\n",
		      i);
	      abort();
	    }
	    gpars.ANNOTATE[i] = c;
	    if(gpars.ANNOTATE[i] == 'C' && gpars.L[i]%3 != 0){
	      fprintf(errfile,"length of all coding loci must be multiples \
of 3, locus %ld has length %ld\n",i, gpars.L[i]);
	      abort();
	    }
	    
	    /* now get selDistType information */
	    while((c=fgetc(annfile)) != EOF && (c == ',' || isspace(c)));
	    ppars[0].selDistType[i] = c-'0';
	    if(c == '0' || c == '4'); /* no parameters */
	    else{ /* get all parameters */
	      j = 0;
	      while((c=fgetc(annfile)) != EOF){
		if(c != ',')
		  line[j++] = c;
		else
		  break;
	      }
	      line[j] = '\0';

	      if(ppars[0].selDistType[i] == 1){
		if(sscanf(line,"%lf %lf %lf",&ppars[0].GAMMA[i],
			  &ppars[0].ProbPos[i], &ppars[0].ProbDel[i]) != 3){
		  fprintf(errfile,"must include 3 parameters to selDistType 1 \
in locus %ld in annotation file\n",i);
		  abort();
		}
		ppars[0].ProbNeut[i]=1-ppars[0].ProbPos[i]-ppars[0].ProbDel[i];
	      }
	      else if(ppars[0].selDistType[i] == 2){
		if(sscanf(line,"%lf %lf %lf %lf %lf",&ppars[0].ProbPosGamma[i],
			  &ppars[0].alphaP[i], &ppars[0].lambdaP[i], 
			  &ppars[0].alphaN[i], &ppars[0].lambdaN[i]) == 5){
		  ppars[0].ProbNegGamma[i] = 1-ppars[0].ProbPosGamma[i];
		}
		else if(sscanf(line,"%lf %lf %lf %lf %lf %lf", 
			       &ppars[0].ProbPosGamma[i],
			       &ppars[0].alphaP[i], &ppars[0].lambdaP[i], 
			       &ppars[0].alphaN[i], &ppars[0].lambdaN[i],
			       &ppars[0].ProbNegGamma[i]) == 6){}
		else{
		  fprintf(errfile,"must include 5/6 parameters to selDistType 2\
 in locus %ld in annotation file\n",i);
		  abort();
		}
	      }
	      else if(ppars[0].selDistType[i] == 3){
		if(sscanf(line,"%lf %lf",&ppars[0].normMean[i], 
			  &ppars[0].normVar[i]) != 2){
		  fprintf(errfile,"must include 2 parameters to selDistType 3 \
in locus %ld in annotation file\n",i);
		  abort();
		}
	      }
	      else{
		fprintf(errfile,"unrecognized selDistType=%d in locus %ld\n",
			ppars[0].selDistType[i], i);
		abort();
	      }
	    }
	    if(ppars[0].selDistType[i] != 0)
	      for(j=0; j<gpars.NPOP; j++){
		if(ppars[j].stayNeut == 1){
		  fprintf(errfile,"error: when using annotations, cannot specify a neutral population right now... contact Ryan if this is essential.\n");
		  abort();
		}
		ppars[j].neutpop = 0;
	      }
	    
	    /* now get constraint parameter */
	    while((c=fgetc(annfile)) != EOF && (c == ',' || isspace(c)));
	    ungetc(c, annfile);
	    j = 0;
	    while((c=fgetc(annfile)) != EOF && c != ';'){
	      line[j++] = c;
#ifdef VERBOSE_DEBUG
	      fprintf(errfile,"c = %c\n",c);
#endif
	    }
	    line[j] = '\0';
	    if(sscanf(line, "%f", &ppars[0].f0[i]) != 1){
	      fprintf(errfile,"did not read constraint parameter in locus \
%ld (%s)\n", i, line);
	      abort();
	    }
	    if(ppars[0].f0[i] < 1-FLT_MIN){
	      if(ppars[0].stayNeut == 1){
		fprintf(errfile,"error: when using annotations, cannot specify a neutral population right now... contact Ryan if this is essential.\n");
		abort();
	      }
	      ppars[0].neutpop = 0;
	    }

	    if(c != ';'){
	      fprintf(errfile,"locus %ld did not end in a semi-colon in \
annotation file\n", i);
	      abort();
	    }
#ifdef VERBOSE_DEBUG
	    fprintf(errfile,"locus %ld: length=%ld; sequence=%s; ann=%c; \
selDistType=%d;  constraint=%f\n",i,gpars.L[i], popn[0].conSeq[i],
		   gpars.ANNOTATE[i], ppars[0].selDistType[i], ppars[0].f0[i]);
#endif
	  } /* end loop over loci */
	  fclose(annfile);
	  free(line);
#endif
	}
	else{
	  gpars.ANNOTATE[0] = argv[arg++][0];
	  if(gpars.ANNOTATE[0] != 'C' && gpars.ANNOTATE[0] != 'N'){
	    fprintf(errfile,"sfs_code error: annotations must either be \'C\' \
or \'N\'.  Read \'%c\' for locus %d\n",gpars.ANNOTATE[0], 0);
	    abort();
	  }
	  if(arg >= argc || argv[arg][0] == '-')
	    ind = 0;
	  else
	    ind = 1;
	  for(i=1; i<gpars.R; i++){
	    if(ind == 0)
	      gpars.ANNOTATE[i] = gpars.ANNOTATE[0];
	    else if(ind < i)
	      gpars.ANNOTATE[i] = gpars.ANNOTATE[i%ind];
	    else{
	      if(strcmp(argv[arg], "R") == 0){ /* start repeating */
		gpars.ANNOTATE[i] = gpars.ANNOTATE[0];
		arg++;
		continue;
	      }
	      gpars.ANNOTATE[i] = argv[arg++][0];
	      if(gpars.ANNOTATE[i] != 'C' && gpars.ANNOTATE[i] != 'N'){
		fprintf(errfile,"sfs_code error: annotations must either be \
\'C\' or \'N\'.  Read \'%c\' for locus %ld\n",gpars.ANNOTATE[i], i);
		abort();
	      }	    
	      ind++;
	    }
	  }
	  if(arg < argc && argv[arg][0] != '-'){
	    fprintf(errfile,"sfs_code error:  Too many arguments to option \
--annotate (-a).  Either 1 or R=%ld arguments are required.  Use this option \
after -L option.\n",gpars.R);
	    abort();
	  }
	}
      }
      else if(strcmp(args,"-sex") == 0 || strcmp(args,"x") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	gpars.SEX[0] = argv[arg++][0];
	if(gpars.SEX[0] != '0' && gpars.SEX[0] != '1'){
	  fprintf(errfile,"sfs_code error: annotations must either be \'0\' \
or \'1\'.  Read \'%c\' for locus %d\n",gpars.SEX[0], 0);
	  abort();
	}
	if(arg >= argc || argv[arg][0] == '-')
	  ind = 0;
	else
	  ind = 1;
	for(i=1; i<gpars.R; i++){
	  if(ind == 0)
	    gpars.SEX[i] = gpars.SEX[0];
	  else if(ind < i)
	    gpars.SEX[i] = gpars.SEX[i%ind];
	  else{
	    if(strcmp(argv[arg], "R") == 0){ /* start repeating */
	      gpars.SEX[i] = gpars.SEX[0];
	      arg++;
	      continue;
	    }
	    gpars.SEX[i] = argv[arg++][0];
	    if(gpars.SEX[i] != '0' && gpars.SEX[i] != '1'){
	      fprintf(errfile,"sfs_code error: annotations must either be \
\'0\' or \'1\'.  Read \'%c\' for locus %ld\n",gpars.SEX[i], i);
	      abort();
	    }	    
	    ind++;
	  }
	}
	if(arg < argc && argv[arg][0] != '-'){
	  fprintf(errfile,"sfs_code error:  Too many arguments to option \
--sex (-x).  Either 1 or R=%ld arguments are required.  Use this option after \
-L option.\n",gpars.R);
	  abort();
	}
      }
      else if(strcmp(args,"-ploidy") == 0 || strcmp(args,"P") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	gpars.P = atoi(argv[arg++]);
	if(gpars.P != 1 && gpars.P != 2 && gpars.P != 4){
	  fprintf(errfile,
		  "Unfortunately, only P=1,2,4 have been implemented\n");
	  abort();
	}
	if(gpars.P != 2 && devents != NULL){
	  fprintf(errfile,"sfs_code error:  P must be defined before\
evolutionary events so that the time parameter will be correct\n");
	  abort();
	}
      }
      else if(strcmp(args,"-tetraType") == 0 || strcmp(args,"p") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	gpars.tetraType = atoi(argv[arg++]);
      }
      else if(strcmp(args,"-substMod") == 0 || strcmp(args,"M") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	ind = atoi(argv[arg]);
	gpars.substMod = ind;
	if(ind == 0 || ind == 4 || ind == 5) arg++; /* no parameters */
	else if(ind == 1){
	  argcheck(arg+1, arg, argc, argv);
	  arg++;
	  for(i=0; i<gpars.NPOP; i++)  ppars[i].PSI = atof(argv[arg]);
	  arg++;
	}
	else if(ind == 2){
	  argcheck(arg+1, arg, argc, argv);
	  arg++;
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].KAPPA = atof(argv[arg]);
	    for(j=0; j<4; j++){
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : (j==((k+2)%4) ? 
			       ppars[i].KAPPA*ppars[i].baseFreq[k] :
			       ppars[i].baseFreq[k]));
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	  arg++;
	}
	else if(ind == 3){
	  argcheck(arg+1, arg, argc, argv);
	  arg++;
	  argcheck(arg+1, arg-1, argc, argv);
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].KAPPA = atof(argv[arg]);
	    ppars[i].KAPPACpG = atof(argv[arg+1]);
	    ppars[i].PSI = atof(argv[arg+2]);
	    if(ppars[i].PSI < 0 || ppars[i].PSI > 1){
	      fprintf(errfile,"sfs_code error:  PSI must be between 0 and 1.  \
Read %f for population %d.\n",ppars[atoi(argv[arg])].PSI, atoi(argv[arg]));
	    }
	    for(j=0; j<4; j++){
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : (j==((k+2)%4) ? 
			       ppars[i].KAPPA*ppars[i].baseFreq[k] :
			       ppars[i].baseFreq[k]));
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	  arg += 3;
	}
	else if(ind == 6){ /* GTR model */
	  argcheck(arg+1, arg, argc, argv);
	  arg++;
	  for(i=0; i<gpars.NPOP; i++){
	    for(j=0; j<6; j++){
	      argcheck(arg+j, arg-1, argc, argv);
	      ppars[i].GTR[j] = atof(argv[arg+j]);
	    }
	    for(j=0; j<4; j++){
	      for(k=0; k<4; k++){
		int t = j+k-(j==0||k==0 ? 1 : 0);
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : ppars[i].GTR[t]*ppars[i].baseFreq[k]);
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	  arg += 6;
	}
      }
      else if(strcmp(args,"-INF_SITES") == 0 || strcmp(args,"I") == 0){
	gpars.INFSITES = 1;
	arg++;
      }
      else if(strcmp(args,"-seed") == 0 || strcmp(args,"s") == 0){
	if(INITIALSEED == 0){
	  argcheck(arg+1, arg, argc, argv);
	  arg++;
	  gpars.seed = atol(argv[arg++]);
	  if(gpars.seed > 0)
	    gpars.seed *= -1;
	}
	else{ /* do not reset seed after first iteration */
	  arg += 2;
	}
      }
      else if(strcmp(args,"-theta") == 0 || strcmp(args,"t") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  ppars[atoi(argv[arg])].THETA = atof(argv[arg+1]);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++)
	    ppars[i].THETA = atof(argv[arg]);
	  arg++;
	}
      }
      else if(strcmp(args,"Tt") == 0){ /* change theta */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 3;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.THETA = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-rho") == 0 || strcmp(args,"r") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	double trho = -1.0;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  ppars[atoi(argv[arg])].RHO = atof(argv[arg+1]);
	  arg += 2;
	  continue;
	}
	else if(strcmp(argv[arg], "F") == 0){ /* read rec. map file */
	  if(gpars.R > 1 && gpars.LinkDistType != 'p'){
	    fprintf(errfile,"SFS_CODE error:  if using a recombination map with\
 multiple loci, you must specify inter-locus distances in physical distance.  \
Use --linkage (-l) p <distances...>.\n");
	    abort();
	  }
#ifdef WINDOWS
	  fprintf(errfile, "error cannot specify recombination file on web \
cluster\n");
	  abort();
#else
	  long linenum=0;
	  FILE *mapfile;
	  char line[1000];
	  arg++;
	  if((mapfile = fopen(argv[arg++],"r"))==NULL){
	    fprintf(errfile,"cannot open file %s\n",argv[--arg]);
	    abort();
	  }
	  while(fgets(line, sizeof(line), mapfile) != NULL){
	    if(linenum == 0){
	      if(sscanf(line,"%ld\n", &gpars.mapPoints) != 1 &&
		 sscanf(line,"%ld,%lf;", &gpars.mapPoints,&trho) != 1){
		fprintf(errfile,"didn't read first line in mapfile\n");
		abort();
	      }
	      assert(gpars.recMap = malloc(gpars.mapPoints*sizeof(double)));
	      assert(gpars.recMapPos = malloc(gpars.mapPoints*sizeof(long)));
	    }
	    else{
	      if(linenum > gpars.mapPoints){
		fprintf(errfile,"too many lines in mapfile, expecting %ld\n",
			gpars.mapPoints);
		abort();
	      }
	      if(sscanf(line,"%ld %lf", &gpars.recMapPos[linenum-1],
			&gpars.recMap[linenum-1]) != 2){
		fprintf(errfile,"did not read mapfile line %ld correctl\n",
			linenum);
		abort();
	      }
	      if(linenum>1 && 
		 (gpars.recMapPos[linenum-1] <= gpars.recMapPos[linenum-2] ||
		  gpars.recMap[linenum-1] <= gpars.recMap[linenum-2])){
		fprintf(errfile,"recMap or positions not increasing!\n");
		fprintf(errfile,"recMapPos[%ld] = %ld;  recMapPos[%ld] = %ld\n",
			linenum-2, gpars.recMapPos[linenum-2], linenum-1,
			gpars.recMapPos[linenum-1]);
		fprintf(errfile,"recMap[%ld] = %e;  recMap[%ld] = %e\n",
			linenum-2, gpars.recMap[linenum-2], linenum-1,
			gpars.recMap[linenum-1]);
		abort();
	      }
	    }
	    linenum++;
	  }
	  if(linenum != gpars.mapPoints+1){
	    fprintf(errfile, "not enough lines in mapfile.\n");
	    fprintf(errfile, "expecting %ld, read %ld\n", gpars.mapPoints,
		    linenum-1);
	    abort();
	  }
	  if(fabs(1-gpars.recMap[gpars.mapPoints-1]) > DBL_MIN){
	    fprintf(errfile,"recombination map must be CDF (last entry = 1)\n");
	    fprintf(errfile,"read recMap[%ld] = %f\n",gpars.mapPoints-1,
		    gpars.recMap[gpars.mapPoints-1]);
	    abort();
	  }
	  fclose(mapfile);
#endif
	}
	if(arg >= argc || argv[arg][0] == '-'){
	  fprintf(errfile,"must include rho for option --rho (-r)\n");
	  abort();
	}
	for(i=0; i<gpars.NPOP; i++){
	  if(trho < 0.0)
	    ppars[i].RHO = atof(argv[arg]);
	  else
	    ppars[i].RHO = trho;
	}
	arg++;
      }
      else if(strcmp(args,"Tr") == 0){ /* change rho */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 4;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.RHO = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-pMaleRec") == 0 || strcmp(args,"Y") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	argcheck(arg+1, arg-1, argc, argv); /* population/value */
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  ppars[atoi(argv[arg])].pMaleRec = atof(argv[arg+1]);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].pMaleRec = atof(argv[arg]);
	  }
	  arg++;
	}
      }
      else if(strcmp(args,"TY") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 17;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.pMaleRec = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-geneConversion") == 0 || strcmp(args,"H") == 0){
	double bgc=-1.0;
	pop = -1;
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  arg++;
	  pop = atof(argv[arg]);
	  arg++;
	}
	if(strcmp(argv[arg], "B") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* BGC */
	  arg++;
	  bgc = atof(argv[arg]);
	  arg++;
	}
	if(pop >= 0){
	  if(bgc > -0.5)
	    ppars[pop].BGC = bgc;
	  ppars[pop].fGC = atof(argv[arg]);
	  ppars[pop].GCtract = atof(argv[arg+1]);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    if(bgc > -0.5)
	      ppars[i].BGC = bgc;
	    ppars[i].fGC = atof(argv[arg]);
	    ppars[i].GCtract = atof(argv[arg+1]);
	  }
	  arg += 2;
	}
      }
      else if(strcmp(args,"TH") == 0){ /* change gene conversion rates */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 16;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	p1 = -1;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  argcheck(arg+2, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	if(strcmp(argv[arg], "B") == 0){
	  argcheck(arg+1, arg-(p1==-1?2:4), argc, argv);
	  argcheck(arg+2, arg-(p1==-1?2:4), argc, argv);
	  argcheck(arg+3, arg-(p1==-1?2:4), argc, argv);
	  arg++;
	  tmpEvolEvent.newP.BGC = atof(argv[arg++]);
	}
	else{
	  argcheck(arg+1, arg-(p1==-1?2:4), argc, argv);
	  tmpEvolEvent.newP.BGC = 0.5;
	}
	
	tmpEvolEvent.newP.fGC = atof(argv[arg++]);
	tmpEvolEvent.newP.GCtract = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-indel") == 0 || strcmp(args,"u") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* INSRATE */
	  argcheck(arg+3, arg-1, argc, argv); /* DELRATE */
	  argcheck(arg+4, arg-1, argc, argv); /* INDELlength */
	  arg++;
	  ppars[atoi(argv[arg])].INSRATE = atof(argv[arg+1]);
	  ppars[atoi(argv[arg])].DELRATE = atof(argv[arg+2]);
	  ppars[atoi(argv[arg])].INDELlength = atof(argv[arg+3]);
	  arg += 4;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].INSRATE = atof(argv[arg]);
	    ppars[i].DELRATE = atof(argv[arg+1]);
	    ppars[i].INDELlength = atof(argv[arg+2]);
	  }
	  arg += 3;
	}
      }
      else if(strcmp(args,"Tu") == 0){ /* change indel parameters */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 13;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.INSRATE = atof(argv[arg++]);
	tmpEvolEvent.newP.DELRATE = atof(argv[arg++]);
	tmpEvolEvent.newP.INDELlength = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-longIndel") == 0 || strcmp(args,"U") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* longINSRATE */
	  argcheck(arg+3, arg-1, argc, argv); /* longDELRATE */
	  argcheck(arg+4, arg-1, argc, argv); /* longINDELlength */
	  arg++;
	  ppars[atoi(argv[arg])].longINSRATE = atof(argv[arg+1]);
	  ppars[atoi(argv[arg])].longDELRATE = atof(argv[arg+2]);
	  ppars[atoi(argv[arg])].longINDELlength = atof(argv[arg+3]);
	  arg += 4;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].longINSRATE = atof(argv[arg]);
	    ppars[i].longDELRATE = atof(argv[arg+1]);
	    ppars[i].longINDELlength = atof(argv[arg+2]);
	  }
	  arg += 3;
	}
      }
      else if(strcmp(args,"TU") == 0){ /* change longIndel parameters */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 14;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.longINSRATE = atof(argv[arg++]);
	tmpEvolEvent.newP.longDELRATE = atof(argv[arg++]);
	tmpEvolEvent.newP.longINDELlength = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-inversions") == 0 || strcmp(args,"z") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* INVRATE */
	  argcheck(arg+4, arg-1, argc, argv); /* INVlength */
	  arg++;
	  ppars[atoi(argv[arg])].INVRATE = atof(argv[arg+1]);
	  ppars[atoi(argv[arg])].INVlength = atof(argv[arg+2]);
	  arg += 3;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].INVRATE = atof(argv[arg]);
	    ppars[i].INVlength = atof(argv[arg+1]);
	  }
	  arg += 2;
	}
      }
      else if(strcmp(args,"Tz") == 0){ /* change inversion parameters */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 15;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.INVRATE = atof(argv[arg++]);
	tmpEvolEvent.newP.INVlength = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-popSize") == 0 || strcmp(args,"N") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  if(atoi(argv[arg]) == 0 && devents != NULL){
	    fprintf(errfile,"sfs_code error: The ancestral population size\
N must be defined before setting evolutionary events\n");
	    abort();
	  }
	  ind = atoi(argv[arg]);
	  if(atol(argv[arg+1]) >= LONG_MAX/gpars.P){
	    fprintf(errfile, "sfs_code error:  Maximum chromosome size is %d\
, due to memory requirements.  Please contact author if this is not \
sufficient.\n", LONG_MAX);
	    abort();
	  }
	  if(atol(argv[arg+1]) >= LONG_MAX/gpars.P){
	    fprintf(errfile, "sfs_code error:  Maximum chromosome size is %d\
, due to memory requirements.  Please contact author if this is not	\
sufficient.\n", LONG_MAX);
	    abort();
	  }
	  ppars[ind].N = atol(argv[arg+1]);
	  ppars[ind].Nt = ppars[ind].N + 0.0;
	  ppars[ind].MALES = (long)(ppars[ind].N*ppars[ind].pFEMALES+.5);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    if(atol(argv[arg]) >= LONG_MAX/gpars.P){
	      fprintf(errfile, "sfs_code error:  Maximum chromosome size is \
%d, due to memory requirements.  Please contact author if this is not \
sufficient.\n", LONG_MAX);
	      abort();
	    }
	    ppars[i].N = atol(argv[arg]);
	    ppars[i].Nt = ppars[i].N + 0.0;
	    ppars[i].MALES = (long)(ppars[i].N*ppars[i].pFEMALES+.5);
	  }
	  arg++;
	}
      }
      else if(strcmp(args,"TN") == 0){ /* set specific population size */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 0;
	tmpEvolEvent.tau = atof(argv[arg]);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  arg++;
	  if(tmpEvolEvent.popi >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot apply demographic event to \
population %d (max = %u)\n",tmpEvolEvent.popi,gpars.NPOP-1);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.popi = -1;
	}
	tmpEvolEvent.newP.Nt = atof(argv[arg]);
	if(tmpEvolEvent.newP.Nt >= LONG_MAX/gpars.P){
	  fprintf(errfile, "sfs_code error:  Maximum chromosome size is \
%d, due to memory requirements.  Please contact author if this is not \
sufficient.\n", LONG_MAX);
	  abort();
	}
	arg++;
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-sampSize") == 0 || strcmp(args, "n") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	ppars[0].SS = atoi(argv[arg]);
	if(arg == argc-1 || argv[arg+1][0] == '-'){
	  for(i=0; i<gpars.NPOP; i++)
	    ppars[i].SS = ppars[0].SS;
	  arg++;
	}
	else{
	  for(i=1; i<gpars.NPOP; i++){
	    argcheck(arg+i, arg-1, argc, argv);
	    ppars[i].SS = atoi(argv[arg+i]);
	  }
	  arg += gpars.NPOP;
	}
      }
      else if(strcmp(args,"Tn") == 0){ /* sampling times */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 19;
	tmpEvolEvent.tau = atof(argv[arg++]);
	tmpEvolEvent.nu = -1;
	tmpEvolEvent.freq = -1;
	argcheck(arg, arg-2, argc, argv);
	tmpEvolEvent.popi = -1;
	if(strcmp(argv[arg], "P") == 0){
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  arg++;
	  if(tmpEvolEvent.popi >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot apply demographic event to \
population %d (max = %u)\n",tmpEvolEvent.popi,gpars.NPOP-1);
	    abort();
	  }
	}
	if(strcmp(argv[arg], "R") == 0){
	  arg++;
	  tmpEvolEvent.freq = atof(argv[arg]);
	  arg++;
	}
	if(strcmp(argv[arg], "S") == 0){
	  arg++;
	  tmpEvolEvent.nu = atof(argv[arg]);
	  arg++;
	}
	assert(tmpEvolEvent.SS = malloc(gpars.NPOP*sizeof(*tmpEvolEvent.SS)));
	if(tmpEvolEvent.popi != -1){
	  for(i=0; i<gpars.NPOP; i++)
	    tmpEvolEvent.SS[i] = 0;
	  tmpEvolEvent.SS[tmpEvolEvent.popi] = atoi(argv[arg]);
	  arg++;
	  if(arg != argc-1 && argv[arg+1][0] != '-'){
	    fprintf(errfile,"Error in -Tn.  If you specify a population, you can only include a single sample size (all other population sample sizes will be 0).\n");
	  }
	}
	else{
	  tmpEvolEvent.SS[0] = atoi(argv[arg]);
	  if(arg == argc-1 || argv[arg+1][0] == '-'){
	    for(i=1; i<gpars.NPOP; i++)
	      tmpEvolEvent.SS[i] = tmpEvolEvent.SS[0];
	    arg++;
	  }
	  else{
	    for(i=1; i<gpars.NPOP; i++){
	      argcheck(arg+i, arg-1, argc, argv);
	      tmpEvolEvent.SS[i] = atoi(argv[arg+i]);
	    }
	    arg += gpars.NPOP;
	  }
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
	free(tmpEvolEvent.SS);
	tmpEvolEvent.SS = NULL;
	gpars.TIMECOURSE = 1;
      }
      else if(strcmp(args,"-admix") == 0 || strcmp(args,"TJ") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 8;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	tmpEvolEvent.popi = atoi(argv[arg]);
	if(tmpEvolEvent.popi >= gpars.NPOP){
	  fprintf(errfile,"sfs_code error:  population incorrectly defined in --admix command.  Note that population labels start at 0!\n");
	  abort();
	}
	popInit[tmpEvolEvent.popi]++;
	arg++;
	tmpEvolEvent.newP.Nt = atof(argv[arg]);
	argcheck(arg+1, arg-3, argc, argv);
	arg++;
	tmpEvolEvent.nAncPops = atoi(argv[arg]);
	assert(tmpEvolEvent.ancPops =
	       malloc(tmpEvolEvent.nAncPops*sizeof(int)));
	assert(tmpEvolEvent.maleFreqs =
	       malloc(tmpEvolEvent.nAncPops*sizeof(float)));
	assert(tmpEvolEvent.femaleFreqs =
	       malloc(tmpEvolEvent.nAncPops*sizeof(float)));
	for(i=0; i<tmpEvolEvent.nAncPops; i++){
	  argcheck(arg+1, arg-4-i, argc, argv);
	  arg++;
	  tmpEvolEvent.ancPops[i] = atoi(argv[arg]);
	  if(tmpEvolEvent.ancPops[i] >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error:  population incorrectly defined in --admix command.  Note that population labels start at 0!\n");
	    abort();
	  }
	}
	argcheck(arg+1, arg-4-i, argc, argv);
	arg++;
	for(i=0; i<tmpEvolEvent.nAncPops; i++){
	  argcheck(arg, arg-4-tmpEvolEvent.nAncPops-i, argc, argv);
	  tmpEvolEvent.maleFreqs[i] = atof(argv[arg]);
	  tmpEvolEvent.femaleFreqs[i] = tmpEvolEvent.maleFreqs[i];
	  arg++;
	}
	if(arg < argc && strcmp(argv[arg], "F") == 0){
	  /* get female freqs separate from male freqs */
	  arg++;
	  for(i=0; i<tmpEvolEvent.nAncPops; i++){
	    argcheck(arg, arg-6-2*tmpEvolEvent.nAncPops-i, argc, argv);
	    tmpEvolEvent.femaleFreqs[i] = atof(argv[arg]);
	    arg++;
	  }
	}
	float SUM=0.0;
	for(i=0; i<tmpEvolEvent.nAncPops; i++){
	  SUM += tmpEvolEvent.maleFreqs[i];
	}
	if(fabs(1.0-SUM) > FLT_MIN){
	  fprintf(errfile,"male ancestral frequencies do not sum to 1! Check -TJ (--admix) parameters\n");
	  fprintf(errfile,"read %d ancPops: \n", tmpEvolEvent.nAncPops);
	  for(i=0; i<tmpEvolEvent.nAncPops; i++){
	    fprintf(errfile,"ancPop[%ld] = %d; maleFreq=%f; femaleFreq=%f\n",
		    i, tmpEvolEvent.ancPops[i], tmpEvolEvent.maleFreqs[i],
		    tmpEvolEvent.femaleFreqs[i]);
	  }
	  abort();
	}
	SUM = 0.0;
	for(i=0; i<tmpEvolEvent.nAncPops; i++){
	  SUM += tmpEvolEvent.femaleFreqs[i];
	}
	if(fabs(1.0-SUM) > FLT_MIN){
	  fprintf(errfile,"female ancestral frequencies do not sum to 1! Check -TJ (--admix) parameters:\n");
	  for(i=0; i<tmpEvolEvent.nAncPops; i++){
	    fprintf(errfile,"f[%ld] = %f\n", i, tmpEvolEvent.femaleFreqs[i]);
	  }
	  abort();
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
	free(tmpEvolEvent.ancPops);
	tmpEvolEvent.ancPops = NULL;
	free(tmpEvolEvent.maleFreqs);
	tmpEvolEvent.maleFreqs = NULL;
	free(tmpEvolEvent.femaleFreqs);
	tmpEvolEvent.femaleFreqs = NULL;
      }
      else if(strcmp(args,"-migMat") == 0 || strcmp(args,"m") == 0){
	if(gpars.NPOP == 1){
	  fprintf(errfile,"migration only valid for more than one pop!\n");
	  fprintf(errfile,"please specify <npop> > 1\n");
	  abort();
	}
	argcheck(arg+1, arg, argc, argv);
	arg++;
	/* 3 options, either A, P, or L */
	if(strcmp(argv[arg],"A") == 0){ /* specify a single value */
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  if(!isdigit(argv[arg][0]) && argv[arg][0] != '.'){
	    fprintf(errfile, "error at entry %d; non-negative migration \
rate must be specified;  read %s\n",arg, argv[arg]);
	    abort();
	  }
	  for(i=0; i<gpars.NPOP; i++){
	    for(j=0; j<gpars.NPOP; j++){
	      if(i == j)  continue;
	      gpars.mig_mat[i][j] = atof(argv[arg])/(gpars.NPOP-1);
	    }
	  }
	  arg++;
	}
	else if(strcmp(argv[arg],"L") == 0){ /* specify entry to/from all pops*/
	  arg++;
	  p1=0;
	  for(i=0; i<gpars.NPOP; i++){
	    for(j=0; j<gpars.NPOP; j++){
	      if(i == j)  continue;
	      p1++;
	      argcheck(arg, arg-1-p1, argc, argv);
	      if(!isdigit(argv[arg][0]) && argv[arg][0] != '.'){
		fprintf(errfile, "error at entry %d  (M_[%ld, %ld]);  ",arg,i,j);
		fprintf(errfile,"non-negative migration rate must be\
specified;  read %s\n",argv[arg]);
		abort();
	      }
	      gpars.mig_mat[i][j] = atof(argv[arg++]);
	    }
	  }
	}
	else if(strcmp(argv[arg],"P") == 0){/*specify entry for specific comb.*/
	  arg++;
	  argcheck(arg, arg-2, argc, argv);   /* P1 */
	  argcheck(arg+1, arg-3, argc, argv); /* P2 */
	  argcheck(arg+2, arg-4, argc, argv); /* Mij */
	  if(!isdigit(argv[arg][0]) || !isdigit(argv[arg+1][0]) ||
	     (!isdigit(argv[arg+2][0]) && argv[arg+2][0] != '.')){
	    fprintf(errfile,"error:  for -m P, must specify two populations \
and a non-negative migration rate; read -m %s %s %s\n",argv[arg], argv[arg+1],
		    argv[arg+2]);
	    abort();
	  }
	  p1 = atoi(argv[arg++]);
	  p2 = atoi(argv[arg++]);
	  if(!isdigit(argv[arg][0]) && argv[arg][0] != '.'){
	    fprintf(errfile, "error at entry %d  (M_[%d, %d]);  ",arg,p1,p2);
	    fprintf(errfile,"non-negative migration rate must be\
specified;  read %s\n",argv[arg]);
	    abort();
	  }
	  gpars.mig_mat[p1][p2] = atof(argv[arg++]);
	}
	else{
	  fprintf(errfile,"error at element %d; first argument to migration \
option must be \'A\', \'L\', or \'P\', read %s\n",arg,argv[arg]);
	  abort();
	}
      }
      else if(strcmp(args,"Tm") == 0){
	tmpEvolEvent.eventType = 7;
	tmpEvolEvent.parIndex = 0;
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.tau = atof(argv[arg++]);
	if(tmpEvolEvent.newG.mig_mat == NULL){
	  assert(tmpEvolEvent.newG.mig_mat =
		 malloc(gpars.NPOP*sizeof(double *)));
	  for(i=0; i<gpars.NPOP; i++)
	    assert(tmpEvolEvent.newG.mig_mat[i] =
		   malloc(gpars.NPOP*sizeof(double)));
	}
	if(strcmp(argv[arg],"A") == 0){ /* specify a single value */
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  tmpEvolEvent.popi = -1;
	  tmpEvolEvent.popj = -1;
	  if(!isdigit(argv[arg][0]) && argv[arg][0] != '.'){
	    fprintf(errfile, "error at entry %d; non-negative migration \
rate must be specified;  read %s\n",arg, argv[arg]);
	    abort();
	  }
	  for(i=0; i<gpars.NPOP; i++){
	    for(j=0; j<gpars.NPOP; j++){
	      if(i == j)  continue;
	      tmpEvolEvent.newG.mig_mat[i][j] = atof(argv[arg])/(gpars.NPOP-1);
	    }
	  }
	  arg++;
	}
	else if(strcmp(argv[arg],"L") == 0){ /* specify entry to/from all pops*/
	  arg++;
	  tmpEvolEvent.popi = -1;
	  tmpEvolEvent.popj = -1;
	  p1=0;
	  for(i=0; i<gpars.NPOP; i++){
	    for(j=0; j<gpars.NPOP; j++){
	      if(i == j)  continue;
	      p1++;
	      argcheck(arg, arg-1-p1, argc, argv);
	      if(!isdigit(argv[arg][0]) && argv[arg][0] != '.'){
		fprintf(errfile, "error at entry %d  (M_[%ld, %ld]); ",arg,i,j);
		fprintf(errfile,"non-negative migration rate must be\
specified;  read %s\n",argv[arg]);
		abort();
	      }
	      tmpEvolEvent.newG.mig_mat[i][j] = atof(argv[arg++]);
	    }
	  }
	}
	else if(strcmp(argv[arg],"P") == 0){/*specify entry for specific comb.*/
	  arg++;
	  argcheck(arg, arg-2, argc, argv);   /* P1 */
	  argcheck(arg+1, arg-3, argc, argv); /* P2 */
	  argcheck(arg+2, arg-4, argc, argv); /* Mij */
	  if(!isdigit(argv[arg][0]) || !isdigit(argv[arg+1][0]) ||
	     (!isdigit(argv[arg+2][0]) && argv[arg+2][0] != '.')){
	    fprintf(errfile,"error:  for -m P, must specify two populations \
and a non-negative migration rate; read -m %s %s %s\n",argv[arg], argv[arg+1],
		    argv[arg+2]);
	    abort();
	  }
	  p1 = atoi(argv[arg++]);
	  p2 = atoi(argv[arg++]);
	  if(!isdigit(argv[arg][0]) && argv[arg][0] != '.'){
	    fprintf(errfile, "error at entry %d  (M_[%d, %d]);  ",arg,p1,p2);
	    fprintf(errfile,"non-negative migration rate must be\
specified;  read %s\n",argv[arg]);
	    abort();
	  }
	  tmpEvolEvent.popi = p1;
	  tmpEvolEvent.popj = p2;
	  tmpEvolEvent.newG.mig_mat[p1][p2] = atof(argv[arg++]);
	}
	else{
	  fprintf(errfile,"error at element %d; first argument to migration \
option must be \'A\', \'L\', or \'P\', read %s\n",arg,argv[arg]);
	  abort();
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
	for(i=0; i<gpars.NPOP; i++)
	  free(tmpEvolEvent.newG.mig_mat[i]);
	free(tmpEvolEvent.newG.mig_mat);
	tmpEvolEvent.newG.mig_mat = NULL;
      }
      else if(strcmp(args,"-pMaleMig") == 0 || strcmp(args,"y") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	argcheck(arg+1, arg-1, argc, argv); /* population/value */
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  ppars[atoi(argv[arg])].pMaleMig = atof(argv[arg+1]);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].pMaleMig = atof(argv[arg]);
	  }
	  arg++;
	}
      }
      else if(strcmp(args,"Ty") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 11;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.pMaleMig = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-propFemale") == 0 || strcmp(args,"f") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  ppars[atoi(argv[arg])].pFEMALES = atof(argv[arg+1]);
	  if(fabs(ppars[atoi(argv[arg])].pMaleMig-.5) <= DBL_MIN)
	    ppars[atoi(argv[arg])].pMaleMig = 1-ppars[atoi(argv[arg])].pFEMALES;
	  ppars[atoi(argv[arg])].MALES = 
	    (long)(ppars[atoi(argv[arg])].N*ppars[atoi(argv[arg])].pFEMALES);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].pFEMALES = atof(argv[arg]);
	    if(fabs(ppars[i].pMaleMig-0.5) <= DBL_MIN)
	      ppars[i].pMaleMig = 1-ppars[i].pFEMALES;
	    ppars[i].MALES = (long)(ppars[i].N*ppars[i].pFEMALES);
	  }
	  arg++;
	}
      }
      else if(strcmp(args,"-self") == 0 || strcmp(args,"i") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv); /* population */
	  argcheck(arg+2, arg-1, argc, argv); /* value */
	  arg++;
	  ppars[atoi(argv[arg])].SELF = atof(argv[arg+1]);
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++)
	    ppars[i].SELF = atof(argv[arg]);
	  arg++;
	}
      }
      else if(strcmp(args,"Ti") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 5;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.SELF = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-KAPPA") == 0 || strcmp(args,"K") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(args,"P") == 0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  i = atoi(argv[arg]);
	  ppars[i].KAPPA = atof(argv[arg+1]);
	  for(j=0; j<4; j++){
	    for(k=0; k<4; k++){
	      ppars[i].transProb[j][k] = 
		(j==k ? 0 : (j==((k+2)%4) ? 
			     ppars[i].KAPPA*ppars[i].baseFreq[k] :
			     ppars[i].baseFreq[k]));
	      ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	    }
	    for(k=0; k<4; k++){
	      ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	    }
	  }
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].KAPPA = atof(argv[arg]);
	    for(j=0; j<4; j++){
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : (j==((k+2)%4) ? 
			       ppars[i].KAPPA*ppars[i].baseFreq[k] :
			       ppars[i].baseFreq[k]));
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	  arg++;
	}
      }
      else if(strcmp(args,"TK") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 8;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.KAPPA = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-PSI") == 0 || strcmp(args,"C") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(args,"P") == 0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  ppars[atoi(argv[arg])].PSI = atof(argv[arg+1]);
	  if(ppars[atoi(argv[arg])].PSI < 0 || ppars[atoi(argv[arg])].PSI > 1){
	    fprintf(errfile,"sfs_code error:  PSI must be between 0 and 1.  \
Read %f for population %d.\n",ppars[atoi(argv[arg])].PSI, atoi(argv[arg]));
	    abort();
	  }
	  arg += 2;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].PSI = atof(argv[arg]);
	    if(ppars[i].PSI < 0 || ppars[i].PSI > 1){
	      fprintf(errfile,"sfs_code error:  PSI must be between 0 and 1.  \
Read %f for population %ld.\n",ppars[i].PSI, i);
	      abort();
	    }
	  }
	  arg++;
	}
      }
      else if(strcmp(args,"TC") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 9;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.PSI = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-baseFreq") == 0 || strcmp(args,"q") == 0){
	for(i=1; i<=4; i++){
	  argcheck(arg+i, arg, argc, argv);
	}
	arg++;
	if(strcmp(argv[arg],"P") == 0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  i = atoi(argv[arg]);
	  double tmp = 0.0;
	  for(j=1; j<=4; j++){
	    ppars[i].baseFreq[j-1] = atof(argv[arg+j]);
	    tmp += ppars[i].baseFreq[j-1];
	  }
	  if(fabs(tmp-1.0) > FLT_MIN){
	    fprintf(errfile,"error, base frequencies do not sum to 1\n");
	    abort();
	  }
	  for(j=0; j<4; j++){
	    for(k=0; k<4; k++){
	      if(gpars.substMod == 6){
		int t = j+k-(j==0||k==0 ? 1 : 0);
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : ppars[i].GTR[t]*ppars[i].baseFreq[k]);
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      else{
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : (j==((k+2)%4) ? 
			       ppars[i].KAPPA*ppars[i].baseFreq[k] :
			       ppars[i].baseFreq[k]));
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	    }
	    for(k=0; k<4; k++){
	      ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	    }
	  }
	  arg += 5;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    double tmp = 0.0;
	    for(j=0; j<4; j++){
	      ppars[i].baseFreq[j] = atof(argv[arg+j]);
	      tmp += ppars[i].baseFreq[j];
	    }
	    if(fabs(tmp-1.0) > FLT_MIN){
	      fprintf(errfile,"error, base frequencies do not sum to 1\n");
	      abort();
	    }
	    for(j=0; j<4; j++){
	      tmp = 0.0;
	      for(k=0; k<4; k++){
		if(gpars.substMod == 6){
		  int t = j+k-(j==0||k==0 ? 1 : 0);
		  ppars[i].transProb[j][k] =
		    (j==k ? 0 : ppars[i].GTR[t]*ppars[i].baseFreq[k]);
		}
		else{
		  ppars[i].transProb[j][k] = 
		    (j==k ? 0 : (j==((k+2)%4) ? 
				 ppars[i].KAPPA*ppars[i].baseFreq[k] :
				 ppars[i].baseFreq[k]));
		}
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	  arg += 4;
	}
      }
      else if(strcmp(args,"Tq") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 18;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	double tmp = 0.0;
	for(i=0; i<4; i++){
	  tmpEvolEvent.newP.baseFreq[i] = atof(argv[arg++]);
	  tmp += tmpEvolEvent.newP.baseFreq[i];
	}
	if(fabs(tmp-1.0) > FLT_MIN){
	  fprintf(errfile,"error, base frequencies to not sum to 1.\n");
	  abort();
	}
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"-rateClassSites") == 0 || strcmp(args,"V") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  argcheck(arg+1, arg-2, argc, argv);
	  argcheck(arg+2, arg-2, argc, argv);
	  ppars[atoi(argv[arg])].RateClassSites = atoi(argv[arg+1]);
	  ppars[atoi(argv[arg])].RateParamSites = atof(argv[arg+2]);
	  arg += 3;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].RateClassSites = atoi(argv[arg]);
	    ppars[i].RateParamSites = atof(argv[arg+1]);
	  }
	  arg += 2;
	}
      }
      else if(strcmp(args,"-rateClassLoci") == 0 || strcmp(args,"v") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	if(strcmp(argv[arg], "L") == 0){
	  if(gpars.R == 1){
	    fprintf(errfile,"please specify loci using -L before rateClassLoci\n");
	    abort();
	  }
	  /* defining rate classes by locus; initialize to locus length */
	  for(i=0; i<gpars.NPOP; i++){
	    if(ppars[i].RateClassLoci == 1){
	      assert(ppars[i].RateParamLoci = realloc(ppars[i].RateParamLoci, 
						      gpars.R*sizeof(double)));
	      for(j=0; j<gpars.R; j++)
		ppars[i].RateParamLoci[j] = gpars.L[j];
	    }
	    ppars[i].RateClassLoci = 0;
	    if(strcmp(argv[arg+1], "A") == 0){
	      for(j=0; j<gpars.R; j++){
		ppars[i].RateParamLoci[j] = atof(argv[arg+2]);
	      }
	    }
	    else{
	      ppars[i].RateParamLoci[atoi(argv[arg+1])] = atof(argv[arg+2]);
	    }
	  }
	  arg += 3;
	}
	else if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-1, argc, argv);
	  arg++;
	  argcheck(arg+1, arg-2, argc, argv);
	  argcheck(arg+2, arg-2, argc, argv);
	  ppars[atoi(argv[arg])].RateClassLoci = atoi(argv[arg+1]);
	  ppars[atoi(argv[arg])].RateParamLoci[0] = atof(argv[arg+2]);
	  arg += 3;
	}
	else{
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].RateClassLoci = atoi(argv[arg]);
	    ppars[i].RateParamLoci[0] = atof(argv[arg+1]);
	  }
	  arg += 2;
	}
      }
      else if(strcmp(args,"-GenEffect") == 0 || strcmp(args,"G") == 0){
	if(gpars.NPOP == 1){
	  fprintf(errfile,"sfs_code error:  Generation effect only valid when \
more than one population is simulated.\n");
	  abort();
	}
	argcheck(arg+1, arg, argc, argv);
	arg++;
	pop = atoi(argv[arg]);
	arg++;
	ppars[pop].GenEffect = atoi(argv[arg++]);
	if(!(ppars[pop].GenEffect >= 1 || ppars[pop].GenEffect <= -2)){
	  fprintf(errfile,"sfs_code error:  Generation time effect must be \
>= 1 or <= -2.\n");
	  abort();
	}
      }
      else if(strcmp(args,"TG") == 0){
	argcheck(arg+1, arg, argc, argv);
 	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 12;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  argcheck(arg+1, arg-2, argc, argv);
	  arg++;
	  p1 = tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg+1, arg-3, argc, argv);
	  arg++;
	  if(p1 >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot set parameters for population %d \
(max = %u)\n",p1,gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	tmpEvolEvent.newP.GenEffect = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"Td") == 0){ /* discrete demographic event */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 3;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  argcheck(arg, arg-3, argc, argv);
	  arg++;
	  if(tmpEvolEvent.popi >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot apply demographic event to \
population %d (max = %u)\n",tmpEvolEvent.popi,gpars.NPOP-1);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.popi = -1;
	}
	tmpEvolEvent.nu = atof(argv[arg]);
	arg++;
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"Tg") == 0){ /* set exponential growth rate */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 1;
	tmpEvolEvent.tau = atof(argv[arg]);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  arg++;
	  if(tmpEvolEvent.popi >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot apply demographic event to \
population %d (max = %u)\n",tmpEvolEvent.popi,gpars.NPOP-1);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.popi = -1;
	}
	tmpEvolEvent.newP.popAlpha = atof(argv[arg]);
	arg++;
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"Tk") == 0){ /* set logistic growth rate */
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 6;
	tmpEvolEvent.parIndex = 2;
	tmpEvolEvent.tau = atof(argv[arg]);
	arg++;
	if(strcmp(argv[arg], "P") == 0){
	  arg++;
	  tmpEvolEvent.popi = atoi(argv[arg]);
	  arg++;
	  if(tmpEvolEvent.popi >= gpars.NPOP){
	    fprintf(errfile,"sfs_code error: cannot apply demographic event to \
population %d (max = %u)\n",tmpEvolEvent.popi,gpars.NPOP-1);
	    abort();
	  }
	}
	else{
	  tmpEvolEvent.popi = -1;
	}
	tmpEvolEvent.newP.K = atol(argv[arg++]);
	tmpEvolEvent.newP.popAlpha = atof(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"TE") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 5;
	tmpEvolEvent.site = 0;
	tmpEvolEvent.tau = atof(argv[arg++]);
	if(arg < argc && isdigit(argv[arg][0])){
	  tmpEvolEvent.popi = atoi(argv[arg++]);
	  if(tmpEvolEvent.popi > gpars.NPOP){
	    fprintf(errfile,"sfs_code error (-TE):  terminating population \
that does not exist!  pop = %d;  max = %u\n", tmpEvolEvent.popi, gpars.NPOP-1);
	    abort();
	  }
	}
	else
	  tmpEvolEvent.popi = -1;
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"TS") == 0){
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 0;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	tmpEvolEvent.popi = atoi(argv[arg]);
	argcheck(arg+1, arg-2, argc, argv);
	arg++;
	tmpEvolEvent.popj = atoi(argv[arg]);
	popInit[tmpEvolEvent.popj]++;

	if(tmpEvolEvent.popi >= gpars.NPOP || tmpEvolEvent.popj >= gpars.NPOP){
	  fprintf(errfile,"sfs_code error:  NPOP incorrectly defined.  Note \
that population labels start at 0!\n");
	  abort();
	}
	arg++;
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else if(strcmp(args,"TD") == 0){ /* domestication event */
	gpars.DOMESTICATE = 1;
	argcheck(arg+1, arg, argc, argv);
	arg++;
	tmpEvolEvent.eventType = 1;
	tmpEvolEvent.locus = -1;
	tmpEvolEvent.tau = atof(argv[arg]);
	argcheck(arg+1, arg-1, argc, argv);
	arg++;
	tmpEvolEvent.popi = atoi(argv[arg]);
	argcheck(arg+1, arg-2, argc, argv);
	arg++;
	tmpEvolEvent.popj = atoi(argv[arg]);
	popInit[tmpEvolEvent.popj]++;
	argcheck(arg+1, arg-3, argc, argv);
	arg++;
	tmpEvolEvent.freq = atof(argv[arg]);
	argcheck(arg+1, arg-4, argc, argv);
	arg++;
	tmpEvolEvent.newP.Nt = atof(argv[arg]);
	if(tmpEvolEvent.newP.Nt >= LONG_MAX/gpars.P){
	  fprintf(errfile, "sfs_code error:  Maximum population size is \
%d, due to memory requirements.  Please contact author if this is not\
sufficient.\n", LONG_MAX);
	  abort();
	}	
	arg++;
	if(arg < argc && argv[arg][0] != '-')
	  tmpEvolEvent.locus = atol(argv[arg++]);
	if(tmpEvolEvent.tau > gpars.lastEvent){
	  gpars.lastEvent = tmpEvolEvent.tau;
	}
	addEvolEvent(tmpEvolEvent, &devents, NULL);
      }
      else{
	fprintf(errfile,"my apologies. -%s has not been implemented.\n",args);
	exit(-7);
      }
    }
    else{
      fprintf(errfile,"sfs_code error at argument %d/%d  (%s)\n",arg,argc,args);
      abort();
    }
  }

  /* CHECK PARAMETERS... maybe move this to a function */
  for(i=1; i<gpars.NPOP; i++){
    if(popInit[i] == 0){
      fprintf(errfile,"pop %ld must be initialized using -TS option\n",i);
      abort();
    }
    else if(popInit[i] > 1){
      fprintf(errfile,"pop %ld was initialized multiple times using -TS \
or -TD option.  A population can only be initialized once.\n",i);
      abort();
    }
  }
  free(popInit);
}

