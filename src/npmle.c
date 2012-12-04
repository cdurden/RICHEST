#include <R.h>
#include <pari/pari.h>
#include <Rinternals.h>
GEN PoissonMixtureEMSpeciesPerCount(GEN MixingF0, GEN tSpeciesPerCount, int nCounts, int nSupportPoints, int nSpecies) {
  GEN MixingF1 = zerovec(2);
  GEN Pr = zerovec(nCounts);
  GEN PrNorm = zerovec(nCounts);
  GEN PrNorm_j;
  GEN vPi1_i,vLambda1_i;
  pari_sp ltop1;
  int i,j;

  gel(MixingF1,1) = zerovec(nSupportPoints);
  gel(MixingF1,2) = zerovec(nSupportPoints);
  for(j=1;j<=nCounts;j++) {
    gel(Pr,j) = zerovec(nSupportPoints);
    for(i=1;i<=nSupportPoints;i++) {
      ltop1 = avma;
      // generate Pr(i|SpeciesCount_j,vPi0,vLambda0) -- gel(gel(Pr,j),i)
      gel(gel(Pr,j),i) = gerepileupto(ltop1, // construct numerator
        gmul(
          gel(gel(MixingF0,2),i),
          // Pr(vSpeciesCounts_j|i,vPi0,vLambda0)
          gdiv(
            gmul(
              gexp(gsub(gen_0,gel(gel(MixingF0,1),i)), DEFAULTPREC),
              powgi(gel(gel(MixingF0,1),i),gel(gel(tSpeciesPerCount,1),j))
            ),
            mpfactr(itos(gel(gel(tSpeciesPerCount,1),j)),DEFAULTPREC)
          )
        )
      );
    }
    PrNorm_j = gen_0;
    for(i=1;i<=nSupportPoints;i++) { // calc. sum of Pr(i|...
      PrNorm_j = gadd(PrNorm_j,gel(gel(Pr,j),i));
    }
    gel(PrNorm,j) = PrNorm_j;
    for(i=1;i<=nSupportPoints;i++) { // calc. sum of Pr(i|...
      gel(gel(Pr,j),i) = gdiv(gel(gel(Pr,j),i),gel(PrNorm,j));
    }
  }
  for(i=1;i<=nSupportPoints;i++) { // calc. new parameters
    vPi1_i = gen_0;
    vLambda1_i = gen_0;
    for(j=1;j<=nCounts;j++) {
      vPi1_i = gadd(
        vPi1_i,
        gmul(
          gel(gel(tSpeciesPerCount,2),j),
          gel(gel(Pr,j),i)
        )
      );
      vLambda1_i = gadd(
        vLambda1_i,
        gmul(
          gmul(
            gel(gel(tSpeciesPerCount,1),j),
            gel(gel(tSpeciesPerCount,2),j)
          ),
          gel(gel(Pr,j),i)
        )
      );
    }
    vLambda1_i = gdiv(vLambda1_i,vPi1_i);
    vPi1_i = gdiv(vPi1_i,stoi(nSpecies));
    gel(gel(MixingF1,1),i)=vLambda1_i;
    gel(gel(MixingF1,2),i)=vPi1_i;
  }
  return(MixingF1);
}
GEN PoissonMixtureEM(GEN MixingF0, GEN vSpeciesCounts, int nSupportPoints, int nSpecies) {
  GEN MixingF1 = zerovec(2);
  GEN Pr = zerovec(nSpecies);
  GEN PrNorm = zerovec(nSpecies);
  GEN PrNorm_j;
  GEN vPi1_i,vLambda1_i;
  pari_sp ltop1;
  int i,j;

  gel(MixingF1,1) = zerovec(nSupportPoints);
  gel(MixingF1,2) = zerovec(nSupportPoints);
  for(j=1;j<=nSpecies;j++) {
    gel(Pr,j) = zerovec(nSupportPoints);
    for(i=1;i<=nSupportPoints;i++) {
      ltop1 = avma;
      // generate Pr(i|vSpeciesCounts_j,vPi0,vLambda0) -- gel(gel(Pr,j),i)
      gel(gel(Pr,j),i) = gerepileupto(ltop1, // construct numerator
        gmul(
          gel(gel(MixingF0,2),i),
          // Pr(vSpeciesCounts_j|i,vPi0,vLambda0)
          gdiv(
            gmul(
              gexp(gsub(gen_0,gel(gel(MixingF0,1),i)), DEFAULTPREC),
              powgi(gel(gel(MixingF0,1),i),gel(vSpeciesCounts,j))
            ),
            mpfactr(itos(gel(vSpeciesCounts,j)),DEFAULTPREC)
          )
          //
        )
      );
    }
    PrNorm_j = gen_0;
    for(i=1;i<=nSupportPoints;i++) { // calc. sum of Pr(i|...
      PrNorm_j = gadd(PrNorm_j,gel(gel(Pr,j),i));
    }
    gel(PrNorm,j) = PrNorm_j;
    for(i=1;i<=nSupportPoints;i++) { // calc. sum of Pr(i|...
      gel(gel(Pr,j),i) = gdiv(gel(gel(Pr,j),i),gel(PrNorm,j));
    }
  }
  for(i=1;i<=nSupportPoints;i++) { // calc. new parameters
    vPi1_i = gen_0;
    vLambda1_i = gen_0;
    for(j=1;j<=nSpecies;j++) {
      vPi1_i = gadd(vPi1_i,gel(gel(Pr,j),i));
      vLambda1_i = gadd(vLambda1_i,gmul(gel(vSpeciesCounts,j),gel(gel(Pr,j),i)));
    }
    vLambda1_i = gdiv(vLambda1_i,vPi1_i);
    vPi1_i = gdiv(vPi1_i,stoi(nSpecies));
    gel(gel(MixingF1,1),i)=vLambda1_i;
    gel(gel(MixingF1,2),i)=vPi1_i;
  }
  return(MixingF1);
}
GEN CalcLikelihood(GEN vSpeciesCounts, GEN vLambda, GEN vPi, int nSupportPoints, int nSpecies) {
pari_sp ltop00,ltop0,ltop1;
int i,j;
GEN LikelihoodInnerSum;
GEN Likelihood = gen_0;
  ltop00 = avma;
  ltop0 = avma;
  for(j=1;j<=nSpecies;j++) {
    LikelihoodInnerSum = gen_0;
    for(i=1;i<=nSupportPoints;i++) {
      ltop1 = avma;
      LikelihoodInnerSum = gerepileupto(ltop1,
        gadd(
          LikelihoodInnerSum,
          gdiv(
            gmul(
              gmul(
                gel(vPi,i),
                gexp(gsub(gen_0,gel(vLambda,i)), DEFAULTPREC)
              ),
              powgi(gel(vLambda,i),gel(vSpeciesCounts,j))
            ),
            mpfactr(itos(gel(vSpeciesCounts,j)),DEFAULTPREC)
          )
        )
      );
    }
    Likelihood = gerepileupto(ltop0, gadd(Likelihood, glog(LikelihoodInnerSum, DEFAULTPREC)));
  }
  Likelihood = gerepileupto(ltop00, Likelihood);
  return(Likelihood);
}
GEN CalcLikelihoodSpeciesPerCount(GEN tSpeciesPerCount, GEN vLambda, GEN vPi, int nCounts, int nSupportPoints) {
pari_sp ltop0,ltop1;
int i,j;
GEN LikelihoodInnerSum;
GEN Likelihood = gen_0;
  ltop0 = avma;
  for(j=1;j<=nCounts;j++) {
    LikelihoodInnerSum = gen_0;
    for(i=1;i<=nSupportPoints;i++) {
      ltop1 = avma;
      LikelihoodInnerSum = gerepileupto(ltop1,
        gadd(
          LikelihoodInnerSum,
            gdiv(
              gmul(
                gmul(
                  gel(vPi,i),
                  gexp(gsub(gen_0,gel(vLambda,i)), DEFAULTPREC)
                ),
                powgi(gel(vLambda,i),gel(gel(tSpeciesPerCount,1),j))
              ),
              mpfactr(itos(gel(gel(tSpeciesPerCount,1),j)),DEFAULTPREC)
            )
        )
      );
    }
    Likelihood = gerepileupto(ltop0, gadd(Likelihood, gmul(gel(gel(tSpeciesPerCount,2),j),glog(LikelihoodInnerSum, DEFAULTPREC))));
  }
  return(Likelihood);
}

GEN FindMLE(GEN vSpeciesCounts, int nSupportPoints, int nSpecies) {
  if(nSupportPoints > nSpecies) {
    printf("Error: number of support points is greater than the number of species.\n");
    return(0);
  }
  int i, j, it;
  GEN L0 = gen_0;
  GEN L1 = gen_0;
  pari_sp ltop0,lbot;
  ltop0 = avma;
  GEN MixingF = zerovec(2); // vLambda = gel(MixingF,1), vPi = gel(MixingF,2)
  GEN Recip_n = gdiv(gen_1,stoi(nSupportPoints));
  gel(MixingF,1) = zerovec(nSupportPoints);
  gel(MixingF,2) = zerovec(nSupportPoints);
  for (i=1;i<=nSupportPoints;i++) {
    gel(gel(MixingF,1),i) = gel(vSpeciesCounts,i);
    gel(gel(MixingF,2),i) = Recip_n;
  }
  it = 0;
  while(gcmp(gsub(L1,L0),gsub(gen_0,gmul(dbltor(0.00001),L1)))>=0) {
    L0 = L1;
    for(i=1;i<=2;i++) {
      for(j=1;j<=nSupportPoints;j++) {
        printf("%f ",gtodouble(gel(gel(MixingF,i),j)));
      }
      printf("\n");
    }
    MixingF = PoissonMixtureEM(MixingF,vSpeciesCounts,nSupportPoints,nSpecies);
    L1 = CalcLikelihood(vSpeciesCounts,gel(MixingF,1),gel(MixingF,2),nSupportPoints,nSpecies);
    if(it==0) {
      L0 = gsub(L1,gen_1);
    }
    printf("%i: %f\n", it, gtodouble(L1));
    it++;
  }
  lbot = avma;
  MixingF = gerepile(ltop0, lbot, MixingF);
  return(MixingF);
}
GEN FindMLESpeciesPerCount(GEN tSpeciesPerCount, int nCounts, int nSupportPoints, int nSpecies) {
  if(nSupportPoints > nCounts) {
    printf("Error: number of support points is greater than the number of observed species abundances.\n");
    return(gen_0);
  }
  int i, it;
  pari_sp ltop0,lbot;
  GEN L0 = gen_0;
  GEN L1 = gen_0;
  GEN MixingF = zerovec(2); // vLambda = gel(MixingF,1), vPi = gel(MixingF,2)
  GEN Recip_n = gdiv(gen_1,stoi(nSupportPoints));
  gel(MixingF,1) = zerovec(nSupportPoints);
  gel(MixingF,2) = zerovec(nSupportPoints);
  pari_sp limit = stack_lim(avma,1);
  ltop0 = avma;
  for (i=1;i<=nSupportPoints;i++) {
    gel(gel(MixingF,1),i) = gel(gel(tSpeciesPerCount,1),i);
    gel(gel(MixingF,2),i) = Recip_n;
  }
  it = 0;
  while(gcmp(gsub(L1,L0),gsub(gen_0,gmul(dbltor(0.00000001),L1)))>=0 && it <= 500) {
    L0 = L1;
    if(it==0) {
      L0 = CalcLikelihoodSpeciesPerCount(tSpeciesPerCount,gel(MixingF,1),gel(MixingF,2),nCounts,nSupportPoints);
      printf("%i: %f\n", it, gtodouble(L0));
    }
    if (avma < limit) gerepileall(ltop0,2,&L0,&MixingF);
    MixingF = PoissonMixtureEMSpeciesPerCount(MixingF,tSpeciesPerCount,nCounts,nSupportPoints,nSpecies);
    L1 = CalcLikelihoodSpeciesPerCount(tSpeciesPerCount,gel(MixingF,1),gel(MixingF,2),nCounts,nSupportPoints);
    it++;
    printf("%i: %f\n", it, gtodouble(L1));
  }
  lbot = avma;
  MixingF = gerepile(ltop0, lbot, MixingF);
  return(MixingF);
}

SEXP R_FindMLESpeciesPerCount(SEXP R_sSpeciesPerCount1, SEXP R_sSpeciesPerCount2, SEXP R_nCounts, SEXP R_nSupportPoints, SEXP R_nSpecies) {
  int* nCounts = INTEGER(R_nCounts);
  int* nSupportPoints = INTEGER(R_nSupportPoints);
  int* nSpecies = INTEGER(R_nSpecies);
  char* sSpeciesPerCount1 = (char*) CHAR(STRING_ELT(R_sSpeciesPerCount1,0));
  char* sSpeciesPerCount2 = (char*) CHAR(STRING_ELT(R_sSpeciesPerCount2,0));
  SEXP R_MixingF;
  int i;
  GEN MixingF;
  GEN tSpeciesPerCount;

//  printf("%i\n",*nSpecies);
//  printf("%i\n",*nCounts);
//  printf("%i\n",*nSupportPoints);
//  printf("%s\n",sSpeciesPerCount1);
//  printf("%s\n",sSpeciesPerCount2);
  tSpeciesPerCount = zerovec(2);
  gel(tSpeciesPerCount,1) = gp_read_str(sSpeciesPerCount1);
  gel(tSpeciesPerCount,2) = gp_read_str(sSpeciesPerCount2);
  MixingF = FindMLESpeciesPerCount(tSpeciesPerCount,*nCounts,*nSupportPoints,*nSpecies);
  R_MixingF = R_NilValue;
  if(gcmp0(MixingF)) {
  } else {
    PROTECT(R_MixingF = allocVector(REALSXP,2*(*nSupportPoints)));
    for(i=1;i<=*nSupportPoints;i++) {
      REAL(R_MixingF)[i-1] = gtodouble(gel(gel(MixingF,1),i));
      REAL(R_MixingF)[*nSupportPoints+i-1] = gtodouble(gel(gel(MixingF,2),i));
    }
    UNPROTECT(1);
  }
  return(R_MixingF);
}
SEXP R_FindMLE(SEXP R_sSpeciesCounts, SEXP R_nSupportPoints, SEXP R_nSpecies) {
  int* nSpecies = INTEGER(R_nSpecies);
  int* nSupportPoints = INTEGER(R_nSupportPoints);
  char* sSpeciesCounts = (char*) CHAR(STRING_ELT(R_sSpeciesCounts,0));
  SEXP R_MixingF;
  int i;
  GEN MixingF;

//  printf("%i\n",*nSpecies);
//  printf("%s\n",sSpeciesCounts);
  GEN vSpeciesCounts = gp_read_str(sSpeciesCounts);
  MixingF = FindMLE(vSpeciesCounts,*nSupportPoints,*nSpecies);
  R_MixingF = R_NilValue;
  if(MixingF==0) {
  } else {
    PROTECT(R_MixingF = allocVector(REALSXP,2*(*nSupportPoints)));
    for(i=1;i<=*nSupportPoints;i++) {
      REAL(R_MixingF)[i-1] = gtodouble(gel(gel(MixingF,1),i));
      REAL(R_MixingF)[*nSupportPoints+i-1] = gtodouble(gel(gel(MixingF,2),i));
    }
    UNPROTECT(1);
  }
  return(R_MixingF); 
}

SEXP init_pari() {
  pari_init(1000000,0);
  return(R_NilValue);
}
SEXP close_pari() {
  pari_close();
  return(R_NilValue);
}
