/* Only scoreRankGenePairsLimitedNodesC() and mytspC(), which are used in modelfinal.R, are corrected for ranks to be numeric instead of integers!! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#define EPSIL 1E-16

int index_global; /*index_global records how many (disjoint) genepairs there are when max_nodes are limited*/

struct genepair
{
  int gene1;
  int gene2;
  double score1;
  double score2;
  struct genepair *next;
};

struct genepair_bidir
{
  int gene1;
  int gene2;
  double score1;
  double score2;
  struct genepair_bidir *prev;
  struct genepair_bidir *next;
};


SEXP scoreGenePairsC2R(SEXP class1, SEXP class2, SEXP score_tie);
SEXP rankGenePairsC2R(SEXP Rmat, SEXP Rmax_pairs);
int disjointQualify(int num, int *vec, int len);
SEXP scoreRankGenePairsC2R(SEXP Rclass1, SEXP Rclass2, SEXP Rmax_pairs, SEXP Rdisjoint);
struct genepair *scoreRankGenePairsC(int class1[], int class2[], int n1, int n2, int p, int max_pairs, int disjoint);
int disjointQualifyLinkedlist(struct genepair *ptr_advance, struct genepair *head);
SEXP makeCombnC2R(SEXP Rp, SEXP Rscore);
SEXP naiveCombnC2R(SEXP Rp);
SEXP scoreRankGenePairsLimitedNodesC2R(SEXP Rclass1, SEXP Rclass2, SEXP Rmax_pairs, SEXP Rdisjoint, SEXP Rmax_nodes);
struct genepair_bidir *scoreRankGenePairsLimitedNodesC(double class1[], double class2[], int n1, int n2, int p, int max_pairs, int disjoint, int max_nodes);
int disjointQualifyLinkedlist_bidir(struct genepair_bidir *ptr_advance, struct genepair_bidir *head);
SEXP allPairsMajorityVotesC2R(SEXP Rtrain, SEXP Rgrp, SEXP Rtest);



SEXP scoreGenePairsC2R(SEXP class1, SEXP class2, SEXP score_tie)
{ 
  /*
  	class_i is the (integer) rank matrix
  */
  int i,j,k,temp,n1,n2,p,nb_pairs,nb_scores=1,index=0;
  double m1,m2,r1,r2,score1=0,score2=0;
  
  SEXP Rdim1,Rdim2;
  SEXP Rresult;
  
  class1 = coerceVector(class1,INTSXP);
  class2 = coerceVector(class2,INTSXP);
  score_tie = coerceVector(score_tie,LGLSXP);
  
  Rdim1 = getAttrib(class1, R_DimSymbol);
  Rdim2 = getAttrib(class2, R_DimSymbol);
  
  p = INTEGER(Rdim1)[1];
  n1 = INTEGER(Rdim1)[0];
  n2 = INTEGER(Rdim2)[0];
  
  nb_pairs = p*(p-1)/2;
  if(LOGICAL(score_tie)[0] == 1) {
  	nb_scores = 2;
  }
  PROTECT(Rresult=allocMatrix(REALSXP,nb_pairs,nb_scores));
  
  for (i=0; i<p; i++) {
	  for (j=0; j<i; j++) {
	  	
		  /*
		  	For each pair of genes in default order, compute both main and secondary signed scores
		  */
		
		  m1 = 0;
		  r1 = 0;
		  for (k=0; k<n1; k++) {
			  temp = INTEGER(class1)[k+i*n1] - INTEGER(class1)[k+j*n1];
			  if (temp<0) {
				  m1 = m1 + 1;
			  }
			  if (LOGICAL(score_tie)[0] == 1) {
			  	r1 = r1 + temp;
			  }
		  }
		  
		  
		  m2 = 0;
		  r2 = 0;
		  for (k=0; k<n2; k++) {
			  temp = INTEGER(class2)[k+i*n2] - INTEGER(class2)[k+j*n2];
			  if (temp<0) {
				  m2 = m2 + 1;
			  }
			  if (LOGICAL(score_tie)[0] == 1) {
			  	r2 = r2 + temp;
			  }
		  }
		  
		  
		  /*
		  	By default ordering of gene pairs, add scores to result matrix
		  */

		  /*
		  	index = (2*p - 1 - i) * i / 2 + j - i - 1;
		  */
		  
		  score1 = m1/n1 - m2/n2;
		  REAL(Rresult)[index] = score1;
		  
		  if (LOGICAL(score_tie)[0] == 1) {
			score2 = r1/n1 - r2/n2;
		  	REAL(Rresult)[index + nb_pairs] = score2;
		  }
		  
		  index++;
	  }
    
    // Rprintf("Running for gene1idx %d\n",i+1);
  }
  
  UNPROTECT(1);
  return(Rresult);
  
}



SEXP rankGenePairsC2R(SEXP Rmat, SEXP Rmax_pairs)
{
	/*
		Rmat is a 2-col matrix with rows wrt decreasing scores
	*/
	
	int index=0,k=0,nb_pairs,i,j;
	
	SEXP Rdim,Rresult;
	
	Rmat = coerceVector(Rmat, INTSXP);
	Rmax_pairs = coerceVector(Rmax_pairs, INTSXP);
	
    Rdim = getAttrib(Rmat, R_DimSymbol);
    nb_pairs = INTEGER(Rdim)[0];
	
    PROTECT(Rresult=allocMatrix(INTSXP,INTEGER(Rmax_pairs)[0],2));
	
	
	
	/*
		HOW TO INITIALIZE ENTRIES AUTOMATICALLY?
	*/
	
	for(k=0;k<INTEGER(Rmax_pairs)[0];k++)
	{
		INTEGER(Rresult)[k] = 0;
		INTEGER(Rresult)[k+INTEGER(Rmax_pairs)[0]] = 0;
	}
	
	for(k=0,index=0; index<INTEGER(Rmax_pairs)[0] && k<nb_pairs; k++)
	{
		i = INTEGER(Rmat)[k];
		j = INTEGER(Rmat)[k+nb_pairs];
		
		if(disjointQualify(i,INTEGER(Rresult),INTEGER(Rmax_pairs)[0])==1 || disjointQualify(j,INTEGER(Rresult),INTEGER(Rmax_pairs)[0])==1)
		{
			continue;
		}
		
		INTEGER(Rresult)[index] = i;
		INTEGER(Rresult)[index+INTEGER(Rmax_pairs)[0]] = j;
		index++;
    
    // Rprintf("The %d-th gene pair found\n",index);
		
	}
	
	
	
	
	UNPROTECT(1);
	
	return(Rresult);
}



int disjointQualify(int num, int *vec, int len)
{
	int temp,res=0;
	
	for(temp=0; vec[temp] >0 && temp<len; temp++)
	{
		if(num==vec[temp] || num==vec[temp+len])
		{
			res = 1;
			break;
		}
	}
	
	return(res);
}



SEXP scoreRankGenePairsC2R(SEXP Rclass1, SEXP Rclass2, SEXP Rmax_pairs, SEXP Rdisjoint)
{
  
  int p,n1,n2,k;
  struct genepair *res, *ptr;
	
	SEXP Rdim1,Rdim2,Rresult;
	
	
	Rclass1 = coerceVector(Rclass1,INTSXP);
  Rclass2 = coerceVector(Rclass2,INTSXP);
	Rmax_pairs = coerceVector(Rmax_pairs,INTSXP);
	Rdisjoint = coerceVector(Rdisjoint,LGLSXP);
	
	
	Rdim1 = getAttrib(Rclass1, R_DimSymbol);
	Rdim2 = getAttrib(Rclass2, R_DimSymbol);
	p = INTEGER(Rdim1)[1];
	n1 = INTEGER(Rdim1)[0];
	n2 = INTEGER(Rdim2)[0];
	
  
  res = scoreRankGenePairsC(INTEGER(Rclass1),INTEGER(Rclass2),n1,n2,p,INTEGER(Rmax_pairs)[0],LOGICAL(Rdisjoint)[0]);
	
  
  PROTECT(Rresult=allocMatrix(INTSXP,INTEGER(Rmax_pairs)[0],2));
  
  /* Initialization... */
  for(k=0;k<INTEGER(Rmax_pairs)[0];k++)
	{
		INTEGER(Rresult)[k] = 0;
		INTEGER(Rresult)[k+INTEGER(Rmax_pairs)[0]] = 0;
	}
  
  /* Assignment... */
  for(ptr=res->next,k=0; ptr!=0 && k<INTEGER(Rmax_pairs)[0]; res=ptr,ptr=res->next,k++)
  {
    INTEGER(Rresult)[k] = ptr->gene1;
    INTEGER(Rresult)[k+INTEGER(Rmax_pairs)[0]] = ptr->gene2;
    free(res);
  }
	
	UNPROTECT(1);
	
  
	return(Rresult);
}



struct genepair *scoreRankGenePairsC(int class1[], int class2[], int n1, int n2, int p, int max_pairs, int disjoint)
{
    /*
		This function integrates and does everything in C with linked list structure 
    where argument descriptions are seen from above
    while it always uses secondary scores to break ties
    
    This function is super slow and can be further ameliorated by,
    instead of building up the complete linked list of length p*(p-1)/2 and then remaining the top scoring pairs,
    treating differently disjoint==T/F
    a) disjoint==T: an almost complete list is still needed but if the new node with BOTH geneidx identical to some existing nodes in the list 
                    but still a lower score, there's no need to add it in the list; preferably the order of calling i/j for the nested loop can be randomized
                    ??? Q: would there be an UPPER BOUND on the number of pairs needed in the linked list to get a fixed number (max_pairs) of disjoint pairs?
    b) disjoint==F: cf. scoreRankGenePairsLimitedNodesC() below with limited nodes in the (bi-directional) linked list
    */
    
  int i,j,k=0,index=0,temp;
  double m1,m2,r1,r2,score1,score2;
  
  struct genepair *head, *ptr_new, *ptr, *ptr_advance;
  head = (struct genepair *)malloc(sizeof(struct genepair));
  head->gene1=0;head->gene2=0;head->score1=0;head->score2=0;head->next=0;
  
  
  for (i=0; i<p; i++) {
    for (j=0; j<i; j++) {
	  	
		  /*
		  	For each pair of genes in default order, compute both primary and secondary scores to break ties
		  */
      
		  m1 = 0;
      r1 = 0;
		  for (k=0; k<n1; k++) 
      {
        temp = class1[k+i*n1] - class1[k+j*n1];
        r1 = r1 + temp;
			  if (temp < 0) 
        {
				  m1 = m1 + 1;
			  }
		  }
		  
		  m2 = 0;
      r2 = 0;
		  for (k=0; k<n2; k++) 
      {
        temp = class2[k+i*n2] - class2[k+j*n2];
        r2 = r2 + temp;
			  if (temp < 0) 
        {
				  m2 = m2 + 1;
			  }
		  }
		  
  	  score1 = m1/n1 - m2/n2;
      score2 = r1/n1 - r2/n2;
      
      
      /*make new node*/
      
      ptr_new = (struct genepair *)malloc(sizeof(struct genepair));
      
      ptr_new->score1 = fabs(score1);
      ptr_new->score2 = fabs(score2);
      ptr_new->next = 0;
      if(score1>0)
      {
        ptr_new->gene1 = i+1;
        ptr_new->gene2 = j+1;
      }
      else
      {
        ptr_new->gene1 = j+1;
        ptr_new->gene2 = i+1;
      }
      
      
      /* insertion for one-directional sorted linked list*/
      for(ptr=head, ptr_advance=ptr->next; (ptr_advance!=0) && ((ptr_new->score1 < ptr_advance->score1) || ((fabs(ptr_new->score1 - ptr_advance->score1)<EPSIL) && (ptr_new->score2 < ptr_advance->score2)) || ((fabs(ptr_new->score1 - ptr_advance->score1)<EPSIL) && (fabs(ptr_new->score2 - ptr_advance->score2)<EPSIL))); ptr=ptr_advance, ptr_advance=ptr->next);
      ptr->next = ptr_new;
      ptr_new->next = ptr_advance;
      
	  }
    
    // Rprintf("Running for gene1idx %d\n",i+1);
  }
  
  
  
  /* filter out top scoring pairs*/
  if(disjoint == 1)
  {
    for(ptr=head, ptr_advance=ptr->next, index=0; ptr_advance != 0;)
    {
      if(index == max_pairs || disjointQualifyLinkedlist(ptr_advance,head)==1)
      {
        ptr->next = ptr_advance->next;
        free(ptr_advance);
        ptr_advance = ptr->next;
      }
      else
      {
        index++;
        ptr = ptr_advance;
        ptr_advance = ptr->next;
      }
      
    }
    
  }
  else
  {
    for(ptr=head, ptr_advance=ptr->next, index=0; index<max_pairs; ptr=ptr_advance, ptr_advance=ptr->next,index++);
    ptr->next = 0;
    while(ptr_advance!=0)
    {
      ptr = ptr_advance;
      ptr_advance = ptr->next;
      free(ptr);
    }
    
  }
  
  return(head);
}



int disjointQualifyLinkedlist(struct genepair *ptr_advance, struct genepair *head)
{
  struct genepair *ptr_temp;
  int ind=0;
  
  for(ptr_temp=head->next; ptr_temp!=ptr_advance; ptr_temp=ptr_temp->next)
  {
    if(ptr_advance->gene1 == ptr_temp->gene1 || ptr_advance->gene1 == ptr_temp->gene2 || ptr_advance->gene2 == ptr_temp->gene1 || ptr_advance->gene2 == ptr_temp->gene2)
    {
      ind = 1;
      break;
    }
  }
  
  return(ind);
}



SEXP makeCombnC2R(SEXP Rp, SEXP Rscore)
{
  
  /*make all combinations or 1:p with default nested order*/
  
  int i,j,index=0,l;
  
  SEXP Rresult;
  
  Rp = coerceVector(Rp, INTSXP);
  Rscore = coerceVector(Rscore, REALSXP);
  l = (INTEGER(Rp)[0])*(INTEGER(Rp)[0]-1)/2;
  
  PROTECT(Rresult=allocMatrix(INTSXP,l,2));
  
  index = 0;
  for(i=0;i<INTEGER(Rp)[0];i++)
  {
    for(j=(i+1);j<INTEGER(Rp)[0];j++)
    {
      if(REAL(Rscore)[index]>0)
      {
        INTEGER(Rresult)[index] = i+1;
        INTEGER(Rresult)[index+l] = j+1;
      }
      else
      {
        INTEGER(Rresult)[index] = j+1;
        INTEGER(Rresult)[index+l] = i+1;
      }
      
      index++;
    }
  }
  
  UNPROTECT(1);
  
  return(Rresult);
  
}



SEXP naiveCombnC2R(SEXP Rp)
{
    int i,j,index=0,l;
    
    SEXP Rresult;
    
    Rp = coerceVector(Rp,INTSXP);
    l = (INTEGER(Rp)[0])*(INTEGER(Rp)[0]-1)/2;
    
    PROTECT(Rresult=allocMatrix(INTSXP,l,2));
    
    index = 0;
    for(i=0;i<INTEGER(Rp)[0];i++)
    {
        for(j=(i+1);j<INTEGER(Rp)[0];j++)
        {
            INTEGER(Rresult)[index] = i+1;
            INTEGER(Rresult)[index+l] = j+1;
            index++;
        }
    }
    
    UNPROTECT(1);
    
    return(Rresult);
}



SEXP scoreRankGenePairsLimitedNodesC2R(SEXP Rclass1, SEXP Rclass2, SEXP Rmax_pairs, SEXP Rdisjoint, SEXP Rmax_nodes)
{
  
  int p,n1,n2,k;
  struct genepair_bidir *res, *ptr;
  
	SEXP Rdim1,Rdim2,Rresult;
	
	Rclass1 = coerceVector(Rclass1,REALSXP);
  Rclass2 = coerceVector(Rclass2,REALSXP);
	Rmax_nodes = coerceVector(Rmax_nodes,INTSXP);
  Rmax_pairs = coerceVector(Rmax_pairs,INTSXP);
	Rdisjoint = coerceVector(Rdisjoint,LGLSXP);
	
	Rdim1 = getAttrib(Rclass1, R_DimSymbol);
	Rdim2 = getAttrib(Rclass2, R_DimSymbol);
	p = INTEGER(Rdim1)[1];
	n1 = INTEGER(Rdim1)[0];
	n2 = INTEGER(Rdim2)[0];
  
  /* max_nodes is the max number of genepairs to keep track of within the process */
  res = scoreRankGenePairsLimitedNodesC(REAL(Rclass1),REAL(Rclass2),n1,n2,p,INTEGER(Rmax_pairs)[0],LOGICAL(Rdisjoint)[0],INTEGER(Rmax_nodes)[0]);
  
  
  PROTECT(Rresult=allocMatrix(INTSXP,index_global,2));
  
  /* Initialization... */
  for(k=0;k<index_global;k++)
	{
		INTEGER(Rresult)[k] = 0;
		INTEGER(Rresult)[k+index_global] = 0;
	}
  
  /* Assignment... */
  for(ptr=res->next,k=0; ptr!=0 && k<index_global; res=ptr,ptr=res->next,k++)
  {
    INTEGER(Rresult)[k] = ptr->gene1;
    INTEGER(Rresult)[k+index_global] = ptr->gene2;
    free(res);
  }
	
	UNPROTECT(1);
	
  
	return(Rresult);
}



struct genepair_bidir *scoreRankGenePairsLimitedNodesC(double class1[], double class2[], int n1, int n2, int p, int max_pairs, int disjoint, int max_nodes)
{
  /*
    In this script, we run everything with (non-circled) bi-directional linked list
    max_nodes is different from max_pairs in that
    a) for disjoint==F, it's indeed the same
    b) for disjoint==T, it's rather the number of nodes to build within a linked list in the process of running than the number of genepairs to finally output
  */
  
  int i,j,k,index=0,sgn;
  double m1,m2,r1,r2,score1,score2;
  
  struct genepair_bidir *head, *tail, *ptr_new, *ptr, *ptr_advance;
  head = (struct genepair_bidir *)malloc(sizeof(struct genepair_bidir));
  head->gene1=0;head->gene2=0;head->score1=0;head->score2=0;head->prev=0;head->next=0;
  tail = head;
  index_global=0;
  
  /*if disjoint==F, max_pairs=max_nodes=min(max_pairs, max_nodes)*/
  if(disjoint == 0)
  {
    if(max_pairs < max_nodes)
    {
      max_nodes = max_pairs;
      index_global = max_pairs;
    }
    else
    {
      max_pairs = max_nodes;
      index_global = max_nodes;
    }
  }
  
  for (i=0; i<p; i++) {
    for (j=0; j<i; j++) {
    	
		  /*
		  	For each pair of genes in default order, compute both primary and secondary scores to break ties
		  */
      
		  m1 = 0;
      r1 = 0;
		  for (k=0; k<n1; k++) 
      {
        r1 = r1 + (class1[k+i*n1] - class1[k+j*n1]);
			  if (class1[k+i*n1] < class1[k+j*n1])
        {
				  m1 = m1 + 1;
			  }
		  }
		  
		  m2 = 0;
      r2 = 0;
		  for (k=0; k<n2; k++) 
      {
        r2 = r2 + (class2[k+i*n2] - class2[k+j*n2]);
			  if (class2[k+i*n2] < class2[k+j*n2])
        {
				  m2 = m2 + 1;
			  }
		  }
		  
      sgn = m1/n1 > m2/n2;
  	  score1 = fabs(m1/n1 - m2/n2);
      score2 = fabs(r1/n1 - r2/n2);
      
      
      /*if qualified, make this new node*/
      
      if(score1 > (tail->score1) || (fabs(score1 - tail->score1)<EPSIL && score2 > (tail->score2)) || index < max_nodes)
      {
        index++; /*records how many potential nodes are made during the process, possibly being deleted later*/
        
        ptr_new = (struct genepair_bidir *)malloc(sizeof(struct genepair_bidir));
        ptr_new->score1 = score1;
        ptr_new->score2 = score2;
        ptr_new->prev = 0;
        ptr_new->next = 0;
        if(sgn == 1)
        {
          ptr_new->gene1 = i+1;
          ptr_new->gene2 = j+1;
        }
        else if(sgn == 0)
        {
          ptr_new->gene1 = j+1;
          ptr_new->gene2 = i+1;
        }
        
        
        /* insert into a bi-directional sorted linked list*/
        for(ptr=head, ptr_advance=ptr->next; (ptr_advance!=0) && ((ptr_new->score1 < ptr_advance->score1) || ((fabs(ptr_new->score1 - ptr_advance->score1)<EPSIL) && (ptr_new->score2 < ptr_advance->score2)) || ((fabs(ptr_new->score1 - ptr_advance->score1)<EPSIL) && (fabs(ptr_new->score2 - ptr_advance->score2)<EPSIL))); ptr=ptr_advance, ptr_advance=ptr->next);
        ptr->next = ptr_new;
        ptr_new->prev = ptr;
        ptr_new->next = ptr_advance;
        if(ptr_advance != 0)
        {
          ptr_advance->prev = ptr_new;
        }
        else
        {
          tail = ptr_new; /*deal with tail pointer*/
        }
        
        
        /*free tail node and redirect tail pointer*/
        if(index > max_nodes) /*NOTE that index already self-added before!!!*/
        {
          tail = tail->prev;
          free(tail->next);
          tail->next = 0;
        }
        
        /*
        for(ptr=head,k=0;ptr->next!=0;ptr=ptr->next,k++);
        Rprintf("%d\t%d\t",index,k); // k denotes nb of nodes in the linked list
        for(ptr=head,k=0;ptr!=tail;ptr=ptr->next,k++);
        Rprintf("%d\n",k); // k denotes nb of nodes until tail pointer
        */
      }
      
	  }
    
    // Rprintf("Running for gene1idx %d\n",i+1); /*printf in R console*/
  }
  
  
  
  /* filter out disjoint pairs and/or take max_pairs genepairs*/
  if(disjoint == 1)
  {
    for(ptr=head, ptr_advance=ptr->next, index_global=0; index_global<max_pairs && ptr_advance!=0;)
    {
      if(disjointQualifyLinkedlist_bidir(ptr_advance,head)==1)
      {
        ptr->next = ptr_advance->next;
        free(ptr_advance);
        ptr_advance = ptr->next;
      }
      else
      {
        index_global++; /*records how many disjoint nodes are found*/
        ptr = ptr_advance;
        ptr_advance = ptr->next;
      }
      
    }
    
    /* free all nodes left in the list*/
    ptr->next = 0;
    while(ptr_advance!=0)
    {
      ptr = ptr_advance;
      ptr_advance = ptr->next;
      free(ptr);
    }
  }
  
  return(head);
}



int disjointQualifyLinkedlist_bidir(struct genepair_bidir *ptr_advance, struct genepair_bidir *head)
{
  struct genepair_bidir *ptr_temp;
  int ind=0;
  
  for(ptr_temp=head->next; ptr_temp!=ptr_advance; ptr_temp=ptr_temp->next)
  {
    if(ptr_advance->gene1 == ptr_temp->gene1 || ptr_advance->gene1 == ptr_temp->gene2 || ptr_advance->gene2 == ptr_temp->gene1 || ptr_advance->gene2 == ptr_temp->gene2)
    {
      ind = 1;
      break;
    }
  }
  
  return(ind);
}



SEXP allPairsMajorityVotesC2R(SEXP Rtrain, SEXP Rgrp, SEXP Rtest)
{
  /*
    Rtrain (Rtest) denotes matrices of ntr*p (nts*p)
    Rgrp denotes a class vector labeled 1 and 2
  */
  
  int i,j,k,n1,n2,g1,g2;
  int ntr,nts,p,p_temp; /* p_temp controls the number of pairs used for prediction to be odd*/
  double m1,m2;
  
  SEXP Rresult;
  
  Rtrain = coerceVector(Rtrain,REALSXP);
  Rtest = coerceVector(Rtest,REALSXP);
  Rgrp = coerceVector(Rgrp,INTSXP);
  
  /*
  SEXP Rdimtr,Rdimts;
  Rdimtr = getAttrib(Rtrain, R_DimSymbol);
  Rdimts = getAttrib(Rtest, R_DimSymbol);
  p = INTEGER(Rdimtr)[1];
  ntr = INTEGER(Rdimtr)[0];
  nts = INTEGER(Rdimts)[0];
  */
  
  ntr = LENGTH(Rgrp);
  p = LENGTH(Rtrain)/ntr;
  nts = LENGTH(Rtest)/p;
  
  
  for(n1=0,n2=0,k=0; k<ntr; k++)
  {
    if(INTEGER(Rgrp)[k] == 1)
    {
      n1 = n1 + 1;
    }
    else if(INTEGER(Rgrp)[k] == 2)
    {
      n2 = n2 + 1;
    }
  }
  
  /* make sure odd number of candidate pairs are guaranteed, otherwise throw the last pair*/
  p_temp = p;
  if( (p % 2 == 0 && (p/2) % 2 == 0) || ((p-1) % 2 == 0 && ((p-1)/2) % 2 == 0) )
  {
    p_temp = p-1;
  }
  
  PROTECT(Rresult=allocVector(INTSXP,nts));
  
  /* Initialization ... */
  for(k=0; k<nts; k++)
  {
    INTEGER(Rresult)[k] = 0;
  }
  
  /* Train and Test ... */
  for (i=0; i<p_temp-1; i++) { /* careful with upper bound index!!! T.T */
	  for (j=(i+1); j<p; j++) {
      
		  for (m1=0,m2=0,k=0; k<ntr; k++) {
        if(INTEGER(Rgrp)[k] == 1)
        {
          m1 += (REAL(Rtrain)[k+i*ntr] < REAL(Rtrain)[k+j*ntr]);
        }
        else if(INTEGER(Rgrp)[k] == 2)
        {
          m2 += (REAL(Rtrain)[k+i*ntr] < REAL(Rtrain)[k+j*ntr]);
        }
		  }
		  
      if(m1/n1 > m2/n2)
      {
        g1 = i;
        g2 = j;
      }
      else
      {
        g1 = j;
        g2 = i;
      }
      
      /* count class votes for each test sample */
      for(k=0; k<nts; k++)
      {
        if(REAL(Rtest)[k+g1*nts] < REAL(Rtest)[k+g2*nts])
        {
          INTEGER(Rresult)[k] = INTEGER(Rresult)[k] + 1;
        }
        else
        {
          INTEGER(Rresult)[k] = INTEGER(Rresult)[k] - 1;
        }
      }
      
	  }
    
    // Rprintf("Running for gene1idx %d\n",i+1); /*printf in R console*/
  }
  
  /* convert to class prediction */
  for(k=0; k<nts; k++)
  {
    if(INTEGER(Rresult)[k] > 0)
    {
      INTEGER(Rresult)[k] = 1;
    }
    else
    {
      INTEGER(Rresult)[k] = 2;
    }
  }
  
  UNPROTECT(1);
  return(Rresult);
  
}






///////////////////////////////////
///////////////////////////////////
//
//      AVL tree-based TSP       //
//
///////////////////////////////////
///////////////////////////////////


#include "tspavl.h"



struct genepair_avl
{
    int gene1;
    int gene2;
    double score1;
    double score2;
};

int compare_genepair_avl(const void *pa, const void *pb, void *param);
struct avl_table *mytspC(double class1[], double class2[], int n1, int n2, int p, int max_nodes);
SEXP mytspC2R(SEXP Rclass1, SEXP Rclass2, SEXP Rmax_nodes);

int compare_genepair_avl(const void *pa, const void *pb, void *param)
{
    const struct genepair_avl *g1 = pa;
    const struct genepair_avl *g2 = pb;
    
    if ((g1->score1 > g2->score1) || ((g1->score1 == g2->score1)&&(g1->score2 > g2->score2)))
        return +1;
    else if ((g1->score1 < g2->score1) || ((g1->score1 == g2->score1)&&(g1->score2 < g2->score2)))
        return -1;
    else if(g1->gene1 != g2->gene1)
        return (g1->gene1 < g2->gene1)? -1 : +1;
    else if(g1->gene2 != g2->gene2)
        return (g1->gene2 < g2->gene2)? -1 : +1;
    else
        return 0;
}

struct avl_table *mytspC(double class1[], double class2[], int n1, int n2, int p, int max_nodes)
{
    int i,j,k;
    double m1,m2,r1,r2;
    struct avl_table *tree;
    struct avl_traverser *trav;
    struct genepair_avl *newgenepair_avl, *mingenepair_avl;
    
    tree = avl_create (compare_genepair_avl, NULL,NULL);
    trav = (struct avl_traverser *) malloc (sizeof(struct avl_traverser));
    avl_t_init (trav, tree);
    
    for (i=0; i<p; i++)
    {
        for (j=0; j<i; j++)
        {
            
            /*
             For each pair of genes, compute primary and secondary score to order
             */
            
            m1 = 0;
            r1 = 0;
            for (k=0; k<n1; k++)
            {
                r1 = r1 + (class1[k+i*n1] - class1[k+j*n1]);
                if (class1[k+i*n1] < class1[k+j*n1])
                {
                    m1 = m1 + 1;
                }
            }
            
            m2 = 0;
            r2 = 0;
            for (k=0; k<n2; k++)
            {
                r2 = r2 + (class2[k+i*n2] - class2[k+j*n2]);
                if (class2[k+i*n2] < class2[k+j*n2])
                {
                    m2 = m2 + 1;
                }
            }
            
            /* Make new pair */
            newgenepair_avl = (struct genepair_avl *)malloc(sizeof(struct genepair_avl));
            if (m1/n1 > m2/n2)
            {
                newgenepair_avl->gene1 = i+1;
                newgenepair_avl->gene2 = j+1;
            }
            else
            {
                newgenepair_avl->gene1 = j+1;
                newgenepair_avl->gene2 = i+1;
            }
            newgenepair_avl->score1 = fabs(m1/n1 - m2/n2);
            newgenepair_avl->score2 = fabs(r1/n1 - r2/n2);
            
            
            /* Insert if applicable */
            if (tree->avl_count < max_nodes) {
                avl_insert(tree, newgenepair_avl);
                
                mingenepair_avl = avl_t_first(trav, tree);
            }
            else if (compare_genepair_avl(newgenepair_avl,mingenepair_avl,NULL) > 0) {
                avl_insert(tree, newgenepair_avl);
                
                avl_delete(tree, mingenepair_avl);
                free(mingenepair_avl);
                
                mingenepair_avl = avl_t_first(trav, tree);
            }
            else {
                free(newgenepair_avl);
            }
            
        }
    }
    
    free(trav);
    return tree;
    
}


SEXP mytspC2R(SEXP Rclass1, SEXP Rclass2, SEXP Rmax_nodes)
{
    int p,n1,n2,max_nodes,idx=0;
    struct avl_table *ptr_tree;
    struct genepair_avl *ptr_data;
    struct avl_traverser *ptr_trav;
    
    Rclass1 = coerceVector(Rclass1,REALSXP);
    Rclass2 = coerceVector(Rclass2,REALSXP);
    Rmax_nodes = coerceVector(Rmax_nodes,INTSXP);
    
    SEXP Rdim1,Rdim2,Rresult;
    
    Rdim1 = getAttrib(Rclass1, R_DimSymbol);
    Rdim2 = getAttrib(Rclass2, R_DimSymbol);
    p = INTEGER(Rdim1)[1];
    n1 = INTEGER(Rdim1)[0];
    n2 = INTEGER(Rdim2)[0];
    max_nodes = INTEGER(Rmax_nodes)[0];
    
    /* Create avl tree */
    
    ptr_tree = mytspC(REAL(Rclass1),REAL(Rclass2),n1,n2,p,max_nodes);
    
    /* Traversal by decreasing order */
    
    ptr_trav = (struct avl_traverser *) malloc (sizeof(struct avl_traverser));
    avl_t_init (ptr_trav, ptr_tree);
    
    PROTECT(Rresult=allocMatrix(INTSXP,max_nodes,2));
    
    for(ptr_data = avl_t_last(ptr_trav,ptr_tree), idx=0; ptr_data != NULL && idx < max_nodes; ptr_data = avl_t_prev(ptr_trav), idx++)
    {
        INTEGER(Rresult)[idx] = ptr_data->gene1;
        INTEGER(Rresult)[idx + max_nodes] = ptr_data->gene2;
    }
    
    free(ptr_trav);
    avl_destroy (ptr_tree, NULL);
    
    UNPROTECT(1);
    
    return(Rresult);
}

