/*#########################################################################
Karsten Krug, 20090403

functions used by MaxQuant Parser

###########################################################################*/

#include <windows.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <String.h>

#include "maxquantparser.h"

/*#########################################################################

  dependent peptide mass histogramm


##########################################################################*/


void dpMassHist(char **file, double *DPMassDiff, char **DPMod)
{
    FILE *fp;

    /*  pointers to store the columns of 'allPeptides.txt' */
    int *DPMaddPTR, count=0;
    char *DPModPTR;

    /* array to store a single line  */
    char line[1024];

    /* allocate memory for the first element  */
    DPMaddPTR = malloc(1*sizeof(int));
    DPModPTR = malloc(1*sizeof(char));

    /* opne the file */
    fp = fopen( **file, "r"  );


    /* read the file, line by line  */
    while( fgets( line, sizeof line, fp ) != NULL)
    {
	Rprintf("%s\n", line);


	    count++;
    }

    fclose(fp);
}


/*#########################################################################


                        GO term
                        y      n
                    ---------------
     cluster      y |  n11    n12  | nPC

                  n |  n21    n22  |
                    --------------- ---
                       nPGo        | N


		       over representation analysis


 **goIDs                 - pointer to a vector of character representing all annotation terms that were
                           assigned to proteins
 **protWithGoIDs         - pointer to a vector of characters that represent the mapping term <-> protein;
                           if a protein was not annotated by a term the respective array element is the empty string;
                           the length of this vector corresponds to the number of all proteins
 **clustering            - pointer to a vector of the same length as 'protWithGoIDs' illustrating the mapping protein <-> cluster;
                           the order of the elements has to be the same as in 'protWithGoIDs'!
 *nGo                    - integer, number of annotation terms
 *nProt                  - integer, number of all proteins, i.e. also proteins that were not annotated by a term
 *K                      - integer, number of clusters
 **test                  - string representing the test to use, at the moment: ''


########################################################################*/

void ora(char **goIDs, char **protWithGoIDs, char **clustering, int *nGo, int *nProt, int *K, char **test, double result[*nGo][*K], double *result2)
{
    int i=0, j=0, k=0, n11count=0, nPGo=0, nPC=0, n11=0, N=0,kk=0;
    double pVal=0;

    /* determine the number of proteins that were annotated with at least one GO term */
    for(i=0;i<*nProt;i++){
	if( strcmp( protWithGoIDs[i], "" ) != 0  ){
	    N++;
	}
    }

    /* loop over all GO terms */
    for(i=0; i < *nGo; i++){

	/* determine the number of proteins groups that are annotated with the current term */
	nPGo =  mygrep(goIDs[i], protWithGoIDs, *nProt);

	/* loop over the clusters */
	for(k=0; k<*K; k++)
	{
	    n11=0;

	    /* determine the number of proteins in the actual cluster (n1.)
	       as well as n11  */
	    nPC=0;
	    for( j=0; j < *nProt; j++ ){

		/* if the protein belongs to the current cluster
		   and is annotated with a GO term... */
		if( atoi(clustering[j]) == (k+1) &&  strcmp( protWithGoIDs[j], "" ) != 0   ){

		    /* increment number of proteins */
		    nPC++;

		    /* check if the protein is annotated with the current term */
		    if( strstr( protWithGoIDs[j], goIDs[i]  ) != NULL  ){
			n11++;
		    }
		}
	    }

	    /*#########################################################
	       estimate the p-value
	    ##########################################################*/
	    // Rprintf("%s:cluster %d: n11=%d, nPGo=%d, nPC=%d, N=%d\n", goIDs[i], k+1, n11, nPGo, nPC, N);
	    // Rprintf("\ti=%d, k=%d\n", i,k);

	    pVal=0.0;
	    /* exact test */
	    if( strcmp(*test, "exact")==0  )
	    {
		for( n11count=n11; n11count <= nPC; n11count++  )
		{
		    pVal += dhyper((double)n11count, (double)nPGo, (double)(N-nPGo), (double)nPC, 0);

		    //    Rprintf("\tn11count=%d, nPGo=%d, N-nPGo=%d , nPC=%d,dhyper=%f ", n11count, nPGo, N-nPGo, nPC,dhyper((double)n11count, (double)nPGo, (double)(N-nPGo), (double)nPC, 0));
		}
	    }

	    //Rprintf("\n");

	    result[i][k] = pVal;
	    result2[kk] = pVal;
	    kk++;


	    //Rprintf("\tresult[%d][%d]=%f\n", i, k, result[i][k]);

	    //result[i][j] = (double)i + (double)j;
	    //Rprintf("%f\n",result[i][j]);


	}

    }

}
/*####################################################################
  mygrep - the function counts the number of occurences of 'pattern'
           in the vector 'stringVec'

  arguments
	   *pattern     - character, the pattern that is searched
	   **stringVec  - character vector in which is searched

  value
           count  - the number of occurences within the vector
####################################################################*/
int mygrep(char *pattern, char **stringVec, int length)
{
    int i=0, count=0;

    for(i = 0; i < length; i++)
    {
        if(strstr(stringVec[i], pattern) != NULL){
	   count++;
        }
    }
    return count;
}

/*###################################################################
  mywhich
     - given a pattern and a vector the function returns the indices
       where 'pattern' was found in 'vector'

  arguments
  *pattern      - character, the pattern that is searched
  **stringVec   - character vector which is searched
  length        - length of 'stringVec'


######################################################################*/
int *mywhich( char *pattern, char **stringVec, int length)
{
    int i=0, tmp[length], count=0, *pRes;

    /* compare each element of the vector with 'pattern' */
    for(i=0; i<length; i++){

	/* if 'pattern' equals the actual element ...  */
	if( strcmp(pattern, stringVec[i]) == 0 ){
	    tmp[count] = i;
	    count++;

	}
    }

    /* now produce a vector containing only the positions where
       'pattern' was found */
    int result[count];
    for(i=0; i<count;i++){
	result[i] = tmp[i];
    }

    pRes=result;

    return pRes;
}



/*####################################################################

                                    peptidesInExperiments
    for each non-redundant peptide sequence get the number of experiments in which the sequence
    was identified

 arguments:
     **nrSeq        - pointer to a vector contianing all non-redundant peptide
                      sequences
     *nrSeqL        - integer, length of the vector
     **Seq          - pointer to a vector containing all peptide sequences, i.e.
                      the 'Sequence'-column in 'evidence.txt'
     *SeqL          - integer, length of **Seq and **Exp
     **Exp          - pointer to a vector containing the experiments for each
                      peptide sequence, i.e. the 'Experiment'-column in 'evidence.txt'
		    - **Seq and **Exp have to be in the same order!!
     **nrExp        - pointer to a vector containing the experiments non-redundantly
     **nrExpL       - integer, number of non-redundant experiments

#####################################################################*/
void peptidesInExperiments(char **nrSeq, int *nrSeqL, char **Seq, int *SeqL, char **Exp, char **nrExp, int *nrExpL, int *result )
{
    int i=0, j=0, k=0, l=0,countExpPerSeq, nrExpIdx[*nrExpL];

    /* loop over non-redundant peptide sequences */
    for(i=0; i<*nrSeqL; i++)
    {
	/* set all counters to zero */
	countExpPerSeq=0;
	for(l=0; l<*nrExpL; l++){ nrExpIdx[l] = 0;  }


	/* loop over redundant peptide sequences  */
	for(j=0; j < *SeqL; j++)
	{
	    /* identify all redundant peptide sequences */
	    if( strcmp(nrSeq[i], Seq[j] ) == 0){

		/* loop over non-redundnant experiments */
		for(k=0; k < *nrExpL; k++ ){

		    /* if the experiment was not counted before, i.e. 'nrExpIdx = 0', ...  */
		    if( (strcmp( Exp[j], nrExp[k]  )== 0 ) && (nrExpIdx[k]==0) )
		    {
			/* increment the number of experiments in which the actual peptide was identified  */
			countExpPerSeq++;
			/* set the index for the experiment to one in order to flag it as already counted*/
			nrExpIdx[k] = 1;
		    }

		}

	    }
	}

	/* store the number of experiements in which the current non- redundant
	   peptide sequence was identified                                     */
	result[i] = countExpPerSeq;
    }
}
/*##########################################################################################
                                         helper function

 split all multiple protein group IDs matching the same evidence to several rows, e.g.

 Protein.Group.IDs    Experiment
 997;28899              12well

                   to
 997                    12well
 28899                  12well


 this step is necessary to determine the correct number of identified protein groups
 in each experiment.

arguments
  **IDs        - column 'Protein.Group.IDs' of 'evidence.txt'
  *IDsL        - lenght of this vector
  **Exp        - columns 'Experiment' of 'evidence.txt'
  **resultID   -
  **resultExp  -

###########################################################################################*/
void splitProteinIDsInEvidence( char** IDs, int *IDsL, char **Exp,  char **resultID, char**resultExp )
{
    int i=0, j=0;
    char delim[] = ";";
    char *strsplit;


    for(i=0; i<*IDsL; i++){

	/* check if the actual ID consists of multiple IDs*/
	strsplit = strtok( IDs[i], delim );
	if(strsplit != NULL){

	    while(strsplit != NULL){
		resultID[j]  = strsplit;
		resultExp[j] = Exp[i];

		strsplit = strtok(NULL, delim);
		j++;
	    }
	}
	else{

	    resultID[j]=IDs[i];
	    resultExp[j]=Exp[i];

	   j++;
	}


    }
}
/*###########################################################################
                      helper function

 determine the number of all protein group ids in evidence.txt, i.e also consider
 multiple protein group IDs in a single row of evidence.txt

arguments
   IDs     - column 'Protein.Group.IDs' of 'evidence.txt'
   IDsL    - length of this vector
   IDnumb  - result
############################################################################*/
void NumberProteinIDsEvidence( char **IDs, int *IDsL, int *IDnumb)
{
    int i=0, j=0;
    char delim[]=";", *strsplit;

    /* loop over the IDs */
    for(i=0; i<*IDsL; i++)
    {
	/* check if there are multiple IDs  */
	strsplit = strtok(IDs[i], delim);

	/* if this is the case split the IDs */
	if(strsplit != NULL){
	    while(strsplit != NULL){
		/* increment the number */
		j++;

                /* set the pointer to the next ID */
		strsplit = strtok(NULL, delim);
	    }
	}
	else{
	    j++;
	}

    }

    *IDnumb = j;

}




/*

 */
void test4( int *n, int *m, double M[*n][*m] ){

    int i=0, j=0;

    for(i=0; i<*n;i++){
	for(j=0; j<*m; j++){
	    M[i][j] = i+j;
	}
    }

}

void test5( int *n, int *m, double *result ){

    int i=0, j=0, c=0;

    Rprintf("n=%d, m=%d\n", *n, *m );

    for(i=0; i<*n;i++){
	for(j=0; j<*m; j++){
	    result[c] = i+j;

	    c++;
	}

    }


}

/*######################################################################*/
void test(double *x, double *n, double *m, double *k, double *result)
{

    *result = dhyper(*x, *n, *m, *k, 0);

}


SEXP test2(SEXP x, SEXP n, SEXP m, SEXP k)
{
    double xx, nn, mm, kk,*presult;

    SEXP result;

    xx = REAL(coerceVector(x, REALSXP))[0];
    nn = REAL(coerceVector(n, REALSXP))[0];
    mm = REAL(coerceVector(m, REALSXP))[0];
    kk = REAL(coerceVector(k, REALSXP))[0];

    PROTECT( result = NEW_NUMERIC(1));
    UNPROTECT(1);

    presult=REAL(result);

    //*presult =  xx * nn;


    *presult = dhyper(xx, nn, mm, kk, 0);


    return result;
}
/*###################################################################



####################################################################*/
SEXP test3(SEXP goIDs, SEXP protWithGoIDs, SEXP clustering)
{
    //char *pvec;
    int nGo, nProt, nClust, i, *j;

    /* convert the argument vectors*/
    PROTECT( goIDs = coerceVector(goIDs, STRSXP) );
    PROTECT( protWithGoIDs = coerceVector(protWithGoIDs, STRSXP));
    PROTECT( clustering = coerceVector(clustering, INTSXP));
    UNPROTECT(3);

    /* number of GO terms */
    nGo = LENGTH(goIDs);
    printf("number of GO terms: %d\n", nGo);

    /* number of protein groups*/
    nProt = LENGTH(protWithGoIDs);

    // char z ="GO";

    //if( strstr( CHAR(VECTOR_ELT(goIDs,1)) , z  )) printf("jepp\n");

    /* loop over the terms */
    for(i = 0; i < nGo; i++){

	printf(" %s\n", CHAR( VECTOR_ELT(goIDs,i)) );

	/* determine the number of protein groups that are annotated with the current GO term*/


    }

    return(R_NilValue);

}
