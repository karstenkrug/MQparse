/*#########################################################################
Karsten Krug, 20090403

functions used by MaxQuant Parser


changelog:  20120313 cleaned up the code, only function required for the pdf script remained
###########################################################################*/

#include <windows.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <String.h>

#include "maxquantparser.h"


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

