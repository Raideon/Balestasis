/*
	Copyright 2017 罗良逸 ( Luo Liangyi, ルオ．りょういつ)
 
	This file is part of Balestasis

    Balestasis is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Balestasis is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Balestasis.  If not, see <http://www.gnu.org/licenses/>.
*/

/*-----------variables-----------*/
/* the number of variables (last node is the Xnth node, but it's codename is Xn-1) */
extern unsigned int Xn;

/* the number of cases/instances/data points in the data set*/
extern unsigned long N; 

/* the natural log of N */
extern double lnN;

/* the number of pre-loaded orderings */
extern unsigned int On;

/* used as the limits for anytime algorithms' running time*/
extern unsigned int timeA, timeB; 

/* for scoring metrics switching */
extern char metric;

/* for switching the source of ordering */
extern char orderMode;

/* determines pruned search */
extern char pruneMode;

/* determines whether save the network to file (with parameters by MLE) */
extern int netToFile;

/* the score of the above, extern class*/
extern double finalScore;

/*-----------pointers------------*/
/* output path and filename*/
extern char *outputPath;

/* pointer to data */
extern unsigned char *data; //assuming no varible has more than 256 discrete states.

/* pointer to orderings */
extern unsigned int *orderings;

/* pointer to the result ordering, which is needed to recover the network */
extern unsigned int *resultorder;

/* auxiliary extern pointers */
extern double *aux0, *aux1, *aux2;


/*-----------struct--------------*/
/* node information */
struct Node {

	char name[20];
	unsigned char statesCount;
	unsigned int datumAC; // a datum access coded is used for accessing a datum: same as the column number, range: 0 to Xn-1, it is also know as the codename of a node

};
typedef struct Node node;

/* the pointer to the Nodes struct array. extern class */
extern node *Nodes;


/* parent set */
struct ParentSet {

	double score;
	double entropy;
	unsigned int size;
	unsigned int configsCount;
	unsigned int parents[1];

};
typedef struct ParentSet parentSet;

/* the Bayesian network structure for output */
extern parentSet **finalNet;

/* BIC* */
struct BICastk {

	double score;
	parentSet *pis[2];
	unsigned int parlist[1];

};
typedef struct BICastk bicAstk;

/* running log */
struct RuningLog {

	char logFile[64];
	char inputFile[64];
	char orderFile[64];
	char outputFile[64];
	char metricUsed[5];
	char openlistExhausted[3];
	unsigned int nodesNum;
	unsigned long dataNum;
	time_t psiTime;
	time_t soTime;
	time_t soTBT; //not used
	time_t togetherT;
	time_t overallT;
	time_t improveT[256];
	unsigned char psMinMaxSize;
	unsigned char parsetMaxSize; //explored
	unsigned int iCter; //the counter of improvements
	unsigned long trialsTBT; //not used
	unsigned long totaltrials;
	unsigned long tCheckP[256];
	double score;
	double scoreH[256];

};

/* keep a log of key info */
extern struct RuningLog zlog; 
