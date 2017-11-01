/*
	Balestasis - learning Bayesian networks with heuristics (Please check README)
 
	Copyright 2017 罗良逸 ( Luo Liangyi, ラ リョウイツ)
 
	This file is the main part of Balestasis.

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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "externstuff.h"
#include "utility.c"
#include "compute.c"
#include "scoring.c"


/* the number of variables/nodes, extern class*/
unsigned int Xn;

/* the number of instances/cases, extern class */
unsigned long N; 

/* the natural log of N, extern class.*/
double lnN;

/* the size of ordering (chunk), extern class.*/
unsigned int On;

/* the limits for running time (not strict), extern class*/
unsigned int timeA, timeB; 

/* the pointer to the raw data of the input, extern class */
unsigned char *data; 

/* the pointer to orderings, extern class */
unsigned int *orderings;

/* the pointer to the Nodes (struct), extern class */
node *Nodes;

/* learnt Bayesian network structure for the output, extern class */
parentSet **finalNet;

/* the score of the above, extern class*/
double finalScore;

/* the pointer to ordering that leads to finalNet, which is needed to recover the network, extern class */
unsigned int *resultorder;

/* for switching the scoring metric, extern class */
char metric;

/* for switching the source and the mode of ordering, extern class */
char orderMode; 

/* for switching the source of caches or deciding whether to save the caches, extern class */
char cacheMode;

/* turns pruned search on, no point turning it off, extern class*/
char pruneMode;

/* determines whether to save the network to files, extern class*/
int netToFile;

/* output path and filename, extern class*/
char *outputPath;

/* keep a log of the running info, extern class*/
struct RuningLog zlog; 

/* auxiliary extern pointers, extern class */
double *aux0 = NULL, *aux1 = NULL, *aux2 = NULL;

int main(int argc, char *argv[]){
	
	printf("Started. ");

	time_t totalstart, totalfinish, checkpoint0, checkpointA0, checkpointA1, checkpointB0, checkpointB1, checkpointB2; //cpA for IS's actual completion, B for ASOBS's actually completion
	
	totalstart = time(NULL);
	
	unsigned int soRestart = 1;
	
	srand(time(NULL));
	
	zlog.psMinMaxSize = 255;
	zlog.parsetMaxSize = 1;
	strcpy(zlog.openlistExhausted, "YES");
	
	unsigned int argTimeA, argTimeB;

	cacheMode = argv[1][0]; //'d' raw data in, net out; 'h' raw data in, cache out; 'c' cache in, net out; 
	sscanf(argv[2], "%u", &argTimeA); // IS per-node time limit
	sscanf(argv[3], "%u", &argTimeB); // ASOBS time limit
	sscanf(argv[4], "%s", zlog.logFile); //64-char-limit
	sscanf(argv[5], "%s", zlog.inputFile); //same
	sscanf(argv[11], "%s", zlog.orderFile); //same
	sscanf(argv[10], "%s", zlog.outputFile); //same
	orderMode = argv[6][0]; //'u' for random uniform, 'l' for sample learning, 'm','p' for preloaded orders from file
	sscanf(argv[7], "%u", &soRestart); //if 'u' or 'l', then this can be set
	//argv[8] reserved
	sscanf(argv[9], "%d", &netToFile); // 0 not saved or 1 saved

	if(netToFile == 1) {outputPath = zlog.outputFile;}
	
	
	timeA = argTimeA; //This will be the ESTIMATED time limit of learning parent sets for A node.

	timeB = argTimeB; //This will be the ESTIMATED time limit of building structure using ALL nodes.


	metric = 'B'; //'B' for BIC. Currently, only BIC makes sense.
	zlog.metricUsed[0] = metric;
 
	pruneMode = 1; //It seems that 1 here is always better.

	/* loading data below */
 	FILE *fp;
	fp = fopen(zlog.inputFile, "r");

	fscanf(fp, "%u", &Xn);
	zlog.nodesNum = Xn;

/*!!!*/node Variables[Xn];
	
	if(cacheMode != 'c') fscanf(fp, "%lu", &N);

	unsigned int i;
	for(i=0; i<Xn; i++){

		fscanf(fp, "%s", Variables[i].name); //Use less than 20 chars for names !!! Be cautious of overflow! 

	}

	for(i=0; i<Xn; i++){

		fscanf(fp, "%hhu", &Variables[i].statesCount); 
		Variables[i].datumAC = i;

	}
	
/*!!!*/Nodes = Variables;
	
	parentSet **cacheofscores[Xn]; //caches declared here

	unsigned long cachesizes[Xn]; //the sizes of the caches
	
	printf("The cache mode is %c.\n", cacheMode);
	
	//If use existing caches loaded through the input,
	if(cacheMode == 'c'){
		
		unsigned int z;
		for(z=0; z<Xn; z++){
			
			unsigned int thatNode;
			
			fscanf(fp, "%u %lu", &thatNode, &cachesizes[z]);
			
			parentSet **closedlist;

			closedlist = malloc( cachesizes[z] * sizeof(parentSet *)); //Don't free!!!
			
			unsigned int i;
			for(i=0; i<cachesizes[z]; i++){
				
				double score;
				
				fscanf(fp, "%lf", &score);
				
				unsigned int parSize = 0;
				
				fscanf(fp, "%u", &parSize);
				
				unsigned int temp[parSize];
				
				unsigned int a;
				for(a=0; a<parSize; a++) fscanf(fp, "%u", &temp[a]);

				closedlist[i] = createParentSet(parSize*sizeof(unsigned int), temp);
				closedlist[i]->score = score;
				
			}
			
			cacheofscores[z] = closedlist;

		}
		fclose(fp);
		
	}
	
	printf("System initialisation.\n");
	if(orderMode == 'p'){
		
		printf("Using preloaded orderings.\n");
		
		soRestart = 1; //ASOBS should run only once when orderings are preloaded
	
		fp = fopen(zlog.orderFile, "r");

		unsigned int posts;

		fscanf(fp, "%u", &posts);

		if(posts != Xn){

			printf("Ordering data mismatch.\n");
			exit(0);

		}


		fscanf(fp, "%u", &On);

		printf("There are %u preloaded orderings.\n", On);

		orderings = malloc(Xn*On*sizeof(unsigned int));
		
		unsigned int iv;
		for(iv=0; iv<Xn*On; iv++){

			fscanf(fp, "%u", &orderings[iv]);

		}
		fclose(fp);
		
		printf("Loading of orderings accomplished.\n");
		
	}
	else if(orderMode == 'u' || orderMode == 'l' ){
		
		On = omp_get_max_threads(); // One ordering per thread at a time.
		orderings = malloc(Xn*On*sizeof(unsigned int));
		
	}
	else{
		
		On = 1;
		orderings = malloc(Xn*sizeof(unsigned int));
		
	}
	
/*Portal Blue*/ if(cacheMode == 'c') goto Portal;	
	
	printf("Loading raw data.\n");
		
	lnN = log(N);
	zlog.dataNum = N;
	
	printf("There are %u variables/nodes and %lu instances/data points.\n", Xn, N);

	data = malloc(Xn*N); //unsigned char: assuming no variable has more than 256 states.

	for(i=0; i<Xn*N; i++){

		fscanf(fp, "%hhu", &data[i]);

	}

	fclose(fp);
	printf("Loading of the data set accomplished.\n");
	
	printf("Loading accomplished.\n");
	
	totalfinish = time(NULL);
	
	checkpoint0 = totalfinish - totalstart;
	
	printf("Time passed = %lds\n", checkpoint0);

/* This marks the end of loading data. Computation ensues. */


/* Next: learning parent sets from raw data (caches construction) */
	printf("\nCommencing cache construction.\n");
	checkpointA0 = time(NULL);

	unsigned int z, completed = 0;
	#pragma omp parallel for
	for(z=0; z<Xn; z++){

		cacheofscores[z] = indieSelection(z, &cachesizes[z]);
		completed += 1; //no reduction because it's non-critical info for calculating and displaying progress
		printf("\r  Cache construction for %u out of %u nodes complete (Progress: %u%%).", completed, Xn, (unsigned int) ((float) completed/(float) Xn*100));
		fflush(stdout);

	}
	checkpointA1 = time(NULL);
	
	zlog.psiTime =  checkpointA1 - checkpointA0;
	
	printf("\n Cache construction complete! Time usage = %lds\n", zlog.psiTime );
	
	totalfinish = time(NULL);

	printf("Time passed = %lds\n", totalfinish - totalstart);

	if(netToFile == 0) { free(data); data = NULL; }
	
/* check against anomalies */
	unsigned int r;
	for(r=0; r<Xn; r++){for(z=0; z<cachesizes[r]; z++) parentSetDisplay(cacheofscores[r][z], 0, r, z);}
	
/* Portal Orange */Portal: printf("Caches are ready.\n"); 


/* Save the caches to a file then stop.*/

	if(cacheMode == 'h'){
		
		cachesToFile(cacheofscores, cachesizes);
		
		exit(0);
		
	}


	
/* Next: building global structure */

	unsigned int exprun = 0;
	
	finalScore = -DBL_MAX; //It is here meaning only the best network of all runs will be saved
	
	resultorder = malloc(Xn*sizeof(unsigned int));

/*Portal Red*/Portal2: 

	zlog.iCter = 0;

	exprun += 1;
	
	zlog.iCter = 0; //reset
	
	if(orderMode != 'p' ) for(i=0; i<Xn*On; i++) orderings[i] = i%Xn; //reset; the starting ordering is always 0,1,3....n if not preloaded
	
	printf("\nCommencing global structure construction.\n");
	checkpointB0 = time(NULL);
	
	//p for preloaded, needs a preloaded order file
	if(orderMode == 'p') {
		
		completed = 0;

		unsigned int (*ordering)[Xn] = (unsigned int (*)[Xn])&orderings[0]; 
		
		printf("ASOBS for %u preloaded orderings." ,On);
		
		unsigned int w;
		#pragma omp parallel for
		for(w=0; w<On; w++){
			
			parentSet **anet = asOBSlite(cacheofscores, cachesizes, ordering[w]);
			completed += 1; //for calculating and displaying progress only, not 100% accurate
				
			double anetscore = 0.0;
	
			unsigned int f;
			for(f=0; f<Xn; f++) anetscore = anetscore + anet[f]->score;
				
			if(anetscore > finalScore){
				
				#pragma omp critical
				{
					finalNet = anet; 
					
					finalScore = anetscore;
					
					unsigned int u;
					for(u=0; u<Xn; u++) resultorder[u] = ordering[w][u];

					zlog.trialsTBT = w; //In cases of preloaded orderings, the index of the best ordering will be logged.
				}
			}
			else free(anet); 
			//Warning!!! Don't free any of the anet[i]. Otherwise the caches will be destoryed.

			printf("\r %u out of %u complete (Progress: %u%%).", completed, On, (unsigned int) ((float) completed/(float) On*100));
			fflush(stdout);

		}
		checkpointB1 = time(NULL);
		
		zlog.totaltrials = On;
		zlog.soTime = checkpointB1 - checkpointB0;
		
	}
	else{
		
		double topScore = -DBL_MAX;

		unsigned long trials = 0,  scorelistCounter = 0, scorelistSize = 2048; //scorelistCounter and scorelistSize are for 'l' mode
		
		unsigned int (*ordering)[Xn] = (unsigned int (*)[Xn])&orderings[0]; 		
		
		unsigned int th;
		#pragma omp parallel for
		for(th=0; th<On; th++){
			
			orderSearcher(&topScore, &trials, ordering[th], cacheofscores, cachesizes, &scorelistCounter, &scorelistSize);
			
		}
		zlog.totaltrials = trials;
		zlog.score = topScore;
		free(aux0); aux0 = NULL;  free(aux1); aux1 = NULL; free(aux2); aux2 = NULL;
		
		printf("\n  Best network score in this run: %lf\n", topScore);
		
	}
	
	checkpointB2 = time(NULL);

	zlog.soTime = checkpointB2 - checkpointB0;
	
	printf("\n  Time usage = %lds\n", zlog.soTime );
	
	zlog.togetherT = zlog.psiTime + zlog.soTime;
	
	totalfinish = time(NULL);
	printf("Time passed = %lds\n", totalfinish - totalstart); //This is the accumulated time since the start

/* Next: save the running log to the file */
	fp = fopen(zlog.logFile, "a+");
	
	fprintf(fp, "\nTest No. %u\n", exprun);
	
	if(cacheMode == 'd') fputs("Mode: raw data in, net out\n", fp);
	else if(cacheMode == 'h') fputs("Mode: raw data in, caches out\n", fp);
	else fputs("Mode: caches in, net(structure only) out\n",fp);
	
	fprintf(fp, "{\n\
	Input: %s ;\n\n\
	Nodes: %u ;\n\n\
	Instances: %lu ;\n\n\
	Preloaded Ordering: %s ;\n\n\
	",  zlog.inputFile, zlog.nodesNum, zlog.dataNum, zlog.orderFile);
			
	fprintf(fp, "Metric: %s ;\n\n\
	Score: %lf ;\n\n\
	Total running time: %ld ;\n\n\
	", zlog.metricUsed, zlog.score, zlog.togetherT);

	if(cacheMode != 'c') fprintf(fp, "IS time limit setting (per node): %u ;\n\n\
	Actual IS running time: %ld ;\n\n\
	Size of the largest parent set explored: %hhu - %hhu ;\n\n\
	IS open lists all exhausted: %s ;\n\n\
	", timeA, zlog.psiTime, zlog.psMinMaxSize, zlog.parsetMaxSize, zlog.openlistExhausted);

	fprintf(fp, "ASOBS time limit setting: %u ;\n\n\
	Actual ASOBS running time: %ld ;\n\n\
	Ordering mode: %c ;\n\n\
	Orderings sampled: %lu ;\n\n\
	Improvement history: %u {\n\n\
	",   timeB, zlog.soTime, orderMode, zlog.totaltrials, zlog.iCter);
	
	for(z=0; z<zlog.iCter; z++) fprintf(fp, " %ld ", zlog.improveT[z]); 
	fputs(";\n\n", fp);


	for(z=0; z<zlog.iCter; z++) fprintf(fp, " %lu %lf %u\n", zlog.tCheckP[z], zlog.scoreH[z], exprun); 
	fputs(";\n	", fp);
		
	
	fputs("		}\n\n}\n\n", fp);
	
	fclose(fp);
	
	zlog.overallT = time(NULL) - totalstart;
	
/*Portal Purple*/if(exprun < soRestart) goto  Portal2;  //The number of soRestart is the number of ASOBS reruns.
	
/* Coda */
	double displayScore = 0.0;

	for(z=0; z<Xn; z++) {displayScore = displayScore + finalNet[z]->score;} //A final check.
	
	printf("\nAll done! Final (best) network score: %lf\n", displayScore);
	
	zlog.score = finalScore;
	
	//for(z=0; z<Xn; z++) {printf("%u ", resultorder[z]);}
		
	if(netToFile == 1 && cacheMode != 'c'){ netToDSCFile(); printf("Network saved to a file: %s\n", outputPath); }
	
	if(netToFile == 1 && cacheMode == 'c'){ 
		
		structureToFile();
		
	}
	
	fp = fopen(zlog.logFile, "a+");
	
	fprintf(fp, "Total time: %ld\nFinal score: %lf\nOutput: %s", zlog.overallT, zlog.score, outputPath);
	
	fclose(fp);
	
	printf("\nExit.\n");
	
	return EXIT_SUCCESS;
	
}

