/*
	Balestasis - software for learning Bayesian networks from data with heuristics (Please check README)
 
	Copyright 2017 罗良逸 ( Luo Liangyi, ルオ．りょういつ)
 
	This file is the main part of Balestasis

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

/* limits for running time (not strict), extern class*/
unsigned int timeA, timeB; 

/* the pointer to data, extern class */
unsigned char *data; 

/* pointer to testing orderings, extern class */
unsigned int *orderings;

/* the pointer to the Nodes (struct), extern class */
node *Nodes;

/* the Bayesian network structure for output, extern class */
parentSet **finalNet;

/* the score of the above, extern class*/
double finalScore;

/* pointer to the result ordering, which is needed to recover the network, extern class */
unsigned int *resultorder;

/* for switching the scoring metric, extern class */
char metric;

/* for switching the source and mode of ordering, extern class */
char orderMode; 

/* turns pruned search on, no point turning it off, extern class*/
char pruneMode;

/* determines whether save the network to file, extern class*/
int netToFile;

/* output path and filename, extern class*/
char *outputPath;

/* keep a log of running info, extern class*/
struct RuningLog zlog; 

/* auxiliary extern pointers, extern class */
double *aux0 = NULL, *aux1 = NULL, *aux2 = NULL;

int main(int argc, char *argv[]){

	time_t totalstart, totalfinish, checkpoint0, checkpointA0, checkpointA1, checkpointB0, checkpointB1, checkpointB2; //cpA for IS's actual completion, B for ASOBS's actually completion
	
	totalstart = time(NULL);
	
	unsigned int soRestart = 1;
	
	unsigned char testingDifferentOrderingStrategyMode = 0; //This (experiment mode switch) decides whether orderMode will change in the middle of multiple runs. 
	
	srand(time(NULL));
	
	zlog.psMinMaxSize = 255;
	zlog.parsetMaxSize = 1;
	strcpy(zlog.openlistExhausted, "YES");
	
	unsigned int argTimeA, argTimeB;

	//argv[1] reserved for metric switching, currently only BIC/MDL score is supported
	sscanf(argv[2], "%u", &argTimeA); // IS per-node time limit
	sscanf(argv[3], "%u", &argTimeB); // ASOBS time limit
	sscanf(argv[4], "%s", zlog.logFile); //64-char-limit
	sscanf(argv[5], "%s", zlog.inputFile); //same
	orderMode = argv[6][0]; //'u' for random uniform, 'l' for sample learning, 'p' for preloaded orders from file
	sscanf(argv[7], "%u", &soRestart); //if 'u' or 'l', then this can be set
	sscanf(argv[8], "%s", zlog.orderFile); //same
	sscanf(argv[9], "%d", &netToFile); // 0 or 1
	sscanf(argv[10], "%s", zlog.outputFile); //same

	if(netToFile == 1) {outputPath = zlog.outputFile;}
	
	
	timeA = argTimeA; //This will be the ESTIMATED time limit of parent set searching for A node.

	timeB = argTimeB; //This will be the ESTIMATED time limit of structure optimisation ALL nodes.


	metric = 'B'; //'B' for BIC. Currently, only using BIC makes sense.
	zlog.metricUsed[0] = metric;
 
	pruneMode = 1; //I put a switch here in the early stage. But it turns out it is always better to set it as ON(1).

	/* loading data below */
 	FILE *fp;
	fp = fopen(zlog.inputFile, "r");

	fscanf(fp, "%u", &Xn);

	fscanf(fp, "%lu", &N);
	
	lnN = log(N);

	printf("There are %u variables/nodes and %lu instances/data points.\n", Xn, N);
	zlog.nodesNum = Xn;
	zlog.dataNum = N;

/*!!!*/node Variables[Xn];

	unsigned int i;
	for(i=0; i<Xn; i++){

		fscanf(fp, "%s", Variables[i].name); //Be cautious of overflow! Less than 20 chars allowed.

	}

	for(i=0; i<Xn; i++){

		fscanf(fp, "%hhu", &Variables[i].statesCount); 
		Variables[i].datumAC = i;

	}

/*!!!*/Nodes = Variables;

	data = malloc(Xn*N); //unsigned char: assuming no variable has more than 256 states.

	for(i=0; i<Xn*N; i++){

		fscanf(fp, "%hhu", &data[i]);

	}
	fclose(fp);
	printf("Loading of the data set accomplished.\n");

	if(orderMode == 'p'){
		
		soRestart = 1; //make sure that ASOBS only run once when orderings are preloaded
	
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

		for(i=0; i<Xn*On; i++){

			fscanf(fp, "%u", &orderings[i]);

		}
		fclose(fp);
		
		printf("Loading of orderings accomplished.\n");
		
	}
	else if(orderMode == 'u' || orderMode == 'l' ){
		
		On = omp_get_max_threads();
		orderings = malloc(Xn*On*sizeof(unsigned int));
		
	}
	else{
		
		On = 1;
		orderings = malloc(Xn*sizeof(unsigned int));
		
	}
	
	printf("Loading accomplished.\n");
	
	totalfinish = time(NULL);
	
	checkpoint0 = totalfinish - totalstart;
	
	printf("Time passed = %lds\n", checkpoint0);

/* This marks the end of loading data. Computation ensues. */


/* Next: parent set identification */
	printf("\nCommencing cache construction.\n");
	checkpointA0 = time(NULL);

	parentSet **cacheofscores[Xn];

	unsigned long cachesizes[Xn];

	unsigned int z, completed = 0;
	#pragma omp parallel for
	for(z=0; z<Xn; z++){

		cacheofscores[z] = indieSelection(z, &cachesizes[z]);
		completed += 1; //no reduction because it's non-critical info for calculating progress
		printf("\r  Cache construction for %u out of %u nodes complete (Progress: %u%%).", completed, Xn, (unsigned int) ((float) completed/(float) Xn*100));
		fflush(stdout);

	}
	checkpointA1 = time(NULL);
	
	zlog.psiTime =  checkpointA1 - checkpointA0;
	
	printf("\n Cache construction complete! Time usage = %lds\n", zlog.psiTime );
	
	totalfinish = time(NULL);

	printf("Time passed = %lds\n", totalfinish - totalstart);

	if(netToFile == 0) { free(data); data = NULL; }

/*Check against anomaly*/
	unsigned int r; for(r=0; r<Xn; r++){for(z=0; z<cachesizes[r]; z++) parentSetDisplay(cacheofscores[r][z], 0, r, z);}


/* Next: structure optimisation */

	unsigned int exprun = 0;
	
	finalScore = -DBL_MAX; //It is here; therefore, only the best network of all runs will be saved
	
	resultorder = malloc(Xn*sizeof(unsigned int));

/*Yellow*/Portal: 

	zlog.iCter = 0;

	exprun += 1;
	
	if(testingDifferentOrderingStrategyMode == 1){
		
		if(exprun % 2 == 0) orderMode = 'l'; 
		else orderMode = 'u';
		
	}

	
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
			completed += 1; //for calculating progress only, not 100% accurate
				
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

					zlog.trialsTBT = w; //In cases of preloaded ordering, the index of the best ordering will be logged.
				}
			}
			else free(anet); 
			//Warning!!! To whomever is thinking I didn't free properly here. I did. Don't free any of the anet[i]. If you do, you will destory the caches.

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
	
	fprintf(fp, "Test No. %u\n", exprun);
	
	fprintf(fp, "{\n\
	Data: %s ;\n\n\
	Nodes: %u ;\n\n\
	Instances: %lu ;\n\n\
	Preloaded Ordering: %s ;\n\n\
	",  zlog.inputFile, zlog.nodesNum, zlog.dataNum, zlog.orderFile);
			
	fprintf(fp, "Metric: %s ;\n\n\
	Score: %lf ;\n\n\
	IS+ASOBS running time: %ld ;\n\n\
	", zlog.metricUsed, zlog.score, zlog.togetherT);

	fprintf(fp, "IS time limit setting (per node): %u ;\n\n\
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
	
/*Blue*/if(exprun < soRestart) goto  Portal;  //The number of soRestart is the number of reruns for the ASOBS
	
/* Coda */
	double displayScore = 0.0;

	for(z=0; z<Xn; z++) {displayScore = displayScore + finalNet[z]->score;} //This actually provides a final check
	
	printf("\nAll done! Final (best) network score: %lf\n", displayScore);
	
	zlog.score = finalScore;
	
	//for(z=0; z<Xn; z++) {printf("%u ", resultorder[z]);}
		
	if(netToFile == 1){ netToDSCFile(); printf("Network saved to a file: %s\n", outputPath); } 
	
	fp = fopen(zlog.logFile, "a+");
	
	fprintf(fp, "Total time: %ld\nFinal score: %lf\nOutput network: %s", zlog.overallT, zlog.score, outputPath);
	
	fclose(fp);
	
	printf("\nExit.\n");
	
	return EXIT_SUCCESS;
	
}

