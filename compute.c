/*
	Copyright 2017 罗良逸 ( Luo Liangyi, ラ リョウイツ)
 
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

void localscore(unsigned int Xi, parentSet *Pii, unsigned char mode);

parentSet* emptyParentSet(unsigned int Xi);

bicAstk* BICasterisk(unsigned int Xi, parentSet *PiL, parentSet *PiR, parentSet *empty);

parentSet* BICastCrtSet(unsigned int Xi, bicAstk *bicXiPiLPiR);

void pruner(parentSet **closedlist, unsigned long *cc);

void swap(unsigned int *ordering, unsigned int position);

int comparescores(const void *x, const void *y);

void anceTracer(parentSet **net, unsigned int *list, unsigned int *c, unsigned int *ordering, unsigned int i, unsigned int *ra);

void shuffler(unsigned int *ordering);

parentSet** asOBSlite(parentSet ***cache, unsigned long *cachesize, unsigned int *ordering);


void orderLearner(unsigned long *cct, unsigned long *sasize, double *netscore, unsigned int *ordering){

	if( *cct == 0 ){
						
		aux0 = malloc(*sasize*sizeof(double));
			
		unsigned long varvarsize = Xn*Xn;
					
		aux1 = malloc(varvarsize*sizeof(double));
			
		aux2 = malloc(varvarsize*sizeof(double));
	
		unsigned long a;
		for(a=0; a<*sasize; a++) aux0[a] = -DBL_MAX;
		for(a=0; a<varvarsize; a++){ aux1[a] = 1; aux2[a] = 1; }

	}
	else if( *cct >= *sasize-On){

		unsigned long newsize = *sasize+*sasize;
		aux0 = realloc(aux0, newsize*sizeof(double));
		
		unsigned long a;
		for(a=*sasize; a<newsize; a++) aux0[a] = -DBL_MAX;
		*sasize = newsize;			
	}

	double (*varsr)[Xn] = (double (*)[Xn])&aux1[0]; //varsr[0][1] == 3 means var 0 appeared 3 times as the 1st node in the ordering
	double (*varsro)[Xn] = (double (*)[Xn])&aux2[0]; //varsro[0][1] == 4000 means the sum of score ranking is 4000 for all cases where var 0 has been the 1st node so far
	
	unsigned long ranking; //ranking means the ranking of *netscore among all scores

	*cct += 1;
	
	unsigned long i;
	for(i=0; i<*cct; i++){
		
		if(*netscore > aux0[i]){
			
			if( i < *cct-1 ){
				unsigned long j;
				for(j=*cct-1; j>i; j--){
					
					aux0[j] = aux0[j-1];
					
				}
			}
			aux0[i] = *netscore;
			ranking = i;
			break;
		}
	}
	
	for(i=0; i<Xn; i++){
		
		varsr[ordering[i]][i] += 1;
		varsro[ordering[i]][i] += ranking;
		
	}
	
	unsigned int newRanks[Xn], // ranks in orderings.
				 candidates[Xn],// a list of nodes to be exhausted
				 poselOrder[Xn];// another order for rank(position) selection, aka the auxiliary order


	for(i=0; i<Xn; i++) candidates[i] = poselOrder[i] = i;
	

	shuffler(poselOrder); //pseudo-randomly arrange the order (uniform)

	//The program will draw a uniformly shuffled order or get an order based on patterns with 50/50 chance.
	int which = rand()%100;
	if( which < 50 ) shuffler(ordering); //can be changed
	else{

		for(i=0; i<Xn; i++){
			
			double top = 0.0, current;
			unsigned int candidate, candipos;
				
			unsigned int j;
			for(j=0; j<Xn; j++){
				
				if(candidates[j] < Xn){ 
					
					current = varsr[candidates[j]][poselOrder[i]] / varsro[candidates[j]][poselOrder[i]];

					if(current > top){ top = current; candidate = candidates[j]; candipos = j; }
				}
			}
			newRanks[poselOrder[i]] = candidate;

			candidates[candipos] = Xn; //cross the node who has gotten a rank/position in this round

		}
		
		for(i=0; i<Xn; i++) ordering[i] = newRanks[i]; 

	}
}

void orderSwapper(unsigned int *ordering, unsigned int *position){

	if(*position < Xn-1){

		swapper(ordering, position);
		
		*position += 1;
		
	}
	else{
		
		*position = 0;
		
	}

}

void orderSearcher(double *topScore, unsigned long *trials, unsigned int *ordering, parentSet ***cache, unsigned long *cachesize, unsigned long *scorelistCounter, unsigned long *scorelistSize){
	
	time_t started, finishing;
	started = time(NULL);
	finishing = time(NULL);	
	
	unsigned int position = 0;
	
	double previouScore;
	
	while((finishing - started) < timeB){
			
/*!!!*/ parentSet **anet = asOBSlite(cache, cachesize, ordering);
		finishing = time(NULL);	

		#pragma omp atomic
		*trials+=1;

		unsigned int event = *trials;

		double anetscore = 0.0;
	
		unsigned int f;
		for(f=0; f<Xn; f++) anetscore = anetscore + anet[f]->score;
		
		//only the all time best network will be saved
		#pragma omp critical
		{
			if(anetscore > finalScore && ((finishing - started) < timeB)){
			
			
				finalNet = anet;
				finalScore = anetscore;
				
				unsigned int q;
				for(q=0; q<Xn; q++) resultorder[q] = ordering[q];
			}
			else free(anet); 
		}
		
		//Warning!!! Don't free any of the anet[i]!!!
		
		//Only the all time best is saved, but all running bests' score and relevant info will be logged
		#pragma omp critical
		{
			if(anetscore > *topScore && ((finishing - started) < timeB)){
		

				*topScore = anetscore; 
				
				zlog.improveT[zlog.iCter] = finishing - started;
				zlog.scoreH[zlog.iCter] = anetscore;
				zlog.tCheckP[zlog.iCter] = event;
				zlog.iCter += 1;
				
			}
		}

		
		if( orderMode == 'l' )
			#pragma omp critical
			orderLearner(scorelistCounter, scorelistSize, &anetscore, ordering);
			
		else if( orderMode == 's' || orderMode == 'm' ){
			
			if(anetscore < previouScore && position != 0){ //reverse the swap if worse than before
				
				position-=1;
				
				orderSwapper(ordering, &position);

			}
			
			previouScore = anetscore;
			
			orderSwapper(ordering, &position);
			
			if(position == 0){
				
				if( orderMode == 's') shuffler(ordering);
				else if( orderMode == 'm') orderLearner(scorelistCounter, scorelistSize, &anetscore, ordering);

			} 
			
		}
		
		else shuffler(ordering);
			
		
		finishing = time(NULL);
		printf("\r  Remaining time: %lds; %lu orderings tried. Current score: %lf ", (long) (timeB-(finishing-started)), *trials, *topScore );
		fflush(stdout);
	}
}


/* a barebone asOBS function suitable for concurrent call
   literature: Scanagatta et al., Learning Bayesian Networks with Thousands of Variables, 2015 */
parentSet** asOBSlite(parentSet ***cache, unsigned long *cachesize, unsigned int *ordering){

	parentSet **net;
	net = malloc(Xn * sizeof(parentSet *));

	unsigned char table[Xn][Xn]; //table[A][B] == 1 means A is B's parent or ancestor
	memset(table, 0, sizeof table); 
	
	unsigned int i;
	for(i=0; i<Xn; i++){

		unsigned char accord = 0;

		unsigned int currentnode = ordering[i];

		unsigned int j = 0;
		while( accord != 1 ){

			unsigned int z;
				
			if(cache[currentnode][j]->size > 0){
				
				for(z=0; z<cache[currentnode][j]->size; z++){

					if(table[currentnode][cache[currentnode][j]->parents[z]] == 1){

						accord = 0; break;

					}
					accord = 1;			
				}
			}
			else{ accord = 1; }	
			j+=1;
		}
		j-=1;

		net[i] = cache[currentnode][j];

		if(net[i]->size > 0){
			
			unsigned int todo[Xn];
			todo[0] = currentnode;
			
			unsigned int t, g, u = 1;
			for(t=0; t<Xn; t++){

				if(table[currentnode][t] == 1){
					
					todo[u] = t; u+=1;

				}
			}

			unsigned int ancestors[Xn];
			unsigned int ancePop = 0;
			anceTracer(net, ancestors, &ancePop, ordering, i, &i);

			for(t=0; t<ancePop; t++){
					
				for(g=0; g<u; g++){

					table[ancestors[t]][todo[g]] = 1;
	
				}
			}
		}
	}
	
	return net;

}


/*  An function for searching parent sets based on independence selection algorithm. 
    literature: Scanagatta et al., Learning Bayesian Networks with Thousands of Variables, 2015 */
parentSet** indieSelection(unsigned int Xi, unsigned long *cachesize){

	time_t started, finishing;

	started = time(NULL);

	unsigned long c, cc;
	c = Xn+Xn;
	cc = 0; //This is the counter of how many elements in closedlist

	unsigned long o, oc;
	o = Xn*Xn;
	oc = 0; //Same, but for openlist

	parentSet **closedlist;

	closedlist = malloc( c * sizeof(parentSet *)); //This shouldn't be freed until the very end, well, just don't free it.

	bicAstk **openlist;
	openlist = malloc( o * sizeof(bicAstk *));

	/*Empty and parentSet of size 1 creating, scoring, and adding to closedlist
	  the first batch of size 1 parent set will match the codename 
	  order of nodes, that is closedlist[3] will point to parent set {3}.
	  closedlist[Xi] will point to empty parent set { }.
	  This will be the case until closedlist is sorted by score.
	  All for the ease of computation.
	*/
	unsigned int temp[1];
	unsigned int n, m;
	for(n=0; n<Xn; n++){

		if(n != Xi){

			temp[0] = n;
			closedlist[cc] = createParentSet(sizeof(temp), temp);
			localscore(Xi, closedlist[cc], metric);
			cc+=1;
		}
		else{
			closedlist[cc] = emptyParentSet(Xi);
			cc+=1;
		}
	}
	
	/* ParentSet of size 2 creating, *scoring, and add to openlist. */
	for(n=0; n<Xn; n++){

		if(n != Xi){

			for(m=n+1; m<Xn; m++){

				if(m != Xi){

					openlist[oc] = BICasterisk(Xi, closedlist[n], closedlist[m], closedlist[Xi]);
					oc+=1;

				}
			}
		}
	}

	finishing = time(NULL);

	if((finishing - started) >= timeA) printf("\n Warning! Allocated running time is insufficient for meaningful result.");

	unsigned int biggest = 2; //Keep a record of the size of the biggest parent set checked

	while( oc != 0 && (finishing - started) < timeA ){


		unsigned long highest, k;
		highest = 0; //reset

		double highestscore = openlist[0]->score; //reset

		unsigned char sizeincrease = 0; //reset

		//printf("Check point 0!\n");

		/* Find the parent set with highest BIC* */
		for(k=1; k<oc; k++){

			if(openlist[k]->score >= highestscore){

				highest = k;
				highestscore = openlist[k]->score;

			}
		}
		//printf("Check point 1!\n");

		/* Compute BIC of the parent set with highest BIC* */


		/* Save it to closedlist. (Pass the address) */
		if(cc < c-1){
			closedlist[cc] = BICastCrtSet(Xi, openlist[highest]);
		}
		else{
			c = c+c;
			closedlist = realloc(closedlist, c*sizeof(parentSet *));
			closedlist[cc] = BICastCrtSet(Xi, openlist[highest]);
		}
		localscore(Xi, closedlist[cc], metric);
		//printf("Check point 2!\n");

		/*Find all expension of closedlist[cc]*/
		for(n=0; n<Xn; n++){

			unsigned char other = 1; //other == 1 means closedlist[cc] has not been explored yet

			if(n != Xi){

				for(m=0; m<closedlist[cc]->size; m++){

					if(n == closedlist[cc]->parents[m]){

						other = 0;
						break;

					}

				}
				if(other == 1){

					bicAstk *candidate;

					candidate = BICasterisk(Xi, closedlist[cc], closedlist[n], closedlist[Xi]);

					/* literature: C. P. de Campos and Q. Ji,  Efficient structure learning of Bayesian networks using constraints, 2011 */
					double prune = 0.0;
					if( pruneMode == 1 ){
						
						unsigned int v, qc = 1;
						for(v=1; v<=candidate->parlist[0]; v++){

							qc = qc*Nodes[candidate->parlist[v]].statesCount;

						}
						qc = qc * Nodes[n].statesCount;
						
						prune = candidate->score+0.5*lnN*(Nodes[Xi].statesCount-1)*qc;
						
					} 
					//This differs from the original condition in Independence Selection that there is no ``N'' multiplier. Removing ``N'' improves performance.
					if((prune <= closedlist[Xi]->entropy) && (prune <= (openlist[highest]->pis[0])->entropy) && (prune <= (openlist[highest]->pis[1])->entropy)){

						if(closedlist[cc]->size == biggest){ 
				
							sizeincrease = 1;
							if(oc < o-1){
								openlist[oc] = candidate;
								oc+=1; 
							}
							else{
								o = o+o;
								openlist = realloc(openlist, o*sizeof(bicAstk *));
								openlist[oc] = candidate;
								oc+=1;
							}

						}
						else /* closedlist[cc]->size < biggest */{

							unsigned int *list = Union(closedlist[cc], closedlist[n]);

							unsigned int parlist[list[0]];

							unsigned int i;
							for(i=0; i<list[0]; i++){

								parlist[i] = list[i+1];

							}

							unsigned char duplicate = 0;

							unsigned long t;
							for(t=0; t<oc; t++){

								if(openlist[t]->parlist[0] == list[0]){

									if(memcmp(parlist, &openlist[t]->parlist[1], list[0] * sizeof(unsigned int))){
										duplicate = 1;
										break;
									}
								}
							}
							for(t=0; t<cc; t++){

								if(closedlist[t]->size == list[0]){

									if(memcmp(parlist, closedlist[t]->parents, list[0] * sizeof(unsigned int))){
										duplicate = 1;
										break;
									}
								}
							}
							free(list);

							if(duplicate == 0){

								if(oc < o-1){
									openlist[oc] = candidate;
									oc+=1; 
								}
								else{
									o = o+o;
									openlist = realloc(openlist, o*sizeof(bicAstk *));
									openlist[oc] = candidate;
									oc+=1;
								}
							}
						}
					}
					else free(candidate);
				}		
			}
		}
		
		if(sizeincrease == 1) biggest+=1; //!!!

		openlist[highest] = openlist[oc-1]; //just let openlist[highest] point to the last record
		oc-=1; //Note that openlist[oc] still points to something for now, until being overwritten later, should be freed if in last round

		//printf("cc %lu\n oc %llu\n", cc, oc);
		
		cc+=1; //printf("Check point 3. One round done!\n");
		finishing = time(NULL);

	}
	cc-=1;

	/* logging the range of biggest parsent set size explored */
	if (biggest < zlog.psMinMaxSize) zlog.psMinMaxSize = biggest;
	if (biggest > zlog.parsetMaxSize) zlog.parsetMaxSize = biggest; 
	
	if(oc > 0) strcpy(zlog.openlistExhausted, "NO");

	/*free space*/
	unsigned long t;
	for(t=0; t<oc; t++) free(openlist[t]);

	free(openlist);

	closedlist = realloc(closedlist, cc*sizeof(parentSet *));

	/* post-selection pruning */
	qsort(closedlist, cc, sizeof(parentSet *), comparesizes);

	pruner(closedlist, &cc); //requires being sorted by sizes

	closedlist = realloc(closedlist, cc*sizeof(parentSet *));

	/* Sort the cache */
	qsort(closedlist, cc, sizeof(parentSet *), comparescores);

	*cachesize = cc;

	return closedlist;

}
