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

void saveOrdering(unsigned int *ordering);

void netParaWriterNoParent(unsigned int Xi, FILE *fp);

void netParaWriterNonEmpty(unsigned int Xi, parentSet *Pii, FILE *fp);

void parentSetDisplay(parentSet *set, unsigned char q, unsigned int c, unsigned long index);

void anceTracer(parentSet **net, unsigned int *list, unsigned int *c, unsigned int *ordering, unsigned int i, unsigned int *ra);

void shuffler(unsigned int *ordering);


void cachesToFile(parentSet ***cacheofscores, unsigned long *cachesizes){
	
	FILE *fp;
	fp = fopen(outputPath, "a+");
	
	fprintf(fp, "%u\n", Xn);
	
	unsigned long i,j,k;
	for(i=0; i<Xn-1; i++) fprintf(fp, "%s ", Nodes[i].name);
	fprintf(fp, "%s\n", Nodes[Xn-1].name);
	
	for(i=0; i<Xn-1; i++) fprintf(fp, "%hhu ", Nodes[i].statesCount);
	fprintf(fp, "%hhu\n", Nodes[Xn-1].statesCount);
	
	for(i=0; i<Xn; i++){
		
		fprintf(fp, "%lu %lu\n", i, cachesizes[i]);
		
		for(k=0; k<cachesizes[i]; k++){
	
			fprintf(fp, "%lf %u ", cacheofscores[i][k]->score, cacheofscores[i][k]->size);
		
			for(j=0; j<cacheofscores[i][k]->size; j++) fprintf(fp, " %u", cacheofscores[i][k]->parents[j]);
			
			fputs("\n", fp);
			
		}
	
	}
	fclose(fp);
	
}

void structureToFile(){
	
	FILE *fp;
	fp = fopen(outputPath, "a+");
	
	fprintf(fp, "%u\n", Xn);
	
	unsigned int i,j;
	for(i=0; i<Xn; i++){
		
		fprintf(fp, "%u\n%lf %u", Nodes[resultorder[i]].datumAC, finalNet[i]->score, finalNet[i]->size);
		
		for(j=0; j<finalNet[i]->size; j++) fprintf(fp, " %u", finalNet[i]->parents[j]);
			
		fputs("\n", fp);
		
	}
	
	fclose(fp);
	
}

void netToDSCFile(){
	
	FILE *fp;
	fp = fopen(outputPath, "a+");
	
	fputs("network \"unknown\" { version is 1; }\n\n", fp);
	
	unsigned int i, j;
	for(i=0; i<Xn; i++){
		
		fprintf(fp, "node %s {\n  type : discrete [ %hhu ] = { " , Nodes[i].name, Nodes[i].statesCount);
		
		for(j=0; j<Nodes[i].statesCount-1; j++){ 
		
			fprintf(fp, "\"%u\", ", j); 
				
		}
		fprintf(fp, "\"%u\" };\n}\n", Nodes[i].statesCount-1); 
		
	}
	
	for(i=0; i<Xn; i++){
		
		if(finalNet[i]->size == 0) netParaWriterNoParent(resultorder[i], fp);
		else
			netParaWriterNonEmpty(resultorder[i], finalNet[i], fp);
	
	}
	fclose(fp);
	
}

void netParaWriterNoParent(unsigned int Xi, FILE *fp){
	
	fprintf(fp, "probability ( %s ) {\n   ", Nodes[Xi].name);

	unsigned char (*dataA)[Xn] = (unsigned char (*)[Xn])&data[0]; 

	unsigned char ri = Nodes[Xi].statesCount;

	unsigned int Nik[ri];
	
	memset(Nik, 0, sizeof Nik);

	unsigned long i; 
	for(i=0; i<N; i++){
	
		Nik[dataA[i][Xi]]+=1;

	}

	unsigned char k;
	double p;
	for(k=0; k<ri-1; k++){

		p = (double) Nik[k]/N;
		fprintf(fp, "%lf, ", p);

	}
	fprintf(fp, "%lf;\n}\n", (double) Nik[ri-1]/N);

}



void netParaWriterNonEmpty(unsigned int Xi, parentSet *Pii, FILE *fp){ 

	fprintf(fp, "probability ( %s | ", Nodes[Xi].name);
	
	unsigned int v, parsize = Pii->size;
	for(v=0; v<(parsize-1); v++){
		
		fprintf(fp, "%s, ", Nodes[Pii->parents[v]].name);
		
	}
	fprintf(fp, "%s ) {\n", Nodes[Pii->parents[v]].name);

	unsigned char (*dataA)[Xn] = (unsigned char (*)[Xn])&data[0]; 

	double Nijk[Nodes[Xi].statesCount];

	double ri, qi, Nij;

	ri = Nodes[Xi].statesCount;

	qi = Pii->configsCount;

	unsigned char config[parsize];
	memset(config, 0, sizeof config); // the starting configuration of the parent set.

	unsigned int c; //counter for running through configurations
	for(c=0; c<(unsigned int) qi; c++){
		
		fputs("  (", fp);
		
		unsigned int q;
		for(q=0; q<(parsize-1); q++){
			
			fprintf(fp, "%hhu, ", config[q]);
			
		}
		fprintf(fp, "%hhu) : ", config[q]);
		
		Nij = 0.0; //reset

		memset(Nijk, 0, sizeof Nijk); //reset

		/* find out Nij and Nijk from data */
		unsigned long i; //the i here actually means the ith instance/case
		for(i=0; i<N; i++){

			unsigned char check = 1; //reset

			unsigned int j; //the j here actually means the jth "digit" of a config
			for(j=0; j<parsize; j++){

				if(config[j] != dataA[i][Pii->parents[j]]){

					check = 0;
					break;

				}
			}

			if(check == 1){

				Nij+=1.0;
				
				Nijk[dataA[i][Xi]]+=1.0;
				
			}

		}
		
		unsigned char k;
		unsigned int riminusOne = (unsigned int) ri-1;
		for(k=0; k<riminusOne; k++){

				
			if(Nij == 0) fprintf(fp, "%lf, ", 1.0/ri);

			else fprintf(fp, "%lf, ", Nijk[k]/Nij);

		}
		if(Nij == 0) fprintf(fp, "%lf;\n", 1.0/ri);
		else fprintf(fp, "%lf;\n", Nijk[k]/Nij);
		
		/* Next configuration. (variables need not be binary.) */
		unsigned int t;
		for(t=1; t<=parsize; t++){

			if( config[parsize - t] < Nodes[Pii->parents[parsize - t]].statesCount - 1){

				config[parsize - t] += 1;
				break;
			}

			else
				config[parsize - t] = 0;

		}

	}
	fputs("}\n", fp);
	
}

/* This can be used to save the orderings for analysis or whatever.*/
void saveOrdering(unsigned int *ordering){ 
	
	FILE *fp;
	fp = fopen("/home/luo/c/saved_ordering.dat", "a+");
	
	unsigned int i;
	for(i=0; i<Xn; i++){ 
	
		fprintf(fp, "%u", ordering[i]);
		if(i<Xn-1) fputs(" ", fp);
		
	}
	fputs("\n", fp);
	fclose(fp);
	
}


void shuffler(unsigned int *ordering){
 
	unsigned int i;
	for(i=0; i<Xn-2; i++){

		unsigned int rj = i + rand() % (Xn-i);
		unsigned int swapper = ordering[rj];
		ordering[rj] = ordering[i];
		ordering[i] = swapper;
	}
}

void swapper(unsigned int *ordering, unsigned int *rj){

	unsigned int swapper = ordering[*rj];
	ordering[*rj] = ordering[*rj+1];
	ordering[*rj+1] = swapper;

}

/* searching all ancestors using DFS with the part of network being already constructed
   literature: Scanagatta et al., Learning Bayesian Networks with Thousands of Variables, 2015 */
void anceTracer(parentSet **net, unsigned int *list, unsigned int *c, unsigned int *ordering, unsigned int i, unsigned int *ra){

	unsigned int k;				
	for(k=0; k<net[i]->size; k++){

		unsigned int z, add = 1;
		for(z=0; z<*c; z++){
				
			if(net[i]->parents[k] == list[z]){ add = 0; break; }
				
		}
		if(add == 1){
			
			list[*c] = net[i]->parents[k];
			*c+=1;

			unsigned int t, p;
			for(t=0; t<*ra; t++){
					
				if(net[i]->parents[k] == ordering[t]){
				
					anceTracer(net, list, c, ordering, t, ra);
			
				}													
			}
		}
	}
}


parentSet* createParentSet(unsigned int s, unsigned int *parentlist){

	struct ParentSet *parentsAddr;

	unsigned int a;

	a = offsetof(parentSet, parents)+s; //because it looks clearer

 	parentsAddr = malloc(a); //because it looks clearer

	s = s/sizeof(unsigned int);

	parentsAddr->size = s;

	unsigned int q = 1;

	unsigned int i;
	for(i=0; i<s; i++){

		parentsAddr->parents[i] = parentlist[i];
		q = q*Nodes[parentlist[i]].statesCount;

	}
	parentsAddr->configsCount = q;

	return parentsAddr;

}

void nodeDisplay(unsigned int i){

	printf("Name: %s\nthe number of discrete states: %hhu\ncodename: %u\n", Nodes[i].name, Nodes[i].statesCount, Nodes[i].datumAC);
	
}

/* This also provides basic check against overflow, corruption etc. */
void parentSetDisplay(parentSet *set, unsigned char q, unsigned int c, unsigned long index){
	
	if(set->size > 100 ) { printf("Warning! Anomaly detected. Parent set %lu from node %u size\n", index, c ); q = 1; }
	if(set->score >= 0 || set->score < -99999999 || isnan(set->entropy) ) { printf("Warning! Anomaly detected. Parent set %lu from node %u score\n", index, c ); q = 1; }
	if(set->entropy <= 0 || set->entropy > 20 || isnan(set->entropy) ) { printf("Warning! Anomaly detected. Parent set %lu from node %u entropy\n", index, c ); q = 1; }

	if(q == 1){
		
		printf("Parent set size is %u\n", set->size);
		printf("Parent set score is %lf\n",set->score);
		printf("Parent set entropy is %lf\n",set->entropy);
		printf("This parent set has the following parents.\n");
		unsigned int i;
		for(i=0; i<set->size; i++) printf("Parent (codename): %u\n", set->parents[i]);
		printf("Address offsets of the parent set.\n");
		printf("%p\n", &set->score);
		printf("%p\n", &set->size);
		printf("%p\n", &set->configsCount);
		printf("%p\n", &set->parents);
		printf("%p\n", &set->parents[set->size-1]);

	}

}

/* In this function, I have implemented the algorithm for solving the union set of two sets. 
   It might be more than necessary, hence sets grow incrementally.
   It is not computationally expensive nonetheless. */
unsigned int* Union(parentSet *set1, parentSet *set2){

	unsigned int *s1Unions2;

	s1Unions2 = malloc((set1->size+set2->size+1)*sizeof(unsigned int));

	unsigned int i=0,j=0,c=1;

	while(i<set1->size && j<set2->size){

		if(set1->parents[i] == set2->parents[j]){

			s1Unions2[c] = set1->parents[i];
			i++; j++; c++;
		
		}
		else if(set1->parents[i] < set2->parents[j]){

			s1Unions2[c] = set1->parents[i];
			i++; c++;

		}
		else{

			s1Unions2[c] = set2->parents[j];
			j++; c++;

		}

	}
	for(; i<set1->size; i++){

		s1Unions2[c] = set1->parents[i];
		c++;

	}
	for(; j<set2->size; j++){

		s1Unions2[c] = set2->parents[j];
		c++;
		
	}
	s1Unions2[0] = c-1;

	return s1Unions2; //s1Unions2[0]: the number of elements

}

int comparescores(const void *x, const void *y){

	parentSet **a = (parentSet **)x;
	parentSet **b = (parentSet **)y;

	return (*a)->score < (*b)->score ? 1 : -1 ;

}

int comparesizes(const void *x, const void *y){

	parentSet **a = (parentSet **)x;
	parentSet **b = (parentSet **)y;

	return (*a)->size > (*b)->size ? 1 : -1 ;

}


/* If a set is worse than its subset, prune it.
   literature: C. P. de Campos, Z. Zeng, and Q. Ji, Structure Learning of Bayesian Networks Using Constraints, 2009*/
void pruner(parentSet **list, unsigned long *cc){

	unsigned long last, current;

	current = *cc - 1;
	last = current;

	while(current > 1){

		unsigned int i;
		for(i=0; i<current; i++){

			if(list[i]->score > list[current]->score){

				unsigned char subset = 1;

				unsigned int k,l;
				for(k=0; k<list[i]->size; k++){

					for(l=0; l<list[current]->size; l++){

						if(list[i]->parents[k] == list[current]->parents[l]){subset = 1; break;}
						subset = 0;

					}
					if(subset == 0) break;

				}
				if(subset == 1){

					list[current] = list[last]; 
					last-=1;
					break;

				}
			}
		}	
		current-=1;
	}
	*cc = last;
}

void diag(){
	
	/* Should this project continue its growth, 
	somebody in the future (including you, the future me) who sees this might as well write a diagnosis routine. */
	
}
