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

/* Literature: Hirotsugu Akaike, A New Look at the Statistical Model Identification, 1974
	       Gideon Schwarz, Estimating the Dimension of a Model, 1978
	       Jorma Rissanen, Modeling By Shortest Data Description, 1978
	       Joe Suzuki, A Construction of Bayesian Networks from Databases Based on an MDL Principle, 1993
	       David Heckerman, A Tutorial on Learning With Bayesian Networks, 1995 */

void localscore(unsigned int Xi, parentSet *Pii, unsigned char mode){

	unsigned char (*dataA)[Xn] = (unsigned char (*)[Xn])&data[0]; 

	double Nijk[Nodes[Xi].statesCount];

	double ll = 0.0;
	
	double H = 0.0; //entropy

	double ri, qi, Nij, pen;

	ri = Nodes[Xi].statesCount;

	qi = Pii->configsCount;

	// Metrics other than BIC/MDL are not proven working with a guiding metirc yet
	switch(mode){

		case 'A' :

			pen = 1.0;
			break;

		case 'B' :
		default :

			pen = 0.5*lnN;
			break;

		case 'L' :

			pen = 0.0;
			break;

	}

	unsigned char config[Pii->size];
	memset(config, 0, sizeof config); // the starting configuration of the parent set.

	unsigned int c; //counter for running through configurations
	for(c=0; c<(unsigned int) qi; c++){

		Nij = 0.0; //reset

		memset(Nijk, 0, sizeof Nijk); //reset

		unsigned int check;

		/* find out Nij and Nijk from data */
		unsigned long i; //the i here actually means the ith instance/case
		for(i=0; i<N; i++){

			check = 1; //reset

			unsigned int j; //the j here actually means the jth "digit" of a config
			for(j=0; j<Pii->size; j++){

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

		//printf("Nij: %lf\n", Nij);

		/* Add to LL of the current configuration */
		if(Nij != 0.0){
			
			double Pj;
			Pj = Nij/N;
			
			H = H + Pj*log(1.0/Pj);
			
			unsigned char k;
			for(k=0; k<(unsigned int) ri; k++){

				if(Nijk[k] != 0.0){

					ll = ll + Nijk[k]*log(Nijk[k]/Nij);

				}

			}

		}
		//printf("%lf\n", ll);

		/* Next configuration. (variables need not be binary.) */
		unsigned int t;
		for(t=1; t<=Pii->size; t++){

			if( config[Pii->size - t] < Nodes[Pii->parents[Pii->size - t]].statesCount - 1){

				config[Pii->size - t] += 1;
				break;
			}

			else
				config[Pii->size - t] = 0;

		}

	}
	
	Pii->score = ll - pen*(ri-1)*qi;
	Pii->entropy = H;
	//printf("The score, the entropy %lf\n %lf\n logN", Pii->score, Pii->entropy);
	
}

parentSet* emptyParentSet(unsigned int Xi){

	unsigned char (*dataA)[Xn] = (unsigned char (*)[Xn])&data[0]; 

	double  score, ll, H;

	unsigned char ri = Nodes[Xi].statesCount;

	score = 0.0;

	ll = 0.0;

	H = 0.0;

	unsigned int Nik[ri];
	
	memset(Nik, 0, sizeof Nik);

	unsigned long i; 
	for(i=0; i<N; i++){
	
		Nik[dataA[i][Xi]]+=1;

	}

	unsigned char k;
	double p;
	for(k=0; k<ri; k++){

		if(Nik[k] != 0){

			p = (double) Nik[k]/N;

			H = H + ((p)*log(1.0/p));

			ll = ll + (Nik[k]*log(p));

		}

	}

	score = ll - 0.5*lnN*((double) ri-1);

	//printf("The LL, the entropy %lf %lf\n logN is %lf\n", ll, H, lnN);

	parentSet *empty;

	empty = malloc(offsetof(parentSet, parents)+sizeof(unsigned int));

	empty->score = score;
	empty->entropy = H; //It is actually the entropy of Xi, not the empty parent set.
	empty->size = 0;
	empty->configsCount = 1;
	empty->parents[0] = 1729; //Because the size is 0, this value will be ignored.

	return empty;

}

bicAstk* BICasterisk(unsigned int Xi, parentSet *PiL, parentSet *PiR, parentSet *empty){

	double ri, qi1, qi2, emptyXi; 

	qi1 = PiL->configsCount;

	qi2 = PiR->configsCount;

	ri = Nodes[Xi].statesCount;

	emptyXi = empty->score;

	unsigned int *list = Union(PiL, PiR);

	bicAstk *bicXiPiLPiR;

	bicXiPiLPiR = malloc(offsetof(bicAstk, parlist)+(list[0]+1)*sizeof(unsigned int));

	unsigned int i;
	for(i=0; i<=list[0]; i++){

		bicXiPiLPiR->parlist[i] = list[i];

	}

	free(list);

	bicXiPiLPiR->pis[0] = PiL;
	bicXiPiLPiR->pis[1] = PiR;

	/* literature: Scanagatta et al., Learning Bayesian Networks with Thousands of Variables, 2015 */
	bicXiPiLPiR->score = PiL->score + PiR->score + 0.5*lnN*(ri-1)*(qi1+qi2-qi1*qi2-1) - emptyXi;
	//printf("%lf ,", bicXiPiLPiR->score);

	return bicXiPiLPiR;

}

parentSet* BICastCrtSet(unsigned int Xi, bicAstk *bicXiPiLPiR){

	unsigned int parlist[bicXiPiLPiR->parlist[0]];

	unsigned int i;
	for(i=0; i<bicXiPiLPiR->parlist[0]; i++){

		parlist[i] = bicXiPiLPiR->parlist[i+1];

	}

	parentSet *Pi = createParentSet(sizeof(parlist), parlist);

	return Pi;

}


