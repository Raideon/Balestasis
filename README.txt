﻿Balestasis -- efficient Bayesian network learning software.

version 1.03

Copyright 2017 罗良逸 ( Luo Liangyi; ラ リョウイツ ) luoliangyi@gmail.com

This software is distributed under GNU General Public License version 3 (GPLv3). Please check the COPYING file for more information.

All glory to the scientists and engineers who have contributed to the cutting-edge metrics, algorithms, and heuristics this software is based on. (Please check the literature references and include due citations when necessary.)

----------------------------

System requirements: 

	almost any processors (not necessarily x86) that support C language and OpenMP, (The OpenMP requirement can be removed by removing the relevant code that enables multithreading.)

	50MB-8GB, which depends: The memory requirement depends on the scale of the network and the corresponding data set; for reference, learning from a 37-var 10000-instance data set (Alarm network) uses about 50MB; 910-var 12000-instance (c20ng network) 300MB; 2000-var 5000-instance (Random2000-4-5000 network) 5GB.
	
Limitations:

	poor hyper-threading support: please turn off hyper-threading or use a virtual machine where each virtual CPU thread corresponds to a core of the physical processor.  
  
Running:

	Please follow the instructions in the scripts：
	raw_to_net.sh for raw data input to network structure output;
	raw_to_cache.sh for raw data input to cache output (scored parent sets);
	cache_to_net.sh for cache input to network structure output;
	generic.sh for more configurable options.
	
----------------------------

Literature references:

	M. Scanagatta, C. P. de Campos, G. Corani, and M. Zaffalon, Learning Bayesian Networks with Thousands of Variables, 2015 (Featured algorithms: IS, ASOBS)
	
	L. Luo, Learning Bayesian Networks using Fast Heuristics, Master's thesis, 2017 (Featured algorithm: RPGw)

	C. P. de Campos, Z. Zeng, and Q. Ji, Structure Learning of Bayesian Networks Using Constraints, 2009 (Pruning procedure)

	C. P. de Campos and Q. Ji, Efficient Structure Learning of Bayesian Networks Using Constraints, 2011 (Pruning procedure)

	H. Akaike, A New Look at the Statistical Model Identification, 1974 (AIC)

	G. Schwarz, Estimating the Dimension of a Model, 1978 (Schwarz criterion=BIC)

	J. Rissanen, Modeling By Shortest Data Description, 1978 (MDL)

	J. Suzuki, A Construction of Bayesian Networks from Databases Based on an MDL Principle, 1993 (MDL for learning Bayesian networks)

	D. Heckerman, A Tutorial on Learning With Bayesian Networks, 1995 (BIC for learning Bayesian networks)

