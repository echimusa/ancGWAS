#!/usr/bin/python

from __future__ import with_statement
from __future__ import division
import numpy as np
import os, sys, fileinput, exceptions, types, time,datetime
import logging, pickle, pylab
import subprocess
import warnings
import rpy2.robjects as R0
import itertools as iter
from math import sqrt,log
from networkx import *
from scipy.optimize import leastsq
from scipy import power

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger.addHandler(console)
warnings.filterwarnings("ignore")
r = R0.r 

class ancUtils:
    '''
    Useful functions
    '''
    def __init__(self):
        return None
      
    def logical(self, Valeur):
        '''
        Convert a string to logical variable
        Valeur: a string, NO or YES
        Return a logical variable
        '''
        if Valeur == "YES":
            return True
        else:
            return False
    
    def safe_open(self,filename):
        '''
        Opens a file in a safe manner or raise an exception if the file doesn't exist
        filename: a string, name of a file
        Return: file object
        '''
        if os.path.exists(filename):
            return open(filename)
        else:
            self.printHead()
            sys.stderr.write("Error => No such file or directory\"%s\"\n"%format(filename))
            return False
            
    def check_params(self,params):
        '''
        Check necessary parameters in a dict
        params: a dictionary, of option:value 
        '''
        for key in params:
            if params[key][0] in ['Faffy','Fanc','Fgeno','Fgwas','Fnetwork','Fsnp']:
                if os.path.exists(params[key][1]) == False: 
                    sys.stderr.write("Error => No such file or directory \"%s\"\n\n" % params[key][-1])
                    sys.exit(1)
        
    def list_comprehesive(self, myList):
        '''
        Compute the average of multiple lists
        Input: List of lists of numeric values
        '''
        tmp = []
        for de in range(1,len(myList[0])):
            tm = 0
            for des in range(len(myList)):
                tm = tm + myList[des][de] 
            tmp.append(tm/len(myList))
        return tmp
    
    def complex_list_comprehesive(self, myList):
        '''
        '''
        tmp = []
        tmp1 = []
        n = float(len(myList))
        for de in range(len(myList[0][0])):
            tm = 0
            tm1 = 0
            ts = 0
            for des in range(len(myList)):
                tm = tm + myList[des][0][de] 
                tm1 = tm1 + myList[des][1][de] 
                ts =ts+(myList[des][0][de]**2)
            tmp.append(tm/n)
            tmp1.append(ts/(n*tm1))
        return tmp,tmp1,n-1
    
    def terminate(self):
        '''
        Terminate the process
        '''
        try:
            subprocess.Popen(['mv', self.logFile, self.outfolder])
            ##pkl_files = [fil for fil in os.listdir(self.outfolder) if fil.endswith('pkl')]
            ##for fil in pkl_files:
            ##  os.remove(self.outfolder+fil)
            log = self.outfolder+self.logFile
        except:
            log = os.getcwd()+'/'+self.logFile
        finally:
            print'Log file generated in '+log+'\nHave a nice day!\n'
            sys.exit(1)

class ancInit(ancUtils):
    '''
    Initialize ancGWAS by 
        - creating the output folder
        - Reading the parameter file
    '''
    
    try:
        logger.removeHandler(logger.handlers[1])
    except:
        pass
    
    logger.setLevel(logging.INFO)
    logFile = 'ancGWAS-'+str(time.time()).split('.')[0]+'.log'
    fh = logging.FileHandler(logFile, mode='w')
    logger.addHandler(fh)
    
    logger.info("\n********************************************************************************************************************")
    logger.info("     ancGWAS: Searching Subnetwork Underlying Ethnic Difference in Genetic Risk in  recently admixed population       ")
    logger.info("                                            Computational Biology Group                                               ")
    logger.info("                                    2012, University of Cape Town, South Africa                                       ")
    logger.info("                                                Verson 1.0 Beta                                                       ")
    logger.info(  "********************************************************************************************************************\n")

    
    def __init__(self, argv, logFile=logFile):
        '''
        Initializing ancGWAS by reading the parameter file
        '''
        self.argv = [argv]
        self.logFile = logFile
        popnames = []
        if len(self.argv) == 0 or self.argv == ['']:
            logger.info('Command line usage: %s <parameter file>  ' % sys.argv[0])
            logger.info('eg: python ancGWAS.py parancGWAS.txt\n')
            sys.exit(1)
        elif len(self.argv) == 1:
            try:
                self.paramFile = self.argv[0]
                if os.path.exists(os.getcwd()+'/'+self.paramFile):
                    inFile = open(os.getcwd()+'/'+self.paramFile)
                elif os.path.exists(self.paramFile):
                    inFile = open(self.paramFile)
                else:
                    logger.info('\nError => Failed to process the input, check the parameters!\n')
                    self.terminate()
                rows = {}
                self.opts = []
                self.opts1 = {}
                myList1 = []
                myList = ['anced','stouffer','wzscore','tklc']
                myList2 = ['LDcutoff','Pcutoff','SNPtoGENE','Netpathstep','TOPScore']
                    
                for line in inFile:
                    data = line.split()
                    if not line.startswith("#"):
                        if len(data) != 0:
                            rows[data[0].split(":")[0]] = data[0].split(":")[1]
                            data1 = data[0].split(":")[0]
                            data2 = data[0].split(":")[1]
                            if data1 in myList:
                                data2 = self.logical(data2)
                            self.opts.append(data2)
                            self.opts1[data1] = data2
                            myList1.append(data1)
                            
                if len(self.opts1) != 19:
                    sys.stderr.write('\nMissing parameters!!\nError => Failed to process the input, check the parameters!\n')
                    self.terminate()
                try:
                    self.outfolder = self.opts1['outfolder']
                    path = os.getcwd()
                    if os.path.exists(self.outfolder):
                        pass
                    else:
                        os.makedirs(self.outfolder)
                except IndexError:
                    sys.stderr.write('Error => Can not create directory. Please create a directory OUT in your working directory\n')
                    sys.exit(1)

                for par in self.opts1:
                    if par in myList2:
                        try:
                            self.opts1[par] = float(self.opts1[par])
                        except:
                            sys.stderr.write('Error => Failed to process the input, check the parameters!\n')
                            self.terminate()
                            sys.exit(1)
                        finally:
                            if (par=='LDcutoff' or par=='Pcutoff') and (0>=self.opts1[par] or self.opts1[par]>=1):
                                sys.stderr.write('Invalid LD or Pvalue cutoff; Must range between 0 and 1!\n')
                                sys.stderr.write('Error => Failed to process the input, check the parameters!\n\n')
                                self.terminate()
                                sys.exit(1)
                            if par == 'SNPtoGENE' and self.opts1[par]<=0:
                                sys.stderr.write('Invalid step path; Must be greater than 0!\n')
                                sys.stderr.write('Error => Failed to process the input, check the parameters!\n\n')
                                self.terminate()
                                sys.exit(1)
                allList = ['Fgwas', 'Faffy', 'Fgeno', 'Fsnp', 'Fnetwork', 'LDcutoff', 'Pcutoff', 'SNPtoGENE', 'Gene_pv', 'anclabel', 'ancprop', 'anced', 'Gene_LD', 'stouffer', 'wzscore', 'Netpathstep', 'TOPScore', 'tklc', 'outfolder']
                List1 = []
                for i in range(len(myList1)):
                    if myList1[i] != allList[i]:
                        List1.append(myList1[i])
                if len(List1) != 0:
                    sys.stderr.write('\nInvalide option: '+','.join(List1))
                    sys.stderr.write('\nError => Failed to process the input, check the parameters!\n')
                    self.terminate()
                    sys.exit()
                    
                for param in myList:
                    rows[param] = self.logical(rows[param])
                i = 1
                self.Params = {}
                for param in myList1:
                    if param in ['anced'] or rows[param] != False:
                        self.Params[i] = [param,rows[param]]
                        i+=1
                self.check_params(self.Params)
                self.assocFile = self.opts1['Fgwas']
                self.affyFile = self.opts1['Faffy']
                self.genoFile = self.opts1['Fgeno']
                self.snpFile = self.opts1['Fsnp']
                self.networkFile = self.opts1['Fnetwork']
                self.Pcutoff = float(self.opts1['Pcutoff']) 
                self.boundary =  int(self.opts1['SNPtoGENE']) 
                self.step = int(self.opts1['Netpathstep'])
                self.TOPScore = int(self.opts1['TOPScore'])
                self.LDCutoff = int(self.opts1['LDcutoff'])
                
                Gene_pvs = ['fisher', 'simes', 'smallest', 'gwbon']
                self.Gene_pv = self.opts1['Gene_pv']
                if self.Gene_pv.lower() not in Gene_pvs:
                    sys.stderr.write('Invalid option for gene pvalue option Gene_pv:'+self.Gene_pv+'\n'+'Error => Failed to process the input, check the parameters!\n\n')
                    self.terminate()
                    sys.exit(1)
                self.anc_label = self.opts1['anclabel'].split(",") 
                self.anc_prop = self.opts1['ancprop']
                self.anced = self.opts1['anced']
                Gene_LDs = ['closest', 'zscore', 'maxscore']
                self.Gene_LD = self.opts1['Gene_LD']
                if self.Gene_LD.lower() not in Gene_LDs:
                    sys.stderr.write('Invalid option for gene-LD option:'+self.Gene_LD+'\n'+'Error => Failed to process the input, check the parameters!\n\n')
                    self.terminate()
                    sys.exit(1)
                self.lampld = self.opts[12]
                self.stouffer = self.opts1['stouffer']
                self.wzscore = self.opts1['wzscore']
                self.tklc = self.opts1['tklc']

            except IndexError, TypeError:
                sys.stderr.write('Error => Failed to process the input, check the parameters!\n\n')
                self.terminate()
                sys.exit(1)
        else:
            logger.info('Command line usage: %s <parameter file>  ' % sys.argv[0])
            logger.info('eg. python ancGWAS.py parancGWAS.txt\n')
            self.terminate()
            sys.exit(1) 

class ancmap():
    '''
    Mapping module, which maps snps to genes and computes the p-value at gene level.
    '''
    
    def readGWAS(self):
        '''
        Reads the GWAS file of format SNP - Pv - Ancentryy proportions if available
        '''
        self.assoc_map = {}    
        logger.info('Reading ... %s' % self.assocFile)
        for line in fileinput.input(self.assocFile): 
            data = line.strip().split() 
            if self.anced:
                if len(data) < 2:
                    sys.stderr.write('\nError => Your GWA file does not match with your parameters\n')
                    self.terminate()
                if len(self.anc_label) != len(data[2:]):
                    sys.stderr.write('\nError => Your ancestry labels in your file don\'t match with your parameters\n')
                    self.terminate()
            else:
                if len(data) > 2:
                    sys.stderr.write('\nError => Your GWA file does not match with your parameters. If your file contains ancestry specific information, please activate the option "\anced\"\n')
                    self.terminate()
            if fileinput.lineno() == 1:
                if data[0] != 'SNP' or data[1] != 'P' :
                        logger.info('\nAssociation file not in valide format!')
                        logger.info('The first row should be \'SNP P\', should only contains the rsID and pvalue, and the file should be SPACE-delimited or TAB-delimited')
                        self.terminate()
            else:
                if len(data) != 0:    
                    snp = data[0]
                    for i in data[1:]:
                        if 0<float(i)>1:
                            logger.info('\nInvalid value in your GWA dataset!\nError => Failed to process the input!\n')
                            self.terminate()
                            sys.exit(1)
                    if self.anced == True:
                        self.assoc_map[snp] = map(float,[data[i] for i in range(1,len(data))]) # Accounting for PV
                    else:
                        self.assoc_map[snp] = [float(data[1])] 
        logger.info('Including ... %s SNPs with associated Pvalue' % str(len(self.assoc_map)))
        output = open(self.outfolder+'assoc_map.pkl','wb')
        pickle.dump(self.assoc_map, output)
        output.close()
    
    def prereadAffy(self):
        '''
        Read Affymetrix file and return the min distance of the the smallest pvalue in Assoc file
        '''
        logger.info('\nAjusting the boundary distance')
        minDist = 100000000
        for line in fileinput.input(self.affyFile):
            ## data1 = AffyID, SNP, Dist, Gene
            data1 = line.split() 
            if fileinput.lineno()>1:
                if len(data1) == 4:
                    snp = data1[1]; gene = data1[-1]; dist = data1[-2]
                elif len(data1) == 3:
                    snp = data1[0]; gene = data1[-1]; dist = data1[-2]                 
                if snp == self.minSNP:
                    if minDist > float(dist):
                        minDist = float(dist)
        if 100000 > minDist < 1000 * self.boundary:
            pass
        else:
            self.boundary = minDist/1000
        ##print 1000 * self.boundary, self.boundary, minDist 
    
    
    def readAffy(self):
        '''
        Reads the affymetrix File of format AffyID - SNP - Dist - Gene
        '''
        gene_map = {} 
        colnames_affy = [] 
        affy_map1 = {} 
        affy_map2 = {} 
        affy_map3 = {} 
        gene_map1 ={}        
        file1 = open(self.outfolder+"gene_map.txt","wt") #Gene - SNP - anc(snp)
        if self.anced: 
            file1.writelines("Gene"+"\t"+"SNP"+"\t"+"P"+"\t"+"\t".join(self.anc_label)+"\n")
        else: 
            file1.writelines("Gene"+"\t"+"SNP"+"\t"+"P"+"\n")
        logger.info('Reading ... %s' % self.affyFile)
        logger.info('Reading ... %s SNPs' % str(len(self.assoc_map)))
        for line in fileinput.input(self.affyFile):
            data1 = line.split() 
            if len(data1) == 4:
                snp = data1[1]; gene = data1[-1]; dist = data1[-2]                 
                if snp in self.assoc_map:
                    if gene in affy_map1:
                        f2 = affy_map2[gene]
                        f2.append([snp,float(dist)])
                        if float(dist) < 1000 * self.boundary:
                            f = gene_map[gene]
                            d(self.assoc_map[snp][0])
                            affy_map1[gene]
                            dd(snp)
                        elif float(self.assoc_map[snp][0]) < self.Pcutoff:
                            f = gene_map[gene]
                            f.add(self.assoc_map[snp][0])
                    else:
                        affy_map2[gene]=[[snp,float(dist)]] 
                        if float(dist) < 1000 * self.boundary:
                            gene_map[gene] = set([self.assoc_map[snp][0]])
                            affy_map1[gene]=set([snp])
                            file1.writelines(str(gene)+"\t"+str(snp)+"\t"+"\t".join([str(p) for p in self.assoc_map[data1[1]]])+"\n")
                        elif float(self.assoc_map[snp][0]) < self.Pcutoff:
                            gene_map[gene]=set([self.assoc_map[snp][0]])
                            file1.writelines(str(gene)+"\t"+str(snp)+"\t"+"\t".join([str(p) for p in self.assoc_map[data1[1]]])+"\n")
            elif len(data1) == 3:
                if fileinput.lineno()>1:
                    snp = data1[0]; gene = data1[-1]; dist = data1[-2]                 
                    if snp in self.assoc_map:
                        if gene in affy_map1:
                            f2 = affy_map2[gene]
                            f2.append([snp,float(dist)])
                            if float(dist) < 1000 * self.boundary:
                                f = gene_map[gene]
                                f.add(self.assoc_map[snp][0])
                                f1 = affy_map1[gene]
                                f1.add(snp)
                            elif float(self.assoc_map[snp][0]) < self.Pcutoff:
                                f = gene_map[gene]
                                f.add(self.assoc_map[snp][0])
                        else :
                            affy_map2[gene]=[[snp,float(dist)]] 
                            if float(dist) < 1000 * self.boundary:
                                gene_map[gene]=set([self.assoc_map[snp][0]])
                                affy_map1[gene]=set([snp])
                                file1.writelines(str(gene)+"\t"+str(snp)+"\t"+"\t".join([str(p) for p in self.assoc_map[data1[0]]])+"\n")
                            elif float(self.assoc_map[snp][0]) < self.Pcutoff:
                                gene_map[gene]=set([self.assoc_map[snp][0]])
                                file1.writelines(str(gene)+"\t"+str(snp)+"\t"+"\t".join([str(p) for p in self.assoc_map[data1[0]]])+"\n")
                else:
                    if data1[0] != "SNP":
                        logger.info('Affymetrix file not in valide format!')
                        logger.info('Should include either AffyID-SNP-Dist-Gene or SNP-Dist-Gene columns')
                        self.terminate()
                        raise SystemExit
        affy_map4 = {}
        for gene in affy_map2:
            tmp = sorted(affy_map2[gene], key=lambda a_entry: float(a_entry[-1]))
            affy_map3[gene] = tmp[0]
        output = open(self.outfolder+'affy_map1.pkl','wb')
        pickle.dump(affy_map1, output)
        output = open(self.outfolder+'affy_map3.pkl','wb')
        pickle.dump(affy_map3, output)
        output = open(self.outfolder+'gene_map.pkl','wb')
        pickle.dump(gene_map, output)
        file1.close(), output.close()
        logger.info('Genes mapped ... %s Genes' % str(len(gene_map)))
        self.gene_map = gene_map
        self.affy_map1 = affy_map1
    
    def simes(self):
        '''
        Computes the gene p-value using the simes method.
        '''
        gene2Pv = {} 
        gene2Pv_ = [] 
        file1 = open(self.outfolder+"gene2Pv.txt","wt")
        file1.writelines("## P-value for the gene using Simes Method ##\n")
        file1.writelines("Gene"+"\t"+"Pv"+"\n")
        for gene in self.gene_map:
            gene_p = [float(i) for i in list(self.gene_map[gene])] 
            tmp = sorted(gene_p)    
            pos = range(1,len(tmp)+1)
            tmp2 = (np.array(tmp)*(len(tmp)+1))/np.array(pos)
            gene2Pv[gene] = min(tmp2)
            gene2Pv_.append([gene, min(tmp2)])
        
        gene2Pv_ = sorted(gene2Pv_, key= lambda a_key: float(a_key[1]), reverse=False)
        for lst in gene2Pv_:
            file1.writelines(lst[0]+"\t"+str(lst[1])+"\n")
        file1.close()
        output = open(self.outfolder+'gene2Pv.pkl','wb')
        pickle.dump(gene2Pv, output)
        gene2Pv = {}
        
    def fisher(self):
        '''
        Computes the gene p-value using fisher's method.
        '''
        gene2Pv = {} 
        gene2Pv_ = [] 
        file1 = open(self.outfolder+"gene2Pv.txt","wt")
        file1.writelines("## P-value for the gene using Fisher's Method ##\n")
        file1.writelines("Gene"+"\t"+"Pv"+"\n")
        for gene in self.gene_map:
            gene_p = [float(i) for i in list(self.gene_map[gene]) if float(i)!=0] 
            chi2 = sum(-2*np.log(np.array(gene_p)))
            t = list(r.pchisq(chi2,2*len(gene_p)))
            p = 1 - t[0]
            gene2Pv[gene] = p
            gene2Pv_.append([gene, p])
        
        gene2Pv_ = sorted(gene2Pv_, key= lambda a_key: float(a_key[1]), reverse=False)
        for lst in gene2Pv_:
            file1.writelines(lst[0]+"\t"+str(lst[1])+"\n")
        file1.close()
        output = open(self.outfolder+'gene2Pv.pkl','wb')
        pickle.dump(gene2Pv, output)
        gene2Pv = {}
        
    def smallest(self):
        '''
        Computes the gene p-value using the smallest p-value of all the snps in a given gene.
        '''
        gene2Pv = {} 
        gene2Pv_ = [] 
        file1 = open(self.outfolder+"gene2Pv.txt","wt")
        file1.writelines("## P-value for the gene using Smallest Method ##\n")
        file1.writelines("Gene"+"\t"+"Pv"+"\n")
        for gene in self.gene_map:
            gene_p = [float(i) for i in list(self.gene_map[gene])] 
            gene2Pv[gene] = min(gene_p)
            gene2Pv_.append([gene, min(gene_p)])
        gene2Pv_ = sorted(gene2Pv_, key= lambda a_key: float(a_key[1]), reverse=False)
        for lst in gene2Pv_:
            file1.writelines(lst[0]+"\t"+str(lst[1])+"\n")
        file1.close()
        output = open(self.outfolder+'gene2Pv.pkl','wb')
        pickle.dump(gene2Pv, output)
        gene2Pv = {}
        
    def FDR(self):
        '''
        Computes the gene p-value using the gene-wise FDR method.
        '''
        gene2Pv = {} 
        gene2Pv_ = [] 
        file1 = open(self.outfolder+"gene2Pv.txt","wt")
        file1.writelines("### Computing the p-value for the gene using gene-wise FDR value ###\n")
        file1.writelines("Gene"+"\t"+"Pv"+"\n")
        for gene in self.gene_map:
            gene_p = [float(i) for i in list(self.gene_map[gene])] 
            tmp = list(R0.r['p.adjust'](R0.FloatVector(gene_p),method='bonferroni'))    
            gene2Pv[gene] = min(tmp)
            gene2Pv_.append([gene, min(tmp)])
        gene2Pv_ = sorted(gene2Pv_, key= lambda a_key: float(a_key[1]), reverse=False)
        for lst in gene2Pv_:
            file1.writelines(lst[0]+"\t"+str(lst[1])+"\n")
        file1.close()
        output = open(self.outfolder+'gene2Pv.pkl','wb')
        pickle.dump(gene2Pv, output)
        gene2Pv = {} 

class ancLD():
    '''
    Estimates the egde weight using the LD between snps in genes and constructs the LD-weighted PPI network.
    '''
    def __init__(self):
        return None
    
    def delta(self):
        '''
        Estimating diffrence in defficiency/excess of ancestry
        '''
        self.anc_prop = self.anc_prop.split(',')
        rows = {}
        logger.info('Estimating defficiency and excess ancestry ...')
        for line in fileinput.input(self.assocFile):
            data = line.split()
            tmp1 = []
            if fileinput.lineno() > 1 :
                if len(self.anc_prop) == len(data[2:]) or sum([ float(i) for i in self.anc_prop]) == 1.0:
                    tmp = [float(d)-float(self.anc_prop[data[2:].index(d)]) for d in data[2:]]
                    for k in tmp:
                        z = r.qnorm(1-(abs(k/2)))[0]
                        tmp1.append(abs(k/float(z)))
                    rows[data[0]]=[tmp,tmp1]
                else:
                    print >> sys.stderr, "Error: Ancestral proportion not sum to 1.", "or NO ancestral proportion provided"
                    print >> sys.stderr, "Please check parameter setting!"
                    self.terminate()
                    sys.exit(1)
            else:
                anc = data[2:]
        logger.info('Including ... %s SNPs with associated delta ancestry' % str(len(rows)))
        output = open(self.outfolder+'Delta.pkl','wb')
        pickle.dump(rows,output)
        rows.clear()
        
    def Gene2Anc(self):
        '''
        Computing the gene ancestry proportions
        '''
        self.Gene2snpAnc = {} 
        Gene2delta = {}
        logger.info('Loading Snp - Delta Ancestry and its variances ...')
        delta_file = open(self.outfolder+'Delta.pkl', 'rb')
        G = pickle.load(delta_file)
        logger.info('Loading Gene - Ancestry proportions ...')
        pkl_file = open(self.outfolder+'affy_map1.pkl', 'rb')
        self.affy_map1 = pickle.load(pkl_file)
        pkl_file = open(self.outfolder+'assoc_map.pkl', 'rb')
        self.assoc_map = pickle.load(pkl_file)
        for gene in self.affy_map1:
            tmp = [self.assoc_map[snp] for snp in self.affy_map1[gene]]
            self.Gene2snpAnc[gene] = [self.assoc_map[snp] for snp in self.affy_map1[gene]]
            Gene2delta[gene] = [G[snp] for snp in self.affy_map1[gene]]
        output = open(self.outfolder+'Gene2snpAnc.pkl','wb')
        pickle.dump(self.Gene2snpAnc, output)
        output = open(self.outfolder+'Gene2snpDelta.pkl','wb')
        pickle.dump(Gene2delta, output)
    
    def Weightedgene2Anc(self):
        '''
        Input: merge_map = [Gene , SNP - Pv of SNP , Anc1 , Anc2 , Anc3 ------]
        Output: list gene with propotions of each ancestry (Gene-Anc1-Anc2-----)
        '''
        file1 = open(self.outfolder+"gene2AncWeight.txt","wt")
        gene2AncWeight = {}
        file1.write("GENE"+"\t"+"\t".join(self.anc_label)+"\n")
        for gene in self.Gene2snpAnc:
            tmp = self.list_comprehesive(self.Gene2snpAnc[gene])
            c = list(tmp)
            gene2AncWeight[gene] = c
            file1.writelines(str(gene))
            for i in c:
                file1.writelines("\t"+str(i))
            file1.writelines("\n")
        file1.close()
        output7 = open(self.outfolder+'gene2AncWeight.pkl','wb')
        pickle.dump(gene2AncWeight, output7)
        gene2AncWeight.clear()
        
    def Deltagene2Anc(self):
        '''
        '''
        fin = open(self.outfolder+"Gene.Excess_Deficiency.ANC.txt","wt")
        fin.write("GENE"+"\t"+"\t".join(self.anc_label)+"\t"+"P"+"\n")
        Delta_geneWeight = {}
        delta_file = open(self.outfolder+'Gene2snpDelta.pkl', 'rb')
        DeltaGene2snpAnc = pickle.load(delta_file)
        L = len(self.anc_label)-1
        for gene in DeltaGene2snpAnc:
            delta = []
            def_var = []
            AV,VA,dw = self.complex_list_comprehesive(DeltaGene2snpAnc[gene])
            delta_av = list(AV)
            delta_var = list(VA)
            for subset in iter.combinations(DeltaGene2snpAnc[gene][0][0],2):
                delta.append(abs(float(subset[0])-float(subset[1]))**2)
            for subset in iter.combinations(DeltaGene2snpAnc[gene][0][1],2):
                def_var.append(abs(subset[0]-subset[1]))
            value = sum(delta)/sum(def_var)
            W = [str(i) for i in delta_av]
            P = list(r.pchisq(value,df=L))[0]
            fin.write("\t".join([str(gene)]+W+[str(P)])+"\n")
            Delta_geneWeight[gene] = delta_av + [P]
        fin.close()
        AV = []
        delta = []
        w = []
        VA = []
        DeltaGene2snpAnc.clear()
        out = open(self.outfolder+'Gene.Excess_Deficiency.ANC.pkl','wb')
        pickle.dump(Delta_geneWeight, out)
        Delta_geneWeight.clear()
    
    def data_comb(self, snpFile, genoFile):
        '''
        Link genotype data to relative SNP and limit possible SNP to use for the program.
        Return: a dict, {snp:genotype}
        '''
        map_tmp = {}
        for line in fileinput.input(snpFile):
            data = line.split()
            map_tmp[fileinput.lineno()] = data
        self.tagget_dict = {}
        pos = {}
        for line in fileinput.input(genoFile):
            data = line.split()
            try:
                m_snp = map_tmp[fileinput.lineno()][0]
                self.tagget_dict[m_snp] = data[0] 
                pos[m_snp] = int(map_tmp[fileinput.lineno()][3])
            except:
                self.terminate()
                raise exceptions.SystemError('Failed to process, key not found ')
        del map_tmp    
        return self.tagget_dict
    
    def maxLDGene(self, networkFile, tagget_dict):
        '''
        Computes gene-gene-LD using the max LD between snps in genes
        networkFile: a string, the PPI network file
        tagget_dict: a dictionary, the tagged population genotype file
        '''
        logger.info('Start creating the network Gene-Gene-LD ...')
        rows_LD = {};myset=set()
        file5 = open(self.outfolder+"gene2geneLD.net","wt")
        for line in fileinput.input(networkFile):
            data2 = line.split()
            gene1 = data2[0]; gene2 = data2[1]            
            if gene1 in self.affy_map1 and gene2 in self.affy_map1 and gene1 != "---" and gene2 != "---":
                snp1 = self.affy_map1[gene1]
                snp2 = self.affy_map1[gene2]
                tmp = []
                for subset in iter.product(list(snp1),list(snp2)):
                    SNP_list = list(subset)
                    if SNP_list[0] in tagget_dict and SNP_list[1] in tagget_dict:
                        LD,cor = self.calc_rsq(tagget_dict[SNP_list[0]],tagget_dict[SNP_list[1]],SNP_list[0],SNP_list[1])
                        if LD != "NA":
                            tmp.append(LD)
                if len(tmp) !=0 and max(tmp)>0:
                    b = round(1 - max(tmp),7)
                    #file5.writelines(gene1+"\t"+gene2+"\t"+str(b)+"\n")    
                elif len(tmp) !=0 and max(tmp) < 0:
                    b = round(1 - abs(max(tmp)),7)
                file5.writelines(gene1+"\t"+gene2+"\t"+str(b)+"\n")
                myset.add(gene1); myset.add(gene2)
        logger.info('Genes included: %s' % str(len(list(myset))))
        file5.close()
        tagget_dict.clear()
        tmp = [] 
                       
    def ZscoreLDGene(self, networkFile, tagget_dict):
        '''
        Computes gene-gene-LD using the zscore of LDs between snps in genes
        networkFile: a string, the PPI network file
        tagget_dict: a dictionary, the tagged population genotype file
        '''
        logger.info('Start creating the network Gene-Gene-LD ...')
        rows_LD = {};myset=set()
        file5 = open(self.outfolder+"gene2geneLD.net","wt")
        for line in fileinput.input(networkFile):
            data2 = line.split()
            gene1 = data2[0]; gene2 = data2[1]            
            if gene1 in self.affy_map1 and gene2 in self.affy_map1 and gene1 != "---" and gene2 != "---":
                snp1 = self.affy_map1[gene1]
                snp2 = self.affy_map1[gene2]
                tmp = []
                for subset in iter.product(list(snp1),list(snp2)):
                    SNP_list = list(subset)
                    if SNP_list[0] in tagget_dict and SNP_list[1] in tagget_dict:
                        LD,cor = self.calc_rsq(tagget_dict[SNP_list[0]],tagget_dict[SNP_list[1]],SNP_list[0],SNP_list[1])
                        if LD != "NA":
                            tmp.append(LD)
                if len(tmp) !=0: 
                    a = np.mean(tmp)
                    if a > 0:
                        b = round(1 - a,7)
                        file5.writelines(gene1+"\t"+gene2+"\t"+str(b)+"\n")
                    else:
                        a = max(tmp)
                        if a > 0:
                            b = round(1 - a,7)
                            #file5.writelines(gene1+"\t"+gene2+"\t"+str(b)+"\n")
                        elif a < 0:
                            b = round(1 - abs(a),5)
                        file5.writelines(gene1+"\t"+gene2+"\t"+str(abs(b))+"\n")
                    myset.add(gene1); myset.add(gene2)
        logger.info('Genes included: %s' % str(len(list(myset))))
        file5.close()
        rows_LD.clear()    
        tagget_dict.clear()
        tmp=[]
    
    def closestLDGene(self, networkFile, tagget_dict):  
        '''
        Computes gene-gene-LD using the LD between the 2 closest snps in 2 given genes
        networkFile: a string, the PPI network file
        tagget_dict: a dictionary, the tagged population genotype file
        ''' 
        logger.info('Start creating the network Gene-Gene-LD ...')
        snp_snp = {}
        snp_snp1 = {}
        tmp = set()
        H_gene = []
        myset=set()
        pkl_file = open(self.outfolder+'affy_map3.pkl', 'rb')
        self.affy_map3 = pickle.load(pkl_file)
        i = 0;j = 0;g = 0
        file5 = open(self.outfolder+"gene2geneLD.net","wt")
        for line in fileinput.input(networkFile):
            data2 = line.split()
            gene1 = data2[0]; gene2 = data2[1]        
            if gene1 in self.affy_map3 and gene2 in self.affy_map3:
                try:
                    snp1 = self.affy_map3[gene1][0]
                    snp2 = self.affy_map3[gene2][0]
                    if snp1 in tagget_dict and snp2 in tagget_dict:
                        LD,cor = self.calc_rsq(tagget_dict[snp1],tagget_dict[snp2],snp1,snp2)
                        if LD != "NA" or LD !="na":
                            if LD > 0:
                                if abs(LD) > 1 :
                                    LD1 = 0.9
                                else:
                                    LD1 = round(abs(1 - abs(LD)),7) 
                                j+=1
                                file5.writelines(str(gene1)+"\t"+str(gene2)+"\t"+str(abs(LD1))+"\n")
                                myset.add(gene1); myset.add(gene2)
                                H_gene.append(gene1)
                                H_gene.append(gene2)
                            elif LD < 0 :
                                if abs(LD) > 1:
                                    LD1 = 0.9
                                else:
                                    LD1 = round(abs(1 - abs(LD)),7) 
                                    file5.writelines(str(gene1)+"\t"+str(gene2)+"\t"+str(abs(LD1))+"\n")
                                    H_gene.append(gene1)
                                    H_gene.append(gene2)
                                    myset.add(gene1); myset.add(gene2)
                except (RuntimeError, TypeError, NameError):
                    pass
            elif gene1 in self.affy_map1 and gene1 not in H_gene :
                SNP = self.affy_map1[gene1]
                temp =[]
                for subset in iter.product(SNP,SNP):
                    SNP_list = list(subset)
                    if SNP_list[0] != SNP_list[1] :
                        if SNP_list[0] in tagget_dict and SNP_list[1] in tagget_dict:
                            LD,cor = self.calc_rsq(tagget_dict[SNP_list[0]],tagget_dict[SNP_list[1]],SNP_list[0],SNP_list[1])
                            if LD != "NA" or LD !="na":
                                temp.append(LD)
                a = np.mean(temp)
                b = a/sqrt(2)
                if b != "NA" or b !="na":
                    try:
                        if abs(b) > 1: 
                            LD1 = 0.9
                            j+=1
                            file5.writelines(str(gene1)+"\t"+str(gene2)+"\t"+str(abs(LD1))+"\n")
                            myset.add(gene1); myset.add(gene2)
                            H_gene.append(gene1)
                        else:
                            if  b > 0:
                                LD1 = round(abs(1-abs(b)),7) 
                                j+=1
                                file5.writelines(str(gene1)+"\t"+str(gene2)+"\t"+str(abs(LD1))+"\n")
                                H_gene.append(gene1) 
                                myset.add(gene1); myset.add(gene2)   
                    except (RuntimeError, TypeError, NameError):
                        pass
            elif gene2 in self.affy_map1 and gene2 not in H_gene:
                H_gene.append(gene2)
                SNP = self.affy_map1[gene2]
                temp =[]
                try:
                    for subset in iter.product(SNP,SNP):
                        SNP_list = list(subset)
                        if SNP_list[0] != SNP_list[1] :
                            if SNP_list[0] in tagget_dict and SNP_list[1] in tagget_dict:
                                LD,cor = self.calc_rsq(tagget_dict[SNP_list[0]],tagget_dict[SNP_list[1]],SNP_list[0],SNP_list[1])
                                if LD != "NA":
                                    temp.append(LD)
                    a = np.mean(temp)
                    b = a/sqrt(2)
                    if b != "NA" or b !="na":
                        if abs(b) > 1: 
                            LD1 = 0.9
                            j+=1
                            file5.writelines(str(gene1)+"\t"+str(gene2)+"\t"+str(abs(LD1))+"\n")
                            myset.add(gene1); myset.add(gene2)
                        else:
                            if  b > 0:
                                LD1 = round(abs(1-abs(b)),7) 
                                j+=1
                                file5.writelines(str(gene1)+"\t"+str(gene2)+"\t"+str(abs(LD1))+"\n")
                                myset.add(gene1); myset.add(gene2)
                except (RuntimeError, TypeError, NameError):
                    pass
            else:
                i+=1
        logger.info('Missing LD SNP: %s' % str(i))
        logger.info('pair LD SNP included: %s' % str(j))
        logger.info('Genes included: %s' % str(len(list(myset))))
        file5.close()
        pkl_file.close()

    def calc_rsq(self, genotypes1, genotypes2, rsid1, rsid2):
        '''
        Computes the correlation between two genotype data.
        genotypes1,genotypes2: lists, lists of gentype
        rsid1,rsid2: strings, SNPs to cpmute the LD
        '''
        gen1,gen2 = self.remove_missing(genotypes1,genotypes2,rsid1,rsid2) 
        snp=[]
        if len(gen1)==0 and len(gen2)==0:
            snp.append(rsid2)
            snp.append(0.0)
            return "NA" # too much missing data to compute r
        elif len(gen1) != 0 and len(gen2) != 0:
            if len(gen1) ==len(gen2):    
                corr = self.get_weighted_r(gen1,gen2)
                y = 0.5*log((1+corr)/(1-corr))
                l = sqrt(len(gen1)-3)*y
                return abs(l),abs(corr)
            else:
                logger.info('genotypes do not match (exiting)')
                raise SystemExit

    def get_weighted_r(self, Y, Z):
        '''
        Correlation between two list of values by weighted scale
        '''
        mY, vY, mZ, vZ, cov = self.bivmom(Y,Z)
        den = sqrt(vY*vZ)
        m = len(Z)
        if float(den) > 0.0:
            r1 = float(1 + cov/sqrt(vY*vZ))
            r2 = float(1 - cov/sqrt(vY*vZ))
            if r2 > 0 and r1 >= 0:
                y = log(r1/float(r2))
                return y/float(sqrt(m-3))
            elif r2 < 0 and r1 <= 0:
                y = log(r1/float(r2))
                return y/float(sqrt(m-3))
            else:
                return 0.0
        else:
            return 0.0   
        
    def bivmom(self, vec0, vec1):
        '''
        Calculate means, variances, the covariance, from two data vectors.
        On entry, vec0 and vec1 should be vectors of numeric values and
        should have the same length.  Function returns m0, v0, m1, v1,
        cov, where m0 and m1 are the means of vec0 and vec1, v0 and v1 are
        the variances, and cov is the covariance.
        '''
        m0 = m1 = v0 = v1 = cov = 0
        for x, y in zip(vec0, vec1):
            m0 += x
            m1 += y
            v0 += x*x
            v1 += y*y
            cov += x*y
        n = len(vec0)
        assert n == len(vec1)
        n = float(n)
        m0 /= n
        m1 /= n
        v0 /= n
        v1 /= n
        cov /= n
        cov -= m0 * m1
        v0 -= m0 * m0
        v1 -= m1 * m1
        return m0, v0, m1, v1, cov 
    
    def remove_missing(self, genotypes1,genotypes2,rsid1,rsid2):
        '''
        Assessing missing genotypes
        '''
        geno = {"0":0,"1":1,"2":2,"?":9,"9":9,"-":9,"NA":9}
        if len(genotypes1) != len(genotypes2):
            logger.info('genotypes should be of the same length (exiting)')
            self.terminate()
            raise SystemExit
        else:
            gen1 = [geno[i] for i in genotypes1]
            gen2 = [geno[i] for i in genotypes2]
        return gen1,gen2

class ancGraph():
    '''
    Read the LD-weighted network, computes centraly measures, break down the PPI network to generate sub-networks 
    '''
    def __init__(self):
        return None

    def readGraph(self):
        '''
        Reading the LD-weighted network from the file gene2geneLD.net, generated in previous steps
        '''
        self.G = Graph()    
        fp = open(self.outfolder+'gene2geneLD.net')
        for line in fp:
            ligne = line.strip()
            if not ligne or ligne.startswith('#'): continue
            ligne = ligne.split()
            if float(ligne[-1]) <= 0.0:
                pass
            else:
                self.G.add_edge(ligne[0],ligne[1], weight=1.0/float(ligne[-1]))
        fp.close()
        return self.G
    
    def nodekLink(self):
        '''
        This function takes a network and returns protein degree distribution
        '''
        self.Degree = {}
        for node in self.G:
            self.Degree[node] = self.G.degree(node)
        NumberVertices_for_Nodes = self.Degree.values()
        Kvertices = list(set(sorted(NumberVertices_for_Nodes)))
        CountNodes = [NumberVertices_for_Nodes.count(n) for n in Kvertices]
        Proportions = [float(CountNodes[i])/sum(CountNodes) for i in xrange(len(CountNodes))]
        self.Gamma = leastsq(self.low_Distance,0,(np.array(Kvertices),np.array(Proportions)))    
        x = np.linspace(1,Kvertices[-1],1000)
        pylab.figure(1,dpi=150)    
        pylab.plot(Kvertices,Proportions,'ro', x, power(x,self.Gamma[0]),'k--',linewidth=2)
        pylab.legend(("Connectivity-Distr.: "+r'$\mathcal{P}\left(k\right)$',"Power-Law: "+r'$\mathcal{P}\left(k\right)\simeq k^{%.2f}$'%(self.Gamma[0])))
        pylab.ylim(-0.001,max(Proportions)+0.001)
        pylab.xlim(0.0,Kvertices[-1])
        pylab.grid(True)
        pylab.xlabel("Detected Protein Degree: "+r'$k$'); pylab.ylabel("Connections Frequency: "+r'$\mathcal{P}\left(k\right)$')
        pylab.savefig(self.outfolder+'PowerLawData.png',dpi=150)
    
    def plotPathDistr(self):
        '''
        Compute shortest paths between all nodes in a weighted graph.
        This function takes a network path length and returns path distribution distribution and
        mean shortest paths needed for other computational purposes
        '''
        paths={}
        LPath = []
        for n in self.G:
            Length = single_source_dijkstra_path(self.G,n)
            for b in Length:
                if n==b: continue
                LPath.append(len(Length[b])-1)
        self.ShortPathMean = np.mean(LPath)
        import matplotlib.pyplot as plt
        plt.figure(2,dpi=150)
        n, bins, patches = plt.hist(LPath, normed=1, facecolor='0.8', alpha= 1.0)
        t = np.linspace(0,max(LPath))
        y = np.exp(-(t-np.mean(LPath))**2/(2*np.std(LPath)))/(np.std(LPath)*np.sqrt(2*np.pi))
        plt.plot(t,y,'r--',linewidth=2)
        plt.grid(True)
        plt.xlabel("Path-Length: "+r'$\ell$')
        plt.ylabel("Frequency: "+r'$\mathcal{P}\left(\ell\right)$')
        plt.savefig(self.outfolder+'PathDistr.png',dpi=150)
    
    def findallHubs(self):
        '''
        This returns all hubs of a given network and all connected components
        '''
        Hubs = []
        C = connected_component_subgraphs(self.G) 
        for cliques in C:
            if cliques.order() <= 2: continue
            for gene in cliques:
                Ctemp = cliques.copy()
                Ctemp.remove_node(gene)
                if number_connected_components(Ctemp) > 1: 
                    Hubs.append(gene) #Gene under consideration is hub
                del Ctemp
        output4 = open(self.outfolder+'Hubs.pkl','wb')
        pickle.dump(Hubs, output4)
        return Hubs
    
    def low_Distance(self, parameters,x,values):
        '''
        This function receives one parameter of power-low gamma model into 
        parameters, an array with the x and an array with the corresponding 
        values for each of the x. 
        our model is x**(-gamma)
        '''
        gamma = parameters
        errors = values - x**gamma
        return errors
    
    def subgraphFinding(self):
        '''
        This function computes the centraly measures, break down the PPI network to generate sub-networks.
        '''
        logger.info('Computing protein degree distribution, it may take time ...') 
        logger.info('1. Now plotting Power-Law ...')
        self.nodekLink()
        logger.info('2. Computing length distribution, it may take time ...')
        self.plotPathDistr()
        logger.info("Shortest path mean "+str(self.ShortPathMean))
        logger.info('3. Searching for network structure hubs and it may take time ...')
        Hubs = self.findallHubs()
        logger.info('4. Computing Node betweenness scores...')
        self.Betw = betweenness_centrality(self.G)
        output5 = open(self.outfolder+'Betw.pkl','wb')
        pickle.dump(self.Betw, output5)
        logger.info('5. Computing Node closeness scores')
        self.Clos = closeness_centrality(self.G)
        output5 = open(self.outfolder+'Clos.pkl','wb')
        pickle.dump(self.Clos, output5)
        logger.info('6. Computing Node eigenvector scores')
        EigenScoreWorks = 0
        A = adj_matrix(self.G)
        try:
            a,b,c = svd(A) # We will need either a[:,0] or c[0]
            if any(c[0] < 0.0): vecmax = -c[0]
            else: vecmax = c[0]
            EigenScoreWorks = 1
        except:
            logger.info('Singularity problem occurs while computing Eigen-vector scores ...')
        self.BetOf = self.ShortPathMean * self.G.order()
        if EigenScoreWorks: EigOf = np.mean(vecmax)
        self.ClosOf = 1.0/self.ShortPathMean
        self.DegOf = 2.0*self.G.size()/self.G.order()
        Node = self.G.nodes()
        self.BetweenNodes = set([node for node in self.G if self.Betw[node] >= self.BetOf])
        self.CloseNodes = set([node for node in self.G if self.Clos[node] >= self.ClosOf])
        if EigenScoreWorks: 
            self.EigenNodes = set([Node[i] for i in xrange(len(Node)) if vecmax[i] >= self.EigOf])
        else: 
            self.EigenNodes = set([node for node in self.G if len(self.G[node]) >= self.DegOf])
        output4 = open(self.outfolder+'ShortPathMean.pkl','wb')
        pickle.dump(self.ShortPathMean, output4)
        output1 = open(self.outfolder+'EigenNodes.pkl','wb')
        pickle.dump(self.EigenNodes, output1)
        output2 = open(self.outfolder+'CloseNodes.pkl','wb')
        pickle.dump(self.CloseNodes, output2)
        output3 = open(self.outfolder+'BetweenNodes.pkl','wb')
        pickle.dump(self.BetweenNodes, output3)
        logger.info('7. Running the last step: searching for subgraphs ...')
        if  len(Hubs) > 1 :
            if  len(self.BetweenNodes & self.CloseNodes & self.EigenNodes)> 1 :
                self.CenterGenes = set(Hubs) & self.BetweenNodes & self.CloseNodes & self.EigenNodes
            elif len(self.BetweenNodes & self.CloseNodes)>1 and len(self.EigenNodes)==0 :
                self.CenterGenes = set(Hubs) & self.BetweenNodes & self.CloseNodes
            elif len(self.BetweenNodes & self.EigenNodes)> 1 and len(self.CloseNodes)==0:
                self.CenterGenes = set(Hubs) & self.BetweenNodes & self.EigenNodes
            elif len(self.EigenNodes & self.CloseNodes)>1 and len(self.BetweenNodes)==0:
                self.CenterGenes = set(Hubs) & self.CloseNodes & self.EigenNodes
            else:
                self.CenterGenes = set(Hubs) 
        elif len(Hubs) == 1:
            logger.info('\nFatal Error: No CenterGene found!!!\nancGWAS Could not continue\n')
            self.terminate()
            sys.exit(1)
        else:
            if  len(self.BetweenNodes & self.CloseNodes & self.EigenNodes) > 1 :
                self.CenterGenes = self.BetweenNodes & self.CloseNodes & self.EigenNodes
            elif len(self.BetweenNodes & self.CloseNodes)>1 and len(self.EigenNodes)==0 :
                self.CenterGenes = self.BetweenNodes & self.CloseNodes    
            elif len(self.BetweenNodes & self.EigenNodes)> 1 and len(self.CloseNodes)==0:
                self.CenterGenes = self.BetweenNodes & self.EigenNodes
            elif len(self.EigenNodes & self.CloseNodes)>1 and len(self.BetweenNodes)==0:
                self.CenterGenes = self.CloseNodes & self.EigenNodes
        self.Clouds = {}; n = int(self.step)
        for central in self.CenterGenes:
            self.Clouds[central] = set([central])
            Temp = set(self.G.neighbors(central))
            i = 0
            while i < n:
                Current = Temp-self.Clouds[central]; Temp.clear()
                for gene in Current:
                    self.Clouds[central].add(gene)
                    Temp |= set(self.G.neighbors(gene))
                i += 1
        output = open(self.outfolder+'CloudsStep.pkl','wb')
        pickle.dump(self.Clouds, output)

class ancScoring():
    '''
    This module computes the score of each sub-netowk to detect the ones enriched in p-value and particular ancestries
    '''
    
    def __init__(self):
        return None
    
    def Scoring(self, opt):
        '''
        Scoring sub-networks
        '''
        self.opt = opt
        self.outfolder = opt[-1]
        gene2LD = {}
        for line in fileinput.input(self.outfolder+"gene2geneLD.net"):
            data = line.split()
            gene2LD[data[0]+":"+data[1]] = float(data[2])
        logger.info('Reading saved subnetworks ...')
        pkl_file = open(self.outfolder+"CloudsStep.pkl", 'rb')
        Gene2mod = pickle.load(pkl_file)
        logger.info('Loading gene-pvalue ...')
        gene_p = open(self.outfolder+"gene2Pv.pkl", 'rb')
        gene2pval = pickle.load(gene_p)
        logger.info('Statistic model for scoring subnetworks and writing the final result...')
        suM,suM1 = self.module_scoring(Gene2mod,gene2pval,opt[2:])
        if self.anced == True:
            logger.info('Loading excess/deficiency stored statistics ...')
            pkl_prop = open(self.outfolder+"gene2AncWeight.pkl", 'rb')
            anc = pickle.load(pkl_prop)
            pkl_delta = open(self.outfolder+"Gene.Excess_Deficiency.ANC.pkl","rb")
            delta = pickle.load(pkl_delta)
            logger.info('Processing subnetwork excesss/deficiency statistics ...')
            self.module_delta(suM1,anc,delta,opt)
            self.subnetwork_anc(suM, gene2LD, gene2pval, opt, anc, delta)
        else:
            self.subnetwork_anc(suM, gene2LD, gene2pval, opt)
    
    def module_scoring(self, Gene2module, gene2pval, opt):
        '''
        Computes the score of each module and compute the 95% confidence interval
        '''
        tmp = [] 
        score_mod = [] 
        for mod in Gene2module:
            tmp = []
            if mod in gene2pval:
                tmp.append(gene2pval[mod])
                for gene in Gene2module[mod]:
                    if gene not in gene2pval:
                        continue
                    else:
                        tmp.append(float(gene2pval[gene])) 
                if opt[1]:
                    pval,chisquare = self.fisher_method(tmp)
                    score_mod.append([mod,chisquare])
                elif opt[0]:
                    pval,chisquare = stouffer_zscore(tmp)
                    score_mod.append([mod,chisquare ])
        pval_norm = self.normalization_p(score_mod) # Accounting if the obtained chi2 is not by chance and compute the 95% CI for each subnetwork
        rows = set();rows1 = [] 
        fin = open(self.outfolder+"module.score","wt")
        fin.writelines("95% CI"+"\t"+"nScore"+"\t"+"Zscore"+"\t"+"Hubs.Genes"+"\t"+"Lists.Genes"+"\n")
        for des in pval_norm[-int(opt[-2]):]:
            if len(des) != 0:
                if des[1] in Gene2module:
                    tmp = list([des[3]])+list([des[0],des[2],des[1]])+list(Gene2module[des[1]])
                    fin.write("\t".join([str(de) for de in list(tmp)])+"\n")
                    rows |= set(tmp[3:]) 
                    tp = [des[1]]+list(Gene2module[des[1]])
                    rows1.append(tp)
            else:
                continue
        fin = open(self.outfolder+"module.all.score","wt")
        fin.writelines("95% CI"+"\t"+"nScore"+"\t"+"Zscore"+"\t"+"Hubs.Genes"+"\t"+"Lists.Genes"+"\n")
        for des in pval_norm:
            if len(des) != 0:
                if des[1] in Gene2module:
                    tmp = list([des[3]])+list([des[0],des[2],des[1]])+list(Gene2module[des[1]])
                    fin.write("\t".join([str(de) for de in list(tmp)])+"\n")
            else:
                continue
        
        fin.close()
        return rows,rows1
    
    def fisher_method(self, x):
        '''
        Returns the p-value of getting chi2 from a chi-squared distribution.
            chi2: observed chi-squared statistic
            df: degrees of freedom
        '''
        tp = []
        for i in range(len(x)):
            if x[i] <= 0:
                continue
            else:
                tp.append(log(x[i]))
        value = -2 * sum(tp) 
        dof = 2*len(tp)
        tp = []
        t = list(r.pchisq(value,df=dof)) 
        p = 1.0-t[0]
        return p,value
    
    def normalization_p(self, X):
        '''
        Returns the normalized p-value overall all the modules.
        '''
        y2 = []
        Y = [] 
        for des in X:
            Y.append([des[0],des[1]])
            y2.append(des[1])
        X = []
        n = len(y2)
        for de in y2:
            idx=y2.index(de)
            element = de
            y2.remove(de)
            a1 = np.mean(y2)
            a2 = np.var(y2)
            s = sqrt(a2)
            a3 = element-a1
            y2.insert(idx,element)
            Z = a3/float(a2)
            Y[idx].insert(0,Z)
            se = abs(Z/((r.qnorm(0.975)[0]*s)/sqrt(n)))
            left=Z-se
            right=Z+se
            ci="("+str(left)+","+" "+str(right)+")"
            Y[idx].append(ci)
        y2 = []
        Y.sort() 
        return Y

    def module_delta(self, su, anc, delta, opt):
        tmp = [] 
        score_anc = [] 
        score_delta = [] 
        score_mod = [] 
        fix = open(self.outfolder+"subnet.ancestry.out","wt")
        fix.write("\t".join(self.anc_label)+"\t"+"Lists.Genes"+"\n")
        fin = open(self.outfolder+"subnet.delta.ancestry.out","wt")
        fin.write("\t".join(["95%CI","aP","emP"]+self.anc_label)+"\t"+"Lists.Genes"+"\n")
        for mod in su: 
            tmp_anc = [] 
            tmp_delta = [] 
            tmp= [] 
            for gene in mod:
                if gene in anc:
                    tmp_anc.append(anc[gene])
                if gene in delta:
                    tmp.append(np.log(delta[gene][-1:][0]))
                    tmp_delta.append(np.array(delta[gene][:-1]))
            ANC = self.inverse_dict(tmp_anc)[0]
            DELTA = self.inverse_dict(tmp_delta)[0]
            T = -2 * sum(np.array(tmp))/(sqrt(len(tmp)))
            CHI = list(r.pchisq(T,len(tmp)))[0]
            score_anc.append([mod,ANC])
            score_delta.append([mod,DELTA])
            score_mod.append([mod,CHI])
        pval_norm = self.normalization_p(score_mod)
        for a in score_anc:
            fix.write("\t".join([str(de) for de in list(a[1])]+list(a[0]))+"\n")
        for des in pval_norm:
            if len(des) != 0:
                idx = pval_norm.index(des)   
                tmp = list([des[3]])+list([des[0],des[2]])+list(score_delta[idx][1])+list(des[1])
                fin.write("\t".join([str(de) for de in list(tmp)])+"\n")
            else:
                continue
        fin.close()
        fix.close()
    
    def inverse_dict(self, list1):
        '''
        Return the inverse of a dictionary vlaue:key
        '''
        list2 = []
        list3 = []
        for i in range(len(list1[0])):
            tmp =[]
            for j in range(len(list1)):
                tmp.append(list1[j][i])
            list2.append(np.mean(tmp))
            list3.append(np.var(tmp))
        return list2,list3
    
    def subnetwork_anc(self, sub_module, LD, Pva, opt, anc = '', delta_anc = ''):
        '''
        This function take a biological network and a significant modules to create a subnetwork
        '''
        ld_cutoff = opt[1]
        fi = open(self.outfolder+"sub.network.out","wt")
        fi.writelines("GeneA"+"\t"+"GeneB"+"\t"+"P"+"\t"+"\n")
        tmp = []
        genes = []
        anc_av = []
        if self.anced == True:
            for i in range(len(self.anc_label)):
                tmp.append([])
                anc_av.append([self.anc_label[i]])
        PPI = set()
        for line in fileinput.input(opt[0]):
            data = line.split()
            if data[0] in sub_module and data[1] in sub_module:
                if data[0]+":"+data[1] in LD:
                    GeneLD = data[0]+":"+data[1]
                    if LD[data[0]+":"+data[1]] >= ld_cutoff :
                        if data[1]+":"+data[0] not in PPI:        
                            PPI.add(GeneLD+":"+str(LD[data[0]+":"+data[1]]))
                elif data[1]+":"+data[0] in LD:
                    GeneLD = data[1]+":"+data[0]
                    if LD[data[1]+":"+data[0]]>=ld_cutoff:
                        if data[0]+":"+data[1] not in PPI:
                            PPI.add(GeneLD+":"+str(LD[data[1]+":"+data[0]]))
        for des in list(PPI):
            dim = des.split(":")
            fi.write("\t".join(dim)+"\n")
            genes +=dim
        fi.close()
        if self.anced == True:
            for de in genes:
                if de in anc:
                    for j in range(len(anc[de])):
                        tmp[j].append(float(anc[de][j]))
            if len(anc_av) != 0:
                for k in range(len(anc_av)):
                    anc_av[k].append(np.mean(tmp[k]))
                logger.info('Summary Ancestry Network')
                for i  in range(len(self.anc_label)):
                        logger.info("          "+self.anc_label[i]+": "+str(anc_av[i][1]))
        if opt[-3]:
            process = subprocess.Popen(['Rscript', os.getcwd()+'/Plotting.R', self.outfolder])
            process.wait()

class ancGWAS(ancInit, ancUtils, ancmap, ancLD, ancGraph, ancScoring):
    '''
    Performs the overrall steps of ancGWAS as described in the ancGWAS method (Chimusa et al 2013)
    '''
    def anc(self):
        '''
        Running ancGWAS
        '''
        logger.info('Starting at time:%s' % str(datetime.datetime.today()))
        self.res = ancGWAS(argv)      
        logger.info("Loading parameters from %s ..."%os.path.abspath(argv))
        logger.info("Options in effect:")
        for param in sorted(self.res.Params):
            logger.info('             '+str(self.res.Params[param][0])+': '+str(self.res.Params[param][1]))
        print ''
        self.res.readGWAS()
        self.res.readAffy()
        gpv = {'Simes':[self.res.simes(),'Method used is Simes.'], 'Smallest':[self.res.smallest(), 'Method used is smallest model.'], 'Fisher':[self.res.fisher(), 'Method used is Fisher.'], 'Gwbon':[self.res.FDR(), 'Method used is gene-wise FDR.']}                
        logger.info('Computing Gene pvalues ...')
        logger.info('Method used is %s.'%self.Gene_pv)
        gpv[self.Gene_pv.capitalize()][0] 
        try:
            if self.anced == True:
                self.res.delta()
                logger.info('Computing Average proportions of ancestries for each Gene...')
                self.res.Gene2Anc()
                self.res.Weightedgene2Anc()
                logger.info('Computing excess/deficiency ancestry statistic for each gene...')
                self.res.Deltagene2Anc()
                logger.info('Gene Ancestry proportion and Excess/Deficiency ancestry test are done ...')
        except AttributeError:
            sys.stderr.write('Bad option. Missing ancestry proportions. Check your parameters')
            sys.exit(1)
        tagget_dict = self.res.data_comb(self.snpFile, self.genoFile)
        logger.info('Writing into a file the Gene-Gene LD ...')
        logger.info("Using the %s method for both subnetwork significance and LD ..."%self.Gene_LD)
        if self.Gene_LD.lower() == 'zscore':
            self.res.ZscoreLDGene(self.networkFile, tagget_dict)
        elif self.Gene_LD.lower() == 'closest':
            self.res.closestLDGene(self.networkFile, tagget_dict)
        elif self.Gene_LD.lower() == 'maxscore':
            self.res.maxLDGene(self.networkFile, tagget_dict)
        currenttime1 = time.asctime(time.localtime())
        logger.info('Finish weighting the network at %s' % str(currenttime1))
        self.res.readGraph()
        currenttime = time.asctime(time.localtime())
        logger.info('Start searching at %s' % str(currenttime))
        self.res.subgraphFinding() 
        currenttime1 = time.asctime(time.localtime())
        logger.info('Scoring generated subnetwork at %s' % str(currenttime1))
        self.res.Scoring([self.networkFile, self.LDCutoff, self.stouffer, self.wzscore, self.Pcutoff, self.step, self.lampld, self.tklc, self.TOPScore, self.outfolder])
        self.res.terminate()
        logger.info("Finish at time:%s"%str(datetime.datetime.today()))
        
if __name__ == '__main__':
    try:
        argv = sys.argv[1]
    except IndexError:
        argv = ''
    finally:
        run = ancGWAS(argv)
        run.anc()
