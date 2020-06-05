#!/opt/python/3.6.3/bin/python3

#==============================================================================                                                                                                       
__author__ = 'Rutendo F. Sigauke'
__credits__ = ['Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Rutendo F. Sigauke'
__email__ = 'Rutendo.Sigauke@colorado.edu' #or 'Rutendo.Sigauke@CUAnschutz.edu'                                                                                                       
#==============================================================================                                                                                                       
import os
import sys
import time
import numpy as np


def load_bed(path_to_bed):
    '''Loads a bedfile and returns a list of
    coordinates [[chrm, start, stop, ...]]
                           
    Parameters                   
    ----------                
    path_to_bed : str                 
        path to the input bed file                                                                                              
        
    Returns            
    -------  
    bed_list : list
        list of coordinates from a bed file
    '''

    bed_list = []
    with open(path_to_bed) as bed:
        for line in bed:
            line = line.strip()
            coord = line.split('\t')
            bed_list.append(coord)
            
    return bed_list


def load_fimo(fimo_file):
    '''load fimo tsv file and return a list
    of all results.
    
    Parameters                   
    ---------- 
    fimo_file : str
        paths to fimo .tsv file
    
    Returns            
    ------- 
    fimo_results : list
        entries of Fimo TF motif hits
    
    '''
    print('---------------------------------------------------------------------')
    print("--------------------- STORING FIMO RESULTS FILE ---------------------")
    fimo_results = []
    with open(fimo_file) as motifs:
        next(motifs)
        for row in motifs:
            row = row.strip()
            if not row:  
                continue
            if row.startswith("#"):  
                continue
            results = row.split('\t')
            fimo_results.append(results)
            
    print('-> Number of Fimo Calls = ' + str(len(fimo_results)))       
    
    return fimo_results


def filter_fimo(fimo_res, alpha = 1e-04):

    '''Filter fimo hits based on significance
    
    Parameters                        
    ----------                        
    fimo_res : str             
        path to fimo tsv output                                                              
                                    
    alpha : float                   
        significance cut-off (default 1e-04)
                                                                                                                                
    Returns   
    -------
    significant_results : list of lists
        fimo hits with signicicant matches at states threshold
             
    '''
    
    print("--------- FILTERING MOTIF HITS AT AN ALPHA LEVEL OF {} ----------".format(alpha))
    significant_results = []
    for res in fimo_res:
        if float(res[8]) <= float(alpha):
            significant_results.append(res)
            
    return significant_results


def unique_sig_fimo(sig_results):
    '''function takes list of fimo results and returns
    a list of unique hits per region (input from filter fimo)
    
    Parameters
    ----------
    sig_results : list of lists
        fimo hits in list
    
    Returns
    -------
    unique_fimo_hits_list : list of list
        unique regions with at-least 1 fimo hits
    
    '''
    
    print("------------------------- UNIQUE FIMO HITS --------------------------")
    unique_fimo_hits = {}
    unique_fimo_hits_list = []

    for sig in sig_results:
        hit_id = sig[2]
        chrom = str.split(str(sig[2]),':')[0]
        start = int(str.split((str.split(str(sig[2]),':'))[1],'-')[0]) + 150
        stop = int(str.split((str.split(str(sig[2]),':'))[1],'-')[1]) - 150
        significance = sig[8]
        if hit_id not in unique_fimo_hits:
            unique_fimo_hits[hit_id] = [[chrom, start, stop, significance]]
            unique_fimo_hits_list.append([chrom, start, stop, significance])
            
    return unique_fimo_hits_list


def fetch_merged_regions(merged_region):
    '''Takes merged bed files and returns the center of the called
    width of region, center if region as well as coordinates.
    
    Parameters 
    ----------
    merged_region : str
        path to merged bed file
    
    Returns
    -------
    width : list
        widths for bed regions 
        
    '''
    
    width = []
    mu = []
    region = []
    with open(merged_region) as regions:
        for row in regions:
            row = row.strip()
            results = row.split('\t')
            w_ = int(results[2]) - int(results[1])
            mu_ = int(w_/2)
            id_ = str(results[0])+':'+ str(results[1])+'-'+str(results[2])
            width.append(w_)
            mu.append(mu_)
            region.append(results)
            
    return width #, mu, region

def sig_summit_peaks(summit_peak, unique_fimo_list):
    '''Finds peak regions with motif hits
    
    Parameters
    ----------
    summit_peak : str
        path to summit peak bed files
        
    unique_fimo_list : list of list
        fimo TF motif hits list 
        
    Returns
    -------
    peaks_with_motifs : list of lists
        peak regions with TF motifs
         
    '''

    print("-------------- SIGNIFICANT PEAKS WITH FIMO MOTIF HIT ----------------")

    macs_results = []
    peaks_with_motifs = []
    
    with open(summit_peak) as peaks:
        for row in peaks:
            row = row.strip()
            macs_peaks = row.split('\t')
            macs_results.append(macs_peaks)

    for sig in unique_fimo_list:
        for macs in macs_results:
            if str(sig[0]) == str(macs[0]) and str(sig[1]) == str(macs[1]) and str(sig[2]) == str(macs[2]):
                peaks_with_motifs.append([sig[0], sig[1], sig[2],macs[3]])
                
    print('-> Macs peaks called = ' + str(len(macs_results)))
    
    return peaks_with_motifs

def final_narrow_peaks(narrow_peaks, peaks_with_motifs, outdir, outfile):
    '''Returns the original _narrowPeaks with motif peaks, 
    since motifs were scanned using a summit regions
    
    Parameters
    ----------
    narrow_peaks : str
        path to macs _narrowPeaks peak calling files
        
    peaks_with_motifs : list of lists
        peak regions with significant motif hits
        
    outdir : str
        directory to store new bed files
        
    outfile : str
        name of the bed file
        
    Returns
    -------
    writes a bed file in specified location
    
    '''

    print("---------- SAVE NARROW PEAKS WITH SIGNIFICANT MOTIF HITS ------------")
    
    sig_fimo_narrow = []
    
    with open(narrow_peaks) as peaks:
        for row in peaks:
            row = row.strip()
            macs_peaks = row.split('\t')
            for hits_id in peaks_with_motifs:
                if hits_id[3] == macs_peaks[3]:
                    sig_fimo_narrow.append([macs_peaks[0], macs_peaks[1], macs_peaks[2], macs_peaks[3]])
                    
    print('-> Peaks with significant hits = ' + str(len(sig_fimo_narrow)))
    
    np.savetxt(outdir + outfile,
               np.array(sig_fimo_narrow), delimiter="\t",fmt='%s')
    print('------------------- DONE RUNNING {} -------------------'.format(outfile))
    print('---------------------------------------------------------------------')

            