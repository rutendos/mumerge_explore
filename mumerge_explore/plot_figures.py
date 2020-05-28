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
import numpy as np
import matplotlib.pyplot as plt

def compare_widths(merge, sample, merge_label, sample_label, plot_title):
    
    plt.figure(figsize=(18,7))
    gs1 = plt.GridSpec(1,1)

    ax1 = plt.subplot(gs1[0])

    #bins=np.histogram(np.hstack((merge)), bins=60)[1]
    bins=60
    ax1.hist(merge, bins, normed=1, 
             facecolor='#d7191c',edgecolor='black', alpha=0.65, 
             linewidth=1.2, label=merge_label)
    ax1.hist(sample, bins, normed=1,
             facecolor='#fdae61',edgecolor='black', alpha=0.60,
             linewidth=1.2, label=sample_label)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    plt.xlabel('Region Width (bp)',fontsize=18,fontweight='bold')
    plt.ylabel('Probability',fontsize=18,fontweight='bold')
    plt.title(plot_title,fontsize=25,fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1),fontsize = 15, loc=2, borderaxespad=0.)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.show()
    