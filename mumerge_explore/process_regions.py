#!/opt/python/3.6.3/bin/python3

#==============================================================================                                                                                                       
__author__ = 'Rutendo F. Sigauke'
__credits__ = ['Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Rutendo F. Sigauke'
__email__ = 'Rutendo.Sigauke@colorado.edu' #or 'Rutendo.Sigauke@CUAnschutz.edu'                                                                                                       
#============================================================================== 


def merged_center(in_bed):
    
    bed_newmu = []
    
    with open(in_bed) as regions:
        for row in regions:
            row = row.strip()
            coord = row.split('\t')
            mid = round((int(coord[2]) - int(coord[1]))/2)
            center = int(coord[1]) + mid
            bed_newmu.append([coord[0], coord[1], coord[2], center])

    return bed_newmu

def get_center_coord(infile):
    
    center_coord = []
    with open(infile) as coord:
        for row in coord:
            row = row.strip()
            coords = row.split('\t')
            center = coords[1]
            center_coord.append(center)
        
    return center_coord

def get_center_coord_list(infile):
    
    center_coord = []
    
    for coord in infile:
        center = coord[-1]
        center_coord.append(center)
        
    return center_coord

def overlapping_regions(bed_a, bed_b):
    beda_overlapping = []
    
    for i in bed_a:
        for j in bed_b:
            if str(i[0]) == str(j[0]):
                check_overlap = (int(i[1]) - int(j[2])) * (int(i[2]) - int(j[1]))
                if check_overlap < 0:
                    beda_overlapping.append(i)
                    
    return beda_overlapping

def mean_diff(average,center):
    mean_diff_list = []
    for ave, cen in zip(average,center):
        cen_diff = int(ave) - int(cen)
        mean_diff_list.append(cen_diff)
        
    return mean_diff_list


def unique_sig_fimo_hits(sig_results):
    
    unique_fimo_hits = {}

    for sig in sig_results:
        hit_id = sig[2]
        chrom = str.split(str(sig[2]),':')[0]
        score = sig[6]
        significance = sig[8]

        ##start and stop of motif hit
        motif_start = int(str.split((str.split(str(sig[2]),':'))[1],'-')[0]) + int(sig[3])
        motif_stop = int(str.split((str.split(str(sig[2]),':'))[1],'-')[0]) + int(sig[4])

        ##start and stop of region with peak call
        region_start = int(str.split((str.split(str(sig[2]),':'))[1],'-')[0]) 
        region_stop = int(str.split((str.split(str(sig[2]),':'))[1],'-')[1]) 

        ##select most significant motif hit
        if hit_id not in unique_fimo_hits:
            unique_fimo_hits[hit_id] = [chrom, region_start, region_stop, 
                                        significance, motif_start, motif_stop, score]
        elif hit_id in unique_fimo_hits:
            if float(unique_fimo_hits[hit_id][3]) > float(significance): 
                unique_fimo_hits.update({hit_id : [chrom, region_start, region_stop, 
                                                   significance, motif_start, motif_stop, score]})
                
    return unique_fimo_hits


def summit_motif_dist(narrow, summits, unique_fimo_dict):
    motif_peak_distances = []
    motif_significance = []
    
    for macs, summit in zip(narrow, summits):
        for key, value in unique_fimo_dict.items():
            if str(value[0]) == str(macs[0]) and str(value[1]) == str(macs[1]) and str(value[2]) == str(macs[2]):
                motif_center = (int(value[5]-int(value[4]))+1)/2
                #distance = int(summit[1]) - int(value[4])
                distance = int(summit[1]) - (int(value[4]) + int(motif_center))
                motif_peak_distances.append(distance)
                motif_significance.append(float(value[3]))

    print('Sequences with motif hits => ' + str(len(motif_peak_distances)))
    return motif_peak_distances, motif_significance


def merged_center_bed(merged_narrow):
    
    summit_coords = []
    
    for regions in merged_narrow:
        chrm = regions[0]
        start = int(regions[3])
        stop = int(start + 1)
        summit_coords.append([chrm, start, stop])
        
    return summit_coords