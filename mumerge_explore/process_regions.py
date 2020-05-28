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
