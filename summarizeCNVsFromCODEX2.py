import sys
import math

from bokeh.core.properties import value
from bokeh.io import show, output_file
from bokeh.models import ColumnDataSource, FactorRange, Range1d, Legend, LegendItem, axes
from bokeh.plotting import figure
from bokeh.io import export_png, export_svgs
import pandas as pd

DEBUG = True

infile = sys.argv[1]  # CODEX/SS4_results/frac_files.txt
workingDir = sys.argv[2] # CODEX/SS4_results/
dgvfile = sys.argv[3] 
kit = sys.argv[4]
gencodefile = sys.argv[5]
sampleSourceFile = sys.argv[6]
sampNameFile = sys.argv[7]

chr_size = dict()
chr_size['1'] = 248956422
chr_size['2'] = 242193529
chr_size['3'] = 198295559   
chr_size['4'] = 190214555 
chr_size['5'] = 181538259 
chr_size['6'] = 170805979 
chr_size['7'] = 159345973    
chr_size['8'] = 145138636 
chr_size['9'] = 138394717 
chr_size['10'] = 133797422  
chr_size['11'] = 135086622 
chr_size['12'] = 133275309 
chr_size['13'] = 114364328
chr_size['14'] = 107043718
chr_size['15'] = 101991189 
chr_size['16'] = 90338345
chr_size['17'] = 83257441 
chr_size['18'] = 80373285
chr_size['19'] = 58617616 
chr_size['20'] = 64444167 
chr_size['21'] = 46709983
chr_size['22'] = 50818468
chr_size['Y'] = 57227415 
chr_size['X'] = 156040895 

GENCODE_overlap_threshold = 0.5

def buildGENCODEDict(ff):
    fr = open(ff, "r")
    header = fr.readline()
    gencodeStartEndMapping = dict()
    gencodeStartAnnotationMapping = dict() 
    for line in fr:
        tokens = line.strip().split("\t")
        if tokens[0][0:3] == "chr" and ((tokens[2] == "gene") or (tokens[2] == "exon")):        
            chromosome = tokens[0].replace("chr", "")
            start = int(tokens[3])
            end = int(tokens[4])
            if chromosome not in gencodeStartEndMapping:
                gencodeStartEndMapping[chromosome] = dict()
                gencodeStartEndMapping[chromosome][start] = end
            else:
                if start not in gencodeStartEndMapping[chromosome]:
                    gencodeStartEndMapping[chromosome][start] = end

            annotation = tokens[8]
            idx2 = tokens[8].find(";")
            ensid = annotation[3:idx2]
            idx1 = annotation.find("gene_name=")
            idx2 = annotation.find(";", idx1)
            genename = annotation[idx1+10:idx2]
            #print("genename = ", genename)
            key = chromosome + "_" + str(start)
            if key not in gencodeStartAnnotationMapping:
                gencodeStartAnnotationMapping[key] = genename + "|" + ensid 
                #gencodeStartAnnotationMapping[key] = genename + "|" + ensid + "|" + chromosome + "_" + str(start) +"_" + str(end)
        
    #print("gencodeStartEndMapping = ", gencodeStartEndMapping)
    #print("gencodeStartAnnotationMapping = ", gencodeStartAnnotationMapping)
    return gencodeStartEndMapping, gencodeStartAnnotationMapping

def buildDGVDict(ff):
    fr = open(ff, "r")
    header = fr.readline()
    dgv = dict()
    for line in fr:
        tokens = line.strip().split("\t")
        chromosome = tokens[1]
        start = int(tokens[2])
        end = int(tokens[3])
        variant = tokens[4]
        variantSubType = tokens[5]
        variantType = variant + " " + variantSubType
        reference = tokens[6]
        if chromosome not in dgv:
            dgv[chromosome] = dict()
            dgv[chromosome][start] = end
            #dgv[chromosome][start] = dict()
            #dgv[chromosome][start][end] = list()
            #dgv[chromosome][start][end].append(variantType)
            #dgv[chromosome][start][end].append(reference)
        else:
            if start not in dgv[chromosome]:
                dgv[chromosome][start] = end
                #dgv[chromosome][start][end] = list()
                #dgv[chromosome][start][end].append(variantType)
                #dgv[chromosome][start][end].append(reference)
            """    
            else:
                if end not in dgv[chromosome][start]:
                    dgv[chromosome][start][end] = list()
                    dgv[chromosome][start][end].append(variantType)
                    dgv[chromosome][start][end].append(reference)
                else:
                    dgv[chromosome][start][end].append(variantType)
                    dgv[chromosome][start][end].append(reference)
            """
    #print("dgv = ", dgv)
    return dgv

def buildDGVHumanCNVDict(ff):
    fr = open(ff, "r")
    header = fr.readline()
    dgv = dict()
    for line in fr:
        tokens = line.strip().split("\t")
        chromosome = tokens[0][3:]
        print("in buildDGVHumanCNVDict(): chromosome = ", chromosome)
        start = int(tokens[1])
        end = int(tokens[2])
        variant = tokens[3]
        variantSubType = tokens[5] # Gain or Loss
        variantType = variant + " " + variantSubType
        reference = tokens[6]
        if chromosome not in dgv:
            dgv[chromosome] = dict()
            dgv[chromosome][start] = end
        else:
            if start not in dgv[chromosome]:
                dgv[chromosome][start] = end
    return dgv

def buildDGVGoldStandardDict(ff):
    fr = open(ff, "r")
    header = fr.readline()
    dgv = dict()
    for line in fr:
        tokens = line.strip().split("\t")
        chromosome = tokens[0].replace("chr", "")
        variant = tokens[1]
        start = int(tokens[3])
        end = int(tokens[4])
        variantSubType = tokens[8].split(";")[3].replace("variant_sub_type", "")
        variantType = variant + " " + variantSubType
    
        if chromosome not in dgv:
            dgv[chromosome] = dict()
            dgv[chromosome][start] = end
            #dgv[chromosome][start] = dict()
            #dgv[chromosome][start][end] = list()
            #dgv[chromosome][start][end].append(variantType)
            #dgv[chromosome][start][end].append(reference)
        else:
            if start not in dgv[chromosome]:
                dgv[chromosome][start] = end
                #dgv[chromosome][start][end] = list()
                #dgv[chromosome][start][end].append(variantType)
                #dgv[chromosome][start][end].append(reference)
            """    
            else:
                if end not in dgv[chromosome][start]:
                    dgv[chromosome][start][end] = list()
                    dgv[chromosome][start][end].append(variantType)
                    dgv[chromosome][start][end].append(reference)
                else:
                    dgv[chromosome][start][end].append(variantType)
                    dgv[chromosome][start][end].append(reference)
            """
    #print("dgv = ", dgv)
    return dgv

def buildExACDict(ff):
    fr = open(ff, "r")
    header = fr.readline()
    dgv = dict()
    for line in fr:
        tokens = line.strip().split()
        chromosome = tokens[1]
        start = int(tokens[2])
        end = int(tokens[3])
        gene_symbol = tokens[4]

        if chromosome not in dgv:
            dgv[chromosome] = dict()
            dgv[chromosome][start] = end
            #dgv[chromosome][start] = dict()
            #dgv[chromosome][start][end] = list()
            #dgv[chromosome][start][end].append(variantType)
            #dgv[chromosome][start][end].append(reference)
        else:
            if start not in dgv[chromosome]:
                dgv[chromosome][start] = end
                #dgv[chromosome][start][end] = list()
                #dgv[chromosome][start][end].append(variantType)
                #dgv[chromosome][start][end].append(reference)
            """
            else:
                if end not in dgv[chromosome][start]:
                    dgv[chromosome][start][end] = list()
                    dgv[chromosome][start][end].append(variantType)
                    dgv[chromosome][start][end].append(reference)
                else:
                    dgv[chromosome][start][end].append(variantType)
                    dgv[chromosome][start][end].append(reference)
            """
    #print("dgv = ", dgv)
    return dgv

def getTestedDBstart(dgv_sorted_starts, start):
    #print("dgv_sorted_starts = ", dgv_sorted_starts) 
    for i in range(len(dgv_sorted_starts)-1):
        if dgv_sorted_starts[i] == start:
            return dgv_sorted_starts[i],dgv_sorted_starts[i]
        elif dgv_sorted_starts[i+1] == start:
            return dgv_sorted_starts[i+1],dgv_sorted_starts[i+1]
        elif dgv_sorted_starts[i]< start and start <= dgv_sorted_starts[i+1]:
            #print("dgv_sorted_starts[i] = ", dgv_sorted_starts[i], "dgv_sorted_starts[i+1]", dgv_sorted_starts[i+1])
            return dgv_sorted_starts[i], dgv_sorted_starts[i+1]
    return -1,-1

def isOverlapped(chromosome, start, end, overlap_threshold, flag):
   
    if flag == "dgv" and chromosome not in dgv:
        return False, None
    elif flag == "dgv":
        sorted_starts = sorted(dgv[chromosome].keys())
    elif flag == "gencode":
        sorted_starts = sorted(gencodeStartEndMapping[chromosome].keys())
    
    tested_db_left_start, tested_db_right_start = getTestedDBstart(sorted_starts, start)

    #print("in isOverlapped(): chromosome = ", chromosome)
    if tested_db_left_start == -1:
        return False, None

    if flag == "dgv":
        tested_db_left_end = dgv[chromosome][tested_db_left_start]
        tested_db_right_end = dgv[chromosome][tested_db_right_start]
    elif flag == "gencode":
        tested_db_left_end = gencodeStartEndMapping[chromosome][tested_db_left_start]
        tested_db_right_end = gencodeStartEndMapping[chromosome][tested_db_right_start]

    if tested_db_left_start == tested_db_right_start:
        overlap, key, percent = isOverlappedHelper1(chromosome, start, end, overlap_threshold, tested_db_left_start, tested_db_left_end)
        return overlap, key, chromosome, tested_db_left_start, tested_db_left_end
    else:
        oL, kL, pL = isOverlappedHelper1(chromosome, start, end, overlap_threshold, tested_db_left_start, tested_db_left_end)
        oR, kR, pR = isOverlappedHelper1(chromosome, start, end, overlap_threshold, tested_db_right_start, tested_db_right_end)
        if (pL >= pR):
            return oL, kL, chromosome, tested_db_left_start, tested_db_left_end
        else:
            return oR, kR, chromosome, tested_db_right_start, tested_db_right_end
        
def isOverlappedHelper1(chromosome, start, end, overlap_threshold, tdb_start, tdb_end):
    overlap = False
    key = None
    overlap_size = 0
    percent_overlap = 0.0
    tdb_left_end = -1
    in_size = end - start

    if start >= tdb_start and end >= tdb_end: # (4) - (1) 
                                   #            (1)start                    (2)end 
                                   # (3)tdb__start              (4)tdb_end        
        overlap_size = tdb_end - start
        #print("hh1 overlap_size  = ", overlap_size )
    elif start <= tdb_start and end >= tdb_end: # (4) - (3) 
                                   # (1)start                                  (2)end 
                                   #        (3)tdb__start              (4)tdb_end        
        overlap_size = tdb_end - tdb_start
        # *** basically, overlap_size is bounded by tdb_end - tdb_start but, it I relax is and recalculate teh 
        # overlap_size to be end-start in the following, Apr 13, 2019
        overlap_size = end - start
        #print("hh2 overlap_size  = ", overlap_size )
    elif start >= tdb_start and end <= tdb_end: # (2) - (1)
                                   #        (1)start                    (2)end 
                                   #(3)tdb__start                           (4)tdb_end        
        overlap_size = end - start
        #print("hh3 overlap_size  = ", overlap_size )
    elif start <= tdb_start and end <= tdb_end: # (2) - (3) 
                                   #(1)start                    (2)end 
                                   #    (3)tdb__start                           (4)tdb_end        
        overlap_size = end -tdb_start
        #print("hh4 overlap_size  = ", overlap_size )
   
    percent_overlap = overlap_size*1.0/in_size

    if percent_overlap >= overlap_threshold:
        overlap = True
        key = chromosome + "_" + str(tdb_start)
    
    return overlap, key, percent_overlap

def updateCNVevents(chromosome, start, end, CNV_type):
    flatDict = dict()
    key = chromosome + "_" + str(start) + "_" + str(end)
    if CNV_type not in CNVevents:
        CNVevents[CNV_type] = dict()
        CNVevents[CNV_type][key] = 1
    else:
        if key not in CNVevents[CNV_type]:
            CNVevents[CNV_type][key] = 1
        else:
            CNVevents[CNV_type][key] += 1

def updateCNVsizes(chromosome, start, end, CNV_type):
    overlap_threshold = 0.5
    size = end - start
    dd = math.floor(size/1000)
    if dd < 1: dd = 1
    elif dd <5: dd = 5
    elif dd <10: dd = 10
    elif dd <20: dd = 20
    elif dd <30: dd = 30
    elif dd <40: dd = 40
    else:
        print(chromosome, start, end, CNV_type, "dd = ", dd)
        dd = 50

    if dd not in CNVsizes:
        CNVsizes[dd] = dict()
        CNVsizes[dd][CNV_type] = dict()
        CNVsizes[dd][CNV_type]['known'] = 0
        CNVsizes[dd][CNV_type]['novel'] = 0
        if isOverlapped(chromosome, start, end, overlap_threshold, "dgv")[0] == True:
            CNVsizes[dd][CNV_type]['known'] = 1
        else:
            CNVsizes[dd][CNV_type]['novel'] = 1
    else:
        if CNV_type not in CNVsizes[dd]:
            CNVsizes[dd][CNV_type] = dict()
            if isOverlapped(chromosome, start, end, overlap_threshold, "dgv")[0] == True:
                CNVsizes[dd][CNV_type]['known'] = 1
            else:
                CNVsizes[dd][CNV_type]['novel'] = 1
        else:
            if isOverlapped(chromosome, start, end, overlap_threshold, "dgv")[0] == True:
                if 'known' not in CNVsizes[dd][CNV_type]:
                    CNVsizes[dd][CNV_type]['known'] = 1
                else:
                    CNVsizes[dd][CNV_type]['known'] += 1
            else:
                if 'novel' not in CNVsizes[dd][CNV_type]:
                    CNVsizes[dd][CNV_type]['novel'] = 1
                else:
                    CNVsizes[dd][CNV_type]['novel'] +=1

def printCNVsizes():    
    fw = open("findCNVsizeDist_" + kit + ".csv","w")
    print("CNV size(kb)\tCNV_type\tKnown\tNovel")
    fw.write("CNV size(kb)\tCNV_type\tKnown\tNovel\n")    
    keys = CNVsizes.keys()

    x_axis = ["<1", "<1",  "<5", "<5", "<10", "<10", "<20", "<20","<30","<30","<40","<40",">40",">40"]
    dups_colors = ['F2D7D5','A93226','F2D7D5','A93226','F2D7D5','A93226','F2D7D5','A93226','F2D7D5','A93226','F2D7D5','A93226','F2D7D5','A93226']     
    dels_colors = ['D4E6F1','2471A3','D4E6F1','2471A3','D4E6F1','2471A3','D4E6F1','2471A3','D4E6F1','2471A3','D4E6F1','2471A3','D4E6F1','2471A3']
    dups_values = []
    dels_values = [] 
    for s in sorted(keys): # for different sizes
        print("s = ", s)
        #print("CNVsizes[s] = ", CNVsizes[s])
        types = CNVsizes[s].keys()

        for t in sorted(types): # for dels and dups
            print("t = ", t)
            print(str(s) + "\t" + t + "\t" + str(CNVsizes[s][t]['known']) +"\t" + str(CNVsizes[s][t]['novel']))
            if 'novel' in CNVsizes[s][t]:
                fw.write(str(s) + "\t" + t + "\t" + str(CNVsizes[s][t]['known']) +"\t" + str(CNVsizes[s][t]['novel']) + "\n")
            else:
                fw.write(str(s) + "\t" + t + "\t" + str(CNVsizes[s][t]['known']) +"\t" + "0" + "\n")

            # For Bokeh visualization, data preparation
            if t == 'dup':
                if t in CNVsizes[s]:
                    if 'novel' in CNVsizes[s]['dup']:
                        dups_values.append(CNVsizes[s]['dup']['novel'])
                    else:
                        dups_values.append(0)
                    if 'known' in CNVsizes[s]['dup']:
                        dups_values.append(CNVsizes[s]['dup']['known'])
                    else:
                        dups_values.append(0)
                else:
                    dups_values.append(0)
                    dups_values.append(0)
        
            elif t == 'del':
                if t in CNVsizes[s]:
                    if 'novel' in CNVsizes[s]['del']:
                        dels_values.append(CNVsizes[s]['del']['novel'])
                    else:
                        dels_values.append(0)
                    if 'known' in CNVsizes[s]['del']:
                        dels_values.append(CNVsizes[s]['del']['known'])
                    else:
                        dels_values.append(0)
                else:
                    dels_values.append(0)
                    dels_values.append(0)
                
            # Fnd for Bokeh visualization, data preparation     
    fw.close()

    print("dups_values = ", dups_values)
    print("dels_values = ", dels_values)

    # For Bokeh visualization
    df1 = pd.DataFrame([x_axis,dups_values,dups_colors]).T
    print("before pivot df1 = ", df1)
    df1.columns = ['QID','score', 'DocID']
    df1 = df1.pivot(index='QID', columns='DocID', values='score').fillna(0)
    print("after pivot df1 = ", df1)
    df1.index = [(x, 'amp') for x in df1.index]
    print("df1.index = ", df1.index)

    df2 = pd.DataFrame([x_axis,dels_values,dels_colors]).T
    df2.columns = ['QID','score', 'DocID']
    df2 = df2.pivot(index='QID', columns='DocID', values='score').fillna(0)
    print("df2 = ", df2)
    df2.index = [(x, 'del') for x in df2.index]
    print("df2.index = ", df2.index)

    df = pd.concat([df1,df2])
    df = df.fillna(0)
    print("df = ", df)
    print("df.columns = ", df.columns)
    print("df.index = ", df.index)

    column_order = [('<1', 'amp'), ('<5', 'amp'), ('<10', 'amp'), ('<20', 'amp'), ('<30', 'amp'), ('<40', 'amp'), ('>40', 'amp'),  ('<1', 'del'),('<5', 'del'), ('<10', 'del'), ('<20', 'del'), ('<30', 'del'), ('<40', 'del'),('>40', 'del')]

    print("[value(x) for x in df.columns] = ", [value(x) for x in df.columns])
    colors_set = ["#D4E6F1",  "#F2D7D5", "#2471A3", "#A93226"]
    #p = figure(plot_width=1200, x_range=FactorRange(*df.index))

    p = figure(plot_height=350, plot_width=550, x_range=FactorRange(factors=column_order), toolbar_location=None, title="CNV size Distribution")
    print('p = ', p)
    r = p.vbar_stack(df.columns, x='index', width=1.0, fill_color=colors_set,
            line_color=None, source=df)

    """
    li1 = LegendItem(label='Amplification Known', renderers=[r[1]])  
    li2 = LegendItem(label='Amplification Novel', renderers=[r[3]])
    li3 = LegendItem(label='Deletion Known', renderers=[r[0]])
    li4 = LegendItem(label='Deletion Novel', renderers=[r[2]])
    #legend1 = Legend(items=[li1, li2, li3, li4], location='top_left', glyph_height=10, glyph_width=10)
    legend1 = Legend(items=[li1, li2, li3, li4], location='center')
    p.add_layout(legend1) 
    """

    #p.xgrid.grid_line_color = None
    #p.ygrid.grid_line_color = None
    #p.xaxis.subgroup_label = None
    #p.subgroup_text_font_size = 1
    p.xaxis.axis_label_text_font_size = "15pt"
    p.yaxis.axis_label_text_font_size = "15pt"
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"
    #p.legend.label_text_font_size = "10px"
    p.xaxis.group_text_font_size = "12pt"
    p.xaxis.major_tick_line_color = None
    p.xaxis.major_label_text_color = None
    p.xaxis.axis_label = "CNV size (kb)"
    p.yaxis.axis_label = "Number of events"    
    output_file("CNVsizesDistribution_" + kit + ".html")
    export_png(p, filename="CNVsizesDistribution_" + kit + ".png")
    show(p)

def printCNVevents(CNVeventsByFreq):
    fw = open("findCNVeventDist_" + kit + ".csv","w")
    print("CNV type\tCNV event num bin\tKnown\tNovel")
    fw.write("CNV type\tCNV event num bin\tKnown\tNovel\n")    
    keys = CNVeventsByFreq.keys() # dup and del only
    
    num_events = ["Singleton", "2-5", "6-25", "26-50", "51-100", ">100"]
    KKNN = ["Amplification known", "Amplification novel", "Deletion known", "Deletion novel"]
    KN = ["Known", "Novel"]
    Dups = {"num_events": num_events, 'Known':[], 'Novel': []}
    Dels = {"num_events": num_events, 'Known':[], 'Novel': []}

    #for s in sorted(keys):
    bins = CNVeventsByFreq['dup'].keys()
    #print('bins =', bins)
    for t in sorted(bins): #1, 2, 6, 26, 51, 100
        print('dup' + "\t" + str(t) + "\t" + str(CNVeventsByFreq['dup'][t]['known']) +"\t" + str(CNVeventsByFreq['dup'][t]['novel']))
        
        fw.write(str(t) + "\t" + 'dup' + "\t" + str(CNVeventsByFreq['dup'][t]['known']) +"\t" + str(CNVeventsByFreq['dup'][t]['novel']) + "\n")
        #print('del' + "\t" + str(t))
        #print(str(CNVeventsByFreq['del'][t]['known']))
        #print(str(CNVeventsByFreq['del'][t]['novel']))
        if t in CNVeventsByFreq['del']:
            fw.write(str(t) + "\t" + 'del' + "\t" + str(CNVeventsByFreq['del'][t]['known']) +"\t" + str(CNVeventsByFreq['del'][t]['novel']) + "\n")
        else:
            fw.write(str(t) + "\t" + 'del' + "\t" + "0" +"\t" + "0" + "\n")
            
        if t in CNVeventsByFreq['dup']:
            if 'known' in CNVeventsByFreq['dup'][t]:
                Dups['Known'].append(CNVeventsByFreq['dup'][t]['known'])
            else:
                Dups['Known'].append(0)
            if 'novel' in CNVeventsByFreq['dup'][t]:
                Dups['Novel'].append(CNVeventsByFreq['dup'][t]['novel'])
            else:
                Dups['Novel'].append(0)
        else:
            Dups['Known'].append(0)
            Dups['Novel'].append(0)
        if t in CNVeventsByFreq['del']:
            if 'known' in CNVeventsByFreq['del'][t]:
                Dels['Known'].append(-CNVeventsByFreq['del'][t]['known'])
            else:
                Dels['Known'].append(0)
            if 'novel' in CNVeventsByFreq['del'][t]:
                Dels['Novel'].append(-CNVeventsByFreq['del'][t]['novel'])
            else:
                Dels['Novel'].append(0)
        else:
            Dels['Known'].append(0)
            Dels['Novel'].append(0)
        
    fw.close()
    #data = {'num_events': num_events,
    #    'Duplication known' : DupKnowns,
    #    'Duplication novel' : DupNovels,
    #    'Deletion known' : DelKnowns,
    #    'Deletion novel' : DelNovels}
   
    print("Dups = ", Dups)
    print("Dels = ", Dels)
    while len(Dups['Known']) < len(num_events):
        Dups['Known'].append(0)
    while len(Dups['Novel']) < len(num_events):
        Dups['Novel'].append(0)
    while len(Dels['Known']) < len(num_events):
        Dels['Known'].append(0)
    while len(Dels['Novel']) < len(num_events):
        Dels['Novel'].append(0)
    print("Dups = ", Dups)
    print("Dels = ", Dels)
    
    output_file("CNVeventsDistribution_" + kit + ".html")
    
    colors1 = ["#F2D7D5", "#A93226"]
    colors2 = ["#D4E6F1", "#2471A3"]
    #source = ColumnDataSource(data=data)
    #p = figure(x_range=num_events, plot_height=350, title="Number of CNV  Distribution",toolbar_location=None, tools="")
    #p.vbar_stack(KKNN, x='num_events' , width=0.9, color=colors, source=source,
                         #legend=[value(x) for x in KKNN], name=KKNN)
    #p = figure(y_range=num_events, plot_height=350, x_range=(-max(Dups['Known']))-100, max(Dups['Known'])+100), title="Number of CNV Distribution", toolbar_location=None)
    print("max(Dels['Known'])+max(Dels['Novel'])",-(max(Dels['Known'])+max(Dels['Novel']))-100)
    print("Dels['Known'] = ", Dels['Known'])
    print("Dels['Novel'] = ", Dels['Novel'])
    p = figure(y_range=num_events, plot_height=350, x_range=((min(Dels['Known'])+min(Dels['Novel']))-100, max(Dups['Known'])+max(Dups['Novel'])+100), title="Number of CNV Distribution", toolbar_location=None)
    #p = figure(y_range=num_events, plot_height=350, x_range=(-10000, + 10000), title="Number of CNV Distribution", toolbar_location=None)
    p.hbar_stack(KN, y='num_events', height=0.9, color=colors1, source=ColumnDataSource(Dups),legend=["Amplification %s" % x for x in KN])
    p.hbar_stack(KN, y='num_events', height=0.9, color=colors2, source=ColumnDataSource(Dels),legend=["Deletion %s" % x for x in KN])

    #ymax = min(Dups['Known'])
    #ymin = max(Dups['Known'])
    #p.y_range = Range1d(ymin, ymax)
    p.xaxis.major_label_overrides = {-4000: '4000', -2000: '2000'}
    p.xaxis.axis_label_text_font_size = "15pt"
    p.yaxis.axis_label_text_font_size = "15pt"
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"
    #p.legend.label_text_font_size = "8px"
    p.y_range.range_padding = 0.1
    p.ygrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    p.xaxis.axis_label = "Number of CNVs"
    p.yaxis.axis_label = "Number of events"
    p.outline_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "vertical"
    export_png(p, filename="CNVeventsDistribution_" + kit + ".png")
    #p.output_backend = "svg"
    #export_svgs(p, filename="CNVeventsDistribution_" + kit + ".svg")
    show(p)

def sortCNVeventsByFreq(sortedCNVeventsByFreqDUP, fw, genename_mapping, genename_tracking, cnvtype):
    known_genes = dict()
    novel_genes = dict()
    known_genes_wo_event_in_key = dict()
    novel_genes_wo_event_in_key = dict()
    for event, freq in sortedCNVeventsByFreqDUP:
        tokens = event.split("_")
        chromosome = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        genename = "***" + chromosome + "_"+ str(start) + "_" + str(end)
        fgenename = "***" + chromosome + "_"+ str(start) + "_" + str(end)
        overlap, mapped_key, chromosome, gencode_start, gencode_end =  isOverlapped(chromosome, start, end, GENCODE_overlap_threshold, "gencode")
        if mapped_key is not None:
            genename = gencodeStartAnnotationMapping[mapped_key]
            #genename += "GENCODE"
            fgenename = gencodeStartAnnotationMapping[mapped_key]
            fgenename += "[GENCODE:" + chromosome + "_" + str(gencode_start) + "_" + str(gencode_end) + "]"
        overlap, mapped_key, chromsome, dgv_start, dgv_end = isOverlapped(chromosome, start, end, 0.5, "dgv")
        if mapped_key is not None:
            genename += "[known]"
            fgenename += "[known]"
            fgenename += "[DGV:" + chromsome + "_"+ str(dgv_start) + "_" + str(dgv_end) + "]"
        else:
            genename += "[novel]"
            fgenename += "[novel]"
        fw.write(event + "\t" + cnvtype + "\t" + fgenename + "\t" + str(freq) + "\n")
        #fw.write("dup\t" + genename + "\t" + str(freq) + "\n")
        if genename.find("known")>=0:
            genename1 = genename.split("|")[0]
            genename1 = genename1.replace("[known]", "")
            genename1 = genename1.replace("***", "")
            genename_only = genename1
            genename1 = genename1 + "_" + event 
            if genename1 in known_genes:
                known_genes[genename1] += freq
            else:
                known_genes[genename1] = freq

            if genename1 in genename_mapping:
                genename_only = genename_mapping[genename1]
                known_genes_wo_event_in_key[genename_only] = freq
            else:
                # Note that the tracking here represents the CNV with chr, start, end that falls within the same 
                # gene but we treated and count them independently....  as genename, genename.1, genename.2 and so on...
                if genename_only not in genename_tracking:
                    known_genes_wo_event_in_key[genename_only] = freq
                    genename_tracking[genename_only] =  1
                    genename_mapping[genename1] = genename_only
                else:
                    known_genes_wo_event_in_key[genename_only + "." + str(genename_tracking[genename_only])] = freq
                    genename_mapping[genename1] = genename_only + "." + str(genename_tracking[genename_only])
                    genename_tracking[genename_only] += 1
            """
            if genename1 in known_genes:
                if genename_only not in genename_tracking:
                    known_genes_wo_event_in_key[genename_only] = freq
                    genename_tracking[genename_only] =  1
                    genename_mapping[genename1] = genename_only
                else:
                    known_genes_wo_event_in_key[genename_only + "." + str(genename_tracking[genename_only])] = freq
                    genename_mapping[genename1] = genename_only + "." + str(genename_tracking[genename_only]) 
                    genename_tracking[genename_only] += 1
            else:
                known_genes_wo_event_in_key[genename_only] = freq
                genename_tracking[genename_only] =  1
                genename_mapping[genename1] = genename_wo_event
            """
        if genename.find("novel")>=0:
            genename1 = genename.split("|")[0]
            genename1 = genename1.replace("[novel]", "")
            genename1 = genename1.replace("***", "")
            genename_only = genename1
            genename1 = genename1 + "_" + event

            if genename1 in novel_genes:
                novel_genes[genename1] += freq
            else:
                novel_genes[genename1] = freq
            if genename1 in genename_mapping:
                genename_only = genename_mapping[genename1]
                novel_genes_wo_event_in_key[genename_only] = freq
            else:
                if genename_only not in genename_tracking:
                    novel_genes_wo_event_in_key[genename_only] = freq
                    genename_tracking[genename_only] =  1
                    genename_mapping[genename1] = genename_only
                else:
                    novel_genes_wo_event_in_key[genename_only + "." + str(genename_tracking[genename_only])] = freq
                    genename_mapping[genename1] = genename_only + "." + str(genename_tracking[genename_only])
                    genename_tracking[genename_only] += 1
            """
            if genename1 in novel_genes:
                if genename_only not in genename_tracking:
                    novel_genes_wo_event_in_key[genename_only] = freq
                    genename_tracking[genename_only] =  1
                    genename_mapping[genename1] = genename_only
                else:
                    novel_genes_wo_event_in_key[genename_only + "." + str(genename_tracking[genename_only])] = freq
                    genename_mapping[genename1] = genename_only + "." + str(genename_tracking[genename_only])
                    genename_tracking[genename_only] += 1
                    
            else:
                novel_genes_wo_event_in_key[genename_only] = freq
                genename_tracking[genename_only] =  1 
                genename_mapping[genename1] = genename_only
            """
    return known_genes, novel_genes, known_genes_wo_event_in_key, novel_genes_wo_event_in_key
    #return known_genes, novel_genes

def drawFigure(color_1, color_2, genelist_text, genelist, dups, dels,KnownOrNovel,title_1, outfile_name):

    #=====================================================================
    ## Bokeh visualization for Top known genes with both Dups and Dels 
    #=====================================================================
    color1 = color_1
    color2 = color_2
    KN = [KnownOrNovel]
    Dups = {genelist_text: genelist,
            KnownOrNovel: dups}
    Dels = {genelist_text: genelist,
            KnownOrNovel: dels}
    #KK = [KnownOrNovel]
    p = figure(x_range=genelist, plot_height=400, plot_width=550, toolbar_location=None, title=title_1,tools="")
    p.vbar_stack(KN, x=genelist_text, width=0.9, color=color1, source=ColumnDataSource(Dups))
    p.vbar_stack(KN, x=genelist_text, width=0.9, color=color2, source=ColumnDataSource(Dels))

    # START comment by Dec 10, 2019
    #print("genelist_text => == ", genelist_text)
    if genelist_text == 'sorted_known_genes_with_both_dup_and_del':
        p.yaxis.major_label_overrides = {-300: '300', -200: '200', -100: '100'}
    elif genelist_text == 'sorted_novel_genes_with_both_dup_and_del':
        p.yaxis.major_label_overrides = {-100: '100', -50: '50'}
    # END comment by Dec 10, 2019
    #p.y_range.start = 0
    #print("p.yaxis.major_label_overrides = ", p.yaxis.major_label_overrides)
    #p.yaxis.major_label_overrides = {-400: '400', -200: '200'}
    #p.yaxis.major_label_overrides = {-100: '100', -50: '50'}
    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    p.xaxis.axis_label_text_font_size = "15pt"
    p.yaxis.axis_label_text_font_size = "15pt"
    p.xaxis.major_label_text_font_size = "10pt"
    p.yaxis.major_label_text_font_size = "10pt"
    #p.xaxis.axis_label = ""
    p.yaxis.axis_label = "Number of events"
    p.xaxis.major_label_orientation = math.pi/3
    p.outline_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "vertical"

    output_file(outfile_name + kit + ".html")
    export_png(p, filename=outfile_name + kit + ".png")
    show(p)

def mapCNVeventsByGene():
    genename_mapping = dict()
    genename_tracking = dict()
    toporders = 25
    fw = open("mapCNVeventsByGene_" + kit + ".csv", "w")
    
    sortedCNVeventsByFreqDUP = sorted(CNVevents['dup'].items(), key=lambda x:x[1], reverse=True)
    known_genes_with_dup, novel_genes_with_dup, known_genes_with_dup_wo_event, novel_genes_with_dup_wo_event  = sortCNVeventsByFreq(sortedCNVeventsByFreqDUP, fw, genename_mapping, genename_tracking, "dup")
    
    sortedCNVeventsByFreqDEL = sorted(CNVevents['del'].items(), key=lambda x:x[1], reverse=True)
    known_genes_with_del, novel_genes_with_del, known_genes_with_del_wo_event, novel_genes_with_del_wo_event  = sortCNVeventsByFreq(sortedCNVeventsByFreqDEL, fw, genename_mapping, genename_tracking, "del")
    fw.close()

    #===================================================
    # Build file for Top known genes with Dups and Dels
    #===================================================
    """
    known_genes_with_both_dup_and_del = list()
    KnownDupDel = dict() # add June 24, 2018
    for a in known_genes_with_dup.keys():
        if a in known_genes_with_del:
            known_genes_with_both_dup_and_del.append(a)
            KnownDupDel[a] = known_genes_with_dup[a] + known_genes_with_del[a] # add June 24, 2018
    DupKnowns = []
    DelKnowns = []
    #sorted_known_genes_with_both_dup_and_del = []

    sortedKnownDupDel = sorted(KnownDupDel.items(), key=lambda x:x[1], reverse=True) # add June 24, 2018
    # The sorted_known_genes_with_both_dup_and_del here will be used in the Visualization later...
    # However, the gene name comes with the event (chromosome_start_end), too long,
    # So, I build another list which mimick the sorted_known_genes_with_both_dup_and_del but with the adjusted
    # gene names
    sorted_known_genes_with_both_dup_and_del = [x for x,y in sortedKnownDupDel]
    sorted_known_genes_with_both_dup_and_del = sorted_known_genes_with_both_dup_and_del[:toporders]
    #print("sorted_known_genes_with_both_dup_and_del = ", sorted_known_genes_with_both_dup_and_del)

    fw30 = open("TopKnownGenesWithBothDupsAndDels.txt", "w")
    fw30.write("Top known genes w/ CNVs\tCounts\n")
    count = 1
    for g, c in sortedKnownDupDel:
        fw30.write(genename_mapping[g] + "\t" + str(known_genes_with_dup[g]) + "\t" + genename_mapping[g] + "\t" + str(known_genes_with_del[g]) + "\t" + g + "\t"+ str(c) + "\n")
        count += 1
        if count> toporders:
            break
    fw30.close()
     
    mm = 1
    for a, freq in sortedKnownDupDel:
        if mm > toporders: break
        DupKnowns.append(known_genes_with_dup[a])
        DelKnowns.append(-known_genes_with_del[a])
        mm += 1
   
    drawFigure("#F2D7D5","#D4E6F1","sorted_known_genes_with_both_dup_and_del", sorted_known_genes_with_both_dup_and_del, DupKnowns, DelKnowns, "Known", "Top Known Genes with Both Deletions and Amplifications", "CNVKnownGenesWithBothDupsandDels_")
    """
    #================================================
    # FOR SHORTER GENE NAMES 
    #================================================
    known_genes_with_both_dup_and_del = list()
    KnownDupDel = dict() # add June 24, 2018
    #print("HELLO known_genes_with_dup_wo_event = ",known_genes_with_dup_wo_event)
    #print("HELLO known_genes_with_del_wo_event = ", known_genes_with_del_wo_event) 
    for a in known_genes_with_dup_wo_event.keys(): #####
        if a in known_genes_with_del_wo_event: ####
            known_genes_with_both_dup_and_del.append(a)
            KnownDupDel[a] = known_genes_with_dup_wo_event[a] + known_genes_with_del_wo_event[a] #####
    DupKnowns = []
    DelKnowns = []
    #sorted_known_genes_with_both_dup_and_del = []

    sortedKnownDupDel = sorted(KnownDupDel.items(), key=lambda x:x[1], reverse=True) # add June 24, 2018
    # The sorted_known_genes_with_both_dup_and_del here will be used in the Visualization later...
    # However, the gene name comes with the event (chromosome_start_end), too long,
    # So, I build another list which mimick the sorted_known_genes_with_both_dup_and_del but with the adjusted
    # gene names
    sorted_known_genes_with_both_dup_and_del = [x for x,y in sortedKnownDupDel]
    sorted_known_genes_with_both_dup_and_del = sorted_known_genes_with_both_dup_and_del[:toporders]
    #print("sorted_known_genes_with_both_dup_and_del = ", sorted_known_genes_with_both_dup_and_del)
   
    mm = 1
    for a, freq in sortedKnownDupDel:
        if mm > toporders: break
        DupKnowns.append(known_genes_with_dup_wo_event[a]) #####
        DelKnowns.append(-known_genes_with_del_wo_event[a]) #####
        mm += 1

    drawFigure("#F2D7D5","#D4E6F1","sorted_known_genes_with_both_dup_and_del", sorted_known_genes_with_both_dup_and_del, DupKnowns, DelKnowns, "Known", "Top Known Genes with Both Deletions and Amplifications", "CNVKnownGenesWithBothDupsandDels_")

    #================================================
    #### for novel genes with Dups and Dels
    #================================================
    """
    novel_genes_with_both_dup_and_del = list()
    NovelDupDel = dict()
    for a in novel_genes_with_dup.keys():
        if a in novel_genes_with_del:
            novel_genes_with_both_dup_and_del.append(a)
            NovelDupDel[a] = novel_genes_with_dup[a] + novel_genes_with_del[a] # add June 24, 2018
    DupNovels = []
    DelNovels = []
    #sorted_novel_genes_with_both_dup_and_del = []

    #print("novel_genes_with_both_dup_and_del = ", novel_genes_with_both_dup_and_del)
    sortedNovelDupDel = sorted(NovelDupDel.items(), key=lambda x:x[1], reverse=True) # add June 24, 2018
    sorted_novel_genes_with_both_dup_and_del = [x for x,y in sortedNovelDupDel]
    print("long genename : sortedNovelDupDel = ", sortedNovelDupDel)
    fw30 = open("TopNovelGenesWithBothDupsAndDels.txt", "w")
    fw30.write("Top novel genes w/ CNVs\tCounts\n")
    count = 1
    for g, c in sortedNovelDupDel:
        fw30.write(genename_mapping[g] + "\t" + str(novel_genes_with_dup[g]) + "\t" + genename_mapping[g] + "\t" + str(novel_genes_with_del[g]) + "\t" + g + "\t"+ str(c) + "\n")
        count += 1
        if count> toporders:
            break
    fw30.close()

    mm = 1 
    for a, freq in sortedNovelDupDel:
        if mm > toporders: break
        DupNovels.append(novel_genes_with_dup[a])
        DelNovels.append(-novel_genes_with_del[a])
        print(mm, a, freq)
        mm += 1

    drawFigure("#A93226","#2471A3","sorted_novel_genes_with_both_dup_and_del", sorted_novel_genes_with_both_dup_and_del, DupNovels, DelNovels, "Novel", "Top Novel Genes with Both Deletions and Amplifications", "CNVNovelGenesWithBothDupsandDels_")
    """
    #======================================
    # FOR SHORTER GENE NAMES
    #======================================
    novel_genes_with_both_dup_and_del = list()
    NovelDupDel = dict() # add June 24, 2018
    for a in novel_genes_with_dup_wo_event.keys(): #####
        if a in novel_genes_with_del_wo_event: ####
            novel_genes_with_both_dup_and_del.append(a)
            NovelDupDel[a] = novel_genes_with_dup_wo_event[a] + novel_genes_with_del_wo_event[a] #####
    DupNovels = []
    DelNovels = []
    #sorted_known_genes_with_both_dup_and_del = []

    sortedNovelDupDel = sorted(NovelDupDel.items(), key=lambda x:x[1], reverse=True) # add June 24, 2018
    print("short gene name : sortedNovelDupDel = ", sortedNovelDupDel)
    # The sorted_known_genes_with_both_dup_and_del here will be used in the Visualization later...
    # However, the gene name comes with the event (chromosome_start_end), too long,
    # So, I build another list which mimick the sorted_known_genes_with_both_dup_and_del but with the adjusted
    # gene names
    sorted_novel_genes_with_both_dup_and_del = [x for x,y in sortedNovelDupDel]
    sorted_novel_genes_with_both_dup_and_del = sorted_novel_genes_with_both_dup_and_del[:toporders]
    print("sorted_novel_genes_with_both_dup_and_del = ", sorted_novel_genes_with_both_dup_and_del)

    mm = 1
    for a, freq in sortedNovelDupDel:
        if mm > toporders: break
        DupNovels.append(novel_genes_with_dup_wo_event[a]) #####
        DelNovels.append(-novel_genes_with_del_wo_event[a]) #####
        mm += 1

    drawFigure("#A93226","#2471A3","sorted_novel_genes_with_both_dup_and_del", sorted_novel_genes_with_both_dup_and_del, DupNovels, DelNovels, "Novel", "Top Novel Genes with Both Deletions and Amplifications", "CNVNovelGenesWithBothDupsandDels_")

    #### Build graph novel genes with Dups (Graph3)
    DupNovels = []
    genes = []
    sortedNovelDup = sorted(novel_genes_with_dup.items(), key=lambda x:x[1], reverse=True) # add June 24, 2018

    #for a in novel_genes_with_dup:
    mm = 1
    for a,freq in sortedNovelDup:
        if mm > toporders:break
        genes.append(a)
        DupNovels.append(novel_genes_with_dup[a])
        mm += 1
    #print("genes = ", genes)
    #print("DupNovels = ", DupNovels)
    data = {'genes': genes,
            'DupNovels': DupNovels}

    source = ColumnDataSource(data=data)

    #p = figure(x_range=genes, plot_height=350, title="Top novel genes with duplications")
    p = figure(x_range=genes, plot_height=600, plot_width=550, toolbar_location=None, title="Top Novel Genes with Duplications")
    #p.vbar(x='genes', top='DupNovels', width=0.9, color=color1, source=source)
    p.vbar(x=genes, top=DupNovels, width=0.9, color="#A93226")

    p.x_range.range_padding = 0.1
    p.xaxis.axis_label_text_font_size = "15pt"
    p.yaxis.axis_label_text_font_size = "15pt"
    p.xaxis.major_label_text_font_size = "10pt"
    p.yaxis.major_label_text_font_size = "10pt"
    p.xgrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    #p.xaxis.axis_label = ""
    p.yaxis.axis_label = "Number of events"
    p.xaxis.major_label_orientation = math.pi/3
    p.outline_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "vertical"

    output_file("CNVNovelGenesWithDups_" + kit + ".html")
    export_png(p, filename="CNVNovelGenesWithDups_" + kit + ".png")
    #p.output_backend = "svg"
    #export_svgs(p, filename="CNVNovelGenesWithDups_" + kit + ".svg")
    show(p)

    #### Build graph novel genes with Dups (Graph4)
    DelNovels = []
    genes = []
    sortedNovelDel = sorted(novel_genes_with_del.items(), key=lambda x:x[1], reverse=True) # add June 24, 2018
    mm = 1
    for a,freq in sortedNovelDel:
        if mm > toporders:break
        genes.append(a)
        DelNovels.append(novel_genes_with_del[a])
        mm +=1
    #print("genes = ", genes)
    #print("DelNovels = ", DelNovels)
    data = {'genes': genes,
            'DelNovels': DelNovels}

    source = ColumnDataSource(data=data)

    #p = figure(x_range=genes, plot_height=350, title="Top novel genes with deletions")
    p = figure(x_range=genes, plot_height=600, plot_width=550, toolbar_location=None, title="Top Novel Genes with Deletions")
    #p.vbar(x='genes', top='DupNovels', width=0.9, color=color1, source=source)
    p.vbar(x=genes, top=DelNovels, width=0.9, color="#2471A3")

    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    p.xaxis.axis_label_text_font_size = "15pt"
    p.yaxis.axis_label_text_font_size = "15pt"
    p.xaxis.major_label_text_font_size = "10pt"
    p.yaxis.major_label_text_font_size = "10pt"
    #p.xaxis.axis_label = ""
    p.yaxis.axis_label = "Number of events"
    p.xaxis.major_label_orientation = math.pi/3
    p.outline_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "vertical"

    output_file("CNVNovelGenesWithDels_" + kit + ".html")
    export_png(p, filename="CNVNovelGenesWithDels_" + kit + ".png")
    #p.output_backend = "svg"
    #export_svgs(p, filename="CNVNovelGenesWithDels_" + kit + ".svg")
    show(p)

def groupCNVeventsByFrequency():
    overlap_threshold = 0.5
    keys = CNVevents.keys()
    CNVeventsByFreq = dict()
    for k in keys: # here, there are only two keys: 'del' or 'dup'
        print("k = ", k)
        print("CNVevents[" + k +  "].keys() = ", len(CNVevents[k].keys()))
        CNVeventsByFreq[k] = dict()
        for event, freq in CNVevents[k].items():
            print("event = ", event)
            print("freq = ", freq)
            chromosome, start, end = event.split("_")
            start = int(start)
            end = int(end)
            if freq == 1:
                ff = 1
            elif 2<=freq<=5:
                ff = 2
            elif 6 <= freq <= 25:
                ff = 6
            elif 26 <= freq <=50:
                ff = 26
            elif 51 <= freq <= 100:
                ff = 51
            elif freq > 100:
                ff = 100

            if ff not in CNVeventsByFreq[k]:
                CNVeventsByFreq[k][ff] = dict()
                CNVeventsByFreq[k][ff]['known'] = 0
                CNVeventsByFreq[k][ff]['novel'] = 0
                if isOverlapped(chromosome, start, end, overlap_threshold, "dgv")[0] == True:
                    CNVeventsByFreq[k][ff]['known'] = 1
                else:
                    CNVeventsByFreq[k][ff]['novel'] = 1
            else:
                if isOverlapped(chromosome, start, end, overlap_threshold, "dgv")[0] == True:
                    CNVeventsByFreq[k][ff]['known'] += 1
                else:
                    CNVeventsByFreq[k][ff]['novel'] += 1
    print(CNVeventsByFreq)
    return CNVeventsByFreq

def buildSampleToSource(ff):
    sample2source = dict()
    fr = open(ff, 'r')
    for filename in fr:
        source = filename.strip().replace("_sampname", "")
        f1 = open(filename.strip(), "r")
        for line in f1:
            sample = line.strip()
            sample2source[sample] = source
        f1.close()
    fr.close()
    return sample2source
  
def buildSampId2Name(ff):
    sampId2Name = dict()
    fr = open(ff, "r")
    count = 1
    for line in fr:
        sampName = line.strip()
        ID = 'sample' + str(count)
        sampId2Name[ID] = sampName
        print("ID = ", ID, "sampName = ", sampName)
        count += 1
    fr.close()
    return sampId2Name

#=========
# main
#=========
gencodeStartEndMapping, gencodeStartAnnotationMapping = buildGENCODEDict(gencodefile)
#ID2Name = buildSampId2Name(sampNameFile)

if "GS" in dgvfile:
    dgv = buildDGVGoldStandardDict(dgvfile)
elif "exac" in dgvfile:
    dgv = buildExACDict(dgvfile)
elif "Stringent.Gain+Loss" in dgvfile:
    dgv = buildDGVHumanCNVDict(dgvfile)
else:
    dgv = buildDGVDict(dgvfile)

sample2source = buildSampleToSource(sampleSourceFile)

CNVsizes = dict()
CNVevents = dict()
CNVeventsReduced = dict()

samples = []
cnvtypes = dict()
known_count = 0
novel_count = 0
CNV2chr = dict()
sample2CNVs = dict()

chrlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14','15','16','17','18','19','20','21','22', 'Y', 'X']

FULL = True

fr = open(infile, "r")
for line in fr:
    filename = workingDir + line.strip()
    if DEBUG:
        print("filename = ", filename)
    fr1 = open(filename, "r")
    # Skip header
    fr1.readline()
    for fr1_line in fr1:
        content_tokens = fr1_line.strip().split("\t")
        #sample_running = content_tokens[0]
        #sample = ID2Name[sample_running]
        sample = content_tokens[0]
        #sample = sample.split(".")[0] # Just added Jun 26, 2019
        if not sample in samples:
            samples.append(sample)
        chromosome = content_tokens[1]
        chromosome = chromosome.replace("chr", "")
        CNV_type = content_tokens[2]
        if CNV_type not in cnvtypes:
            cnvtypes[CNV_type] = 1
        else:
            cnvtypes[CNV_type] += 1
        #print("CNV_type = ", CNV_type)
        start = int(content_tokens[3])
        end = int(content_tokens[4])
        size = end - start 
        #print("sample, CNV_type, chr, start, end, size = ", sample, CNV_type, chromosome, start, end, size)
        if FULL == True:
            updateCNVsizes(chromosome, start, end, CNV_type)
        if sample not in sample2CNVs:
            sample2CNVs[sample] = dict()
            sample2CNVs[sample][chromosome] = 1
        else:
            if chromosome not in sample2CNVs[sample]:
                sample2CNVs[sample][chromosome] = 1
            else:
                sample2CNVs[sample][chromosome] += 1
        # For counting the overall novels
        """
        overlap, key = isOverlapped(chromosome, start, end, 0.5, "dgv")
        if overlap == True:
            known_count += 1
        else:
            novel_count += 1
        """
        # End For counting the overall novels

        if chromosome not in CNV2chr:
            CNV2chr[chromosome] = 1
        else:
            CNV2chr[chromosome] += 1
        updateCNVevents(chromosome, start, end, CNV_type)

fr.close()

# Counting all CNVevents
fw20 = open("CNVevents_details_" + kit + ".csv", "w")
fw20.write("CNV type\tKey\tChromosome\tStart\tEnd\tKnown or Novel\tDGV\t\tDiff size\tCounts\tGene\GENCODE\n")
total_events = 0
total_events_known = 0
total_events_novel = 0

types = CNVevents.keys()
for t in types:
    keys = CNVevents[t].keys()
    for k in keys:
        string_out = ""
        tokens = k.split("_")
        chromosome = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        returns = isOverlapped(chromosome, start, end, 0.5, "dgv")
        if len(returns) == 5:
            overlap, dgv_mapped_key, chromosome, dgv_start, dgv_end =  isOverlapped(chromosome, start, end, 0.5, "dgv")
            if (overlap == True):
                string_out = t + "\t" + k + "\t" + chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + "Known" + "\t" + str(end-start) + "\t" + str(CNVevents[t][k]) + "\t" + "[DGV:" + chromosome + "_" + str(dgv_start) + "_" + str(dgv_end) + "]\t"
                total_events_known += CNVevents[t][k]
            else:    
                string_out = t + "\t" + k + "\t" + chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + "Novel"  + "\t" + str(end-start) + "\t" + str(CNVevents[t][k]) + "\t" + "None" + "\t"
                total_events_novel += CNVevents[t][k]
            overlap, gencode_mapped_key, chromosome, gencode_start, gencode_end =  isOverlapped(chromosome, start, end, GENCODE_overlap_threshold, "gencode")
            if (overlap == True):
                string_out += gencodeStartAnnotationMapping[gencode_mapped_key] + "\t" + "[DGV:" + chromosome + "_" + str(dgv_start) + "_" + str(dgv_end) + "]\t"
            else:
                string_out += "None\tNone"
            string_out += "\n"
            fw20.write(string_out)
            total_events += CNVevents[t][k]
    
fw20.write("Total events known" + "\t" + str(total_events_known) + "\n")
fw20.write("Total events novel" + "\t" +str(total_events_novel) + "\n")
fw20.write("Total events" + "\t" +str(total_events) + "\n")
fw20.close()
#sys.exit()

if FULL == True:
    printCNVsizes()
    
CNVeventByFreq = groupCNVeventsByFrequency()
print("CNVeventByFreq = ", CNVeventByFreq)
printCNVevents(CNVeventByFreq)
mapCNVeventsByGene()


fw09 = open("CNVstats_" + kit + ".csv", "w")
print("==================================================\n")
fw09.write("==================================================\n")
print("sample count = ", len(samples))
fw09.write("sample count = " + "\t" + str(len(samples)))
total_cnvs = 0 
for c in sorted(cnvtypes):
    print("CNV_type = ", c, "counts = ", cnvtypes[c])
    fw09.write("CNV_type = " +"\t" + c +  "\tcounts = " + str(cnvtypes[c]) + "\n")
    total_cnvs += cnvtypes[c]
    print("average ", c , " per sample = ", cnvtypes[c]/len(samples))
    fw09.write("average " + "\t" + c + "\tper sample = " + str(cnvtypes[c]/len(samples)) + "\n")
print("==================================================\n")
fw09.write("==================================================\n")
print("total_cnvs = ", total_cnvs)
fw09.write("total_cnvs = " + "\t" + str(total_cnvs))
#print("known_dgv_count = ", known_count)
#print("novel_dgv_count = ", novel_count)


keys = CNV2chr.keys()
for k in sorted(keys):
    print("chr" + k + " = " + str(CNV2chr[k]) + " density by chr_size = " + str(chr_size[k]/CNV2chr[k])) 

fw10 = open("sample2cnvStats_" + kit + ".csv", "w")
keys = sample2CNVs.keys()
for k in sorted(keys):
    fw10.write(sample2source[k] + "\t" + k)
    for c in chrlist:
        if c in sample2CNVs[k]:
            fw10.write("\t" + str(sample2CNVs[k][c]))
        else:
            fw10.write("\t" + "0")
    fw10.write("\n")
fw10.close()
