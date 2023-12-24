#%% runcell 0 - workflow basics

import regex
from math import log 

prefix= 'AAGAAAGAAATGGTCTGTGATCCCCC'  # GC1
prefix = prefix.upper()

suffix = 'caTTCCCGGCTACAAGGACCCTTCG'  #GC2
suffix = suffix.upper()

# GC1 CTTCCCAGGCCTGCAGTTTGCCCATC   
# GC2 GCTACAAGGACCCTTCGAgccccgttc

# 018 F tCATGGTGCAGTGTAGCCGGGAATG   # ovo su flankovi u nekim plazmidima
# 018 R GGGGGATCACAGACCATTCTCGGCT 

class Read:

    def __init__( self, sequence, name, quality, no_repeats = None, direction = None, gaps_bool = None, gaps = None, flanks = None, repeats = None, read_accuracy = None, index = None, correctedSequence = None ):
        self.sequence = sequence # cela sekvenca
        self.name = name
        self.quality = quality # phred score
        self.no_repeats = no_repeats 
        self.direction = direction
        self.gaps = gaps
        self.gaps_bool = gaps_bool
        self.flanks = flanks
        self.repeats = repeats
        self.read_accuracy = read_accuracy
        self.index = index
        self.correctedSequence = correctedSequence

    def checkForward(self): 

        seq_temp = self.sequence

        prefix_temp = "(%s){s<=5}"%(prefix)
        res_prefix = regex.findall( prefix_temp, seq_temp )
        
        if len(res_prefix)>1:
            self.direction = None
        if len(res_prefix)<1:
            self.direction = None

        suffix_temp = "(%s){s<=5}"%(suffix)
        res_suffix = regex.findall( suffix_temp, seq_temp )
        if len(res_suffix)>1:
            self.direction = None
        if len(res_suffix)<1:
            self.direction = None

   
        if (len(res_suffix)==1 and len(res_prefix)==1):
            self.direction = 'forward'
            self.flanks = ( res_prefix[0], res_suffix[0] )
            
        if (len(res_suffix)>1 and len(res_prefix)>1):
            self.direction = 'multiple'
            self.flanks = ( res_prefix[0], res_suffix[0] )
            
    def checkReverse(self): 
        if self.direction != 'forward':
            seq_temp = self.sequence
            seq_temp= seq_temp[::-1].translate(str.maketrans('ATCG','TAGC')) #vraca reverzni komplement
            prefix_temp = "(%s){s<=5}"%(prefix)
            res_prefix = regex.findall( prefix_temp, seq_temp )
            
            if len(res_prefix)>1:
                self.direction = None
            if len(res_prefix)<1:
                self.direction = None
        
            suffix_temp = "(%s){s<=5}"%(suffix)
            res_suffix = regex.findall( suffix_temp, seq_temp )
            if len(res_suffix)>1:
                self.direction = None
            if len(res_suffix)<1:
                self.direction = None
        
            
            if (len(res_suffix)==1 and len(res_prefix)==1):
                self.direction = 'reverse'
                flank_prefix_temp = res_suffix[0]
                flank_suffix_temp = res_prefix[0]
                self.flanks = (flank_prefix_temp[::-1].translate(str.maketrans('ATCG','TAGC')) , flank_suffix_temp[::-1].translate(str.maketrans('ATCG','TAGC')) )
                
            if (len(res_suffix)>1 and len(res_prefix)>1):
                self.direction = 'multiple'
                self.flanks = ((res_suffix[0])[::-1].translate(str.maketrans('ATCG','TAGC')), (res_prefix[0])[::-1].translate(str.maketrans('ATCG','TAGC'))) # sufiksu i prefiksu su obrnuta mesta, i takodje se prevode u reverzni komplement

    def findRepeats(self): #pronalazenje regiona sa ponovcima, definisanog kao sekvenca izmedju flankova
        seq_temp = self.sequence
        flank_prefix = self.flanks[0]
        flank_suffix = self.flanks[1]
        pattern = flank_prefix+'(.*)'+flank_suffix
        match = regex.search( pattern, seq_temp )
        if match == None:
            self.repeats = ('',0,0)
        if match != None:
            repeat_seq = regex.findall( pattern, seq_temp)
            repeat_position = match.span() 
            self.repeats = (repeat_seq[0],repeat_position[0],repeat_position[1])

    def countRepeats(self): #prebrojavanje neprekinutih readova i dodeljivanje vrednosti gaps_bool koja oznacava da li read ima prekide ili ne
        if self.direction=='reverse':
            target = 'CTG'
        if self.direction=='forward':
            target = 'CAG'
        repeat_seq = self.repeats[0]
        number_of_repeats = repeat_seq.count(target)
        self.no_repeats = number_of_repeats
        
        if len(self.repeats[0])!=self.no_repeats*3:
            self.gaps_bool = True
        if len(self.repeats[0])==self.no_repeats*3:
            self.gaps_bool = False

    def findGaps(self):
        import re
        repeat_seq = self.repeats[0]
        if self.direction == 'forward':
            pattern = re.compile(r'(CAG)+') # patern za pretragu - bilo koji broj CTG-ova
        if self.direction == 'reverse':
            pattern = re.compile(r'(CTG)+')
               
        # pretrazi match za pattern CTG
        ctg_blokovi = pattern.finditer(repeat_seq) 
        pozicije_ctg = []  # lista u kojoj se cuvaju tuplovi sa pozicijama CTG
        
        for ctg_blok in ctg_blokovi:
            #pronalazi tacne pozicije CTG ponovaka unutar read-a
            range_ctg = ctg_blok.span() # za svaki CTG motiv izbacuje tupl sa range-om gde se on nalazi, zakljucno sa tim brojevima 
            pozicije_ctg.append(range_ctg)
        
        pozicije_i_sekvence_prekida = []
        
        if pozicije_ctg ==[]:
            pozicije_i_sekvence_prekida.append((repeat_seq, 0, len(repeat_seq)))
            self.gaps = pozicije_i_sekvence_prekida
            return
        
        # pozivamo funkciju samo za one koji imaju gaps_bool = True, sto znaci da ima bar jedan prekid
        
        # provera prvog CTG bloka da li je spoljasnji ili unutrasnji
        prvi_ctg = pozicije_ctg[0]
        if prvi_ctg[0]==0:
            pass
        else: 
            prekid = repeat_seq[0:prvi_ctg[0]]
            pozicije_i_sekvence_prekida.append((prekid, 0, prvi_ctg[0]))
    
        # ulazimo u loop za sve ostale
        for i in range(0, len(pozicije_ctg)):
            # da posle ne bismo imali out of bounds gresku
            if i == len(pozicije_ctg)-1:
                break
    
            # pozicije_ctg je lista tuplova sa mestima CTG blokova
            # tupl 1  je prvi CTG, tupl2 je drugi CTG
            # uzmi tuplove jedan za drugim 
            tupl1 = pozicije_ctg[i]
            tupl2 = pozicije_ctg[i+1]
            
            # sekvenca prekida je ona izmedju kraja prvog i pocetka drugog CTG bloka
            sekvenca_prekida = repeat_seq[tupl1[1]:tupl2[0]]
            pozicije_i_sekvence_prekida.append((sekvenca_prekida, tupl1[1], tupl2[0]))
    
        # za poslednji CTG proveravamo da li je spoljasnji ili unutrasnji
        poslednji_ctg = pozicije_ctg[-1]
        if poslednji_ctg[1]==len(repeat_seq):
            pass
        else:
            prekid = repeat_seq[poslednji_ctg[1]:]
            pozicije_i_sekvence_prekida.append((prekid, poslednji_ctg[1], len(repeat_seq)))
        
        self.gaps = pozicije_i_sekvence_prekida  
    
    
    def getReadAccuracy(self): #racunanje preciznosti ocitavanja na osnovu Q scora na nivou baza
        quality_temp = self.quality
        quality_num = [(ord(char)-33) for char in list(quality_temp)]
        self.read_accuracy = -10 * log(sum([10**(q / -10) for q in quality_num]) / len(quality_num), 10)   
    

path_to_file="C:/Users/lana/Desktop/workflow/heterozigoti_pcr/big_fastq_heterozigoti.fastq" #ucitavanje i iteracija kroz sam file
contents = None
with open(path_to_file) as f:
    contents = f.readlines()

contents_metadata = contents[0::4]
contents_reads = contents[1::4]
contents_quality = contents[3::4]

names = [read.split()[0] for read in contents_metadata]

reads_from_fastq=[]
for i in range (0,len(contents_metadata)):
    reads_from_fastq.append((contents_reads[i][:-1],names[i][:-1],contents_quality[i][:-1]))
reads_list=[]

i=0
for sequence, name, quality in reads_from_fastq: # petlja koja prolazi kroz sve readove u fastq fajlu i treba da pozove sve odg. metode poreklom od klase Read
    read_class_instance = Read(sequence, name, quality)
    read_class_instance.checkForward()
    read_class_instance.checkReverse()
    read_class_instance.index=i
    i=i+1
#    read_class_instance.checkReverse()
    if (read_class_instance.direction=='forward' or read_class_instance.direction=='reverse'): 
        read_class_instance.findRepeats()
        read_class_instance.countRepeats()
        read_class_instance.findGaps()
        read_class_instance.getReadAccuracy()
    reads_list.append(read_class_instance) #za svaki read koji je imao flankove, u listu reads_list ce biti dodat gap sa informacijama

########

on_target_reads = [read for read in reads_list if read.flanks!=None]

wt_reads = [read for read in on_target_reads if ((read.repeats is not None)and(len(read.repeats[0])<112))] # 37CTG ili manje
premut_reads = [read for read in on_target_reads if ((read.repeats is not None)and(len(read.repeats[0])>=112)and(len(read.repeats[0])<=150))] # do 50 CTG
mut1_reads = [read for read in on_target_reads if ((read.repeats is not None)and(len(read.repeats[0])>150)and(len(read.repeats[0])<=300))] # do 100
mut2_reads = [read for read in on_target_reads if ((read.repeats is not None)and(len(read.repeats[0])>300)and(len(read.repeats[0])<=1000))] # do 330
mut3_reads = [read for read in on_target_reads if ((read.repeats is not None)and(len(read.repeats[0])>1000)and(len(read.repeats[0])<=2500))] # 830 CTG
over_reads = [read for read in on_target_reads if ((read.repeats is not None)and(len(read.repeats[0])>2500))]

#pravljenje fastq fajlova za svaku kategoriju

def saveCategorizedReads(read_category, reads_to_write):
    indices = [read.index for read in read_category]
    for i in indices:
        reads_to_write.append(contents_metadata[i])
        reads_to_write.append(contents_reads[i])
        reads_to_write.append('+\n')
        reads_to_write.append(contents_quality[i])

# wild type
wt_reads_to_write = []
saveCategorizedReads(wt_reads, wt_reads_to_write)
file = open('wt_reads.fastq', 'w')
file.writelines(wt_reads_to_write)
file.close

# premutation
premut_reads_to_write = []
saveCategorizedReads(premut_reads, premut_reads_to_write)
file = open('premut_reads.fastq', 'w')
file.writelines(premut_reads_to_write)
file.close

# mutation 1
mut1_reads_to_write = []
saveCategorizedReads(mut1_reads, mut1_reads_to_write)
file = open('mut1_reads.fastq', 'w')
file.writelines(mut1_reads_to_write)
file.close

# mutation 2
mut2_reads_to_write = []
saveCategorizedReads(mut2_reads, mut2_reads_to_write)
file = open('mut2_reads.fastq', 'w')
file.writelines(mut2_reads_to_write)
file.close

# mutation 3
mut3_reads_to_write = []
saveCategorizedReads(mut3_reads, mut3_reads_to_write)
file = open('mut3_reads.fastq', 'w')
file.writelines(mut3_reads_to_write)
file.close

# over
over_reads_to_write = []
saveCategorizedReads(over_reads, over_reads_to_write)
file = open('over_reads.fastq', 'w')
file.writelines(over_reads_to_write)
file.close

#%% Import svih paketa
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter


#%% Funkcije za statistiku (i generalno olaksavanje zivota)

def nacrtaj_scatterplot_broja_ponovaka(lista_readova, *ime_readova):
    # scatter plot
    svi_brojevi_ponovaka = [x.no_repeats for x in lista_readova if x.no_repeats != None]
    plt.scatter(range(len(svi_brojevi_ponovaka)),svi_brojevi_ponovaka)
    plt.ylabel('Repeat number')
    plt.xlabel('Read id')
    
    if ime_readova: 
        tekst = f'Repeat numbers {len(svi_brojevi_ponovaka)} in {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = f'Repeat numbers {len(svi_brojevi_ponovaka)} reads'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()


def nacrtaj_scatterplot_duzina_readova(lista_readova, *ime_readova):
    sve_duzine = [len(x.sequence) for x in lista_readova]
    plt.scatter(range(len(sve_duzine)),sve_duzine)
    plt.ylabel('Repeat length')
    plt.xlabel('Read id')
    
    if ime_readova: 
        tekst = f'Read lengths on {len(sve_duzine)} in {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = f'Read lengths on {len(sve_duzine)} on-target reads'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
    
    
def nacrtaj_scatter_ctg_vs_read_length(lista_readova, *ime_readova):
    svi_brojevi_ponovaka = [x.no_repeats for x in lista_readova if x.no_repeats != None]
    sve_duzine = [len(x.sequence) for x in lista_readova if x.no_repeats != None]
    plt.scatter(svi_brojevi_ponovaka,sve_duzine)
    plt.xlabel('Repeat number')
    plt.ylabel('Read lenght')
    
    
    if ime_readova: 
        tekst = f'Scatter plot Repeat number vs repeat lenght in {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = 'Scatter plot Repeat number vs repeat lenght'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
    
    
    
def nacrtaj_counts_histogram(lista_readova, *ime_readova):
    svi_brojevi_ponovaka = [x.no_repeats for x in lista_readova if x.no_repeats != None]
    bins = list(range(round(min(svi_brojevi_ponovaka))-1, round(max(svi_brojevi_ponovaka))+1))
    plt.hist(svi_brojevi_ponovaka, bins)
    plt.xlabel('Repeat number')
    plt.ylabel('Number of reads')
   
    
  
    if ime_readova: 
        tekst = f'Histogram broja ponovaka in {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = 'Histogram broja ponovaka'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
   
  
    
def nacrtaj_density(lista_readova, *ime_readova):
    svi_brojevi_ponovaka = [x.no_repeats for x in lista_readova if x.no_repeats != None]
    # plotting histogram and density 
    sns.distplot(a=svi_brojevi_ponovaka)
    plt.xlabel('Repeat number')
    # visualizing plot using matplotlib.pyplot library
    
    if ime_readova: 
        tekst = f'Density plot in {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = 'Density plot'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
    
    
def nacrtaj_quality_histogram(lista_readova, *ime_readova): 
    
    Q_scores = [x.read_accuracy for x in lista_readova if x.read_accuracy != None]
    bins = list(range(round(min(Q_scores))-1, round(max(Q_scores))+1))
    plt.hist(Q_scores, bins)
    plt.xlabel('Quality score')
    
    if ime_readova: 
        tekst = f'Histogram skorova kvaliteta za {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = 'Histogram skorova kvaliteta'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
    
    
    
def nacrtaj_bar_plot_broja_ponovaka(lista_readova, *ime_readova):

    brojevi_ponovaka = [x.no_repeats for x in lista_readova]
    brojevi_ponovaka.sort()

    element_counts = Counter(brojevi_ponovaka)
    elements = list(element_counts.keys())
    counts = list(element_counts.values())
    
    # Create a bar plot
    plt.bar(elements, counts)
    plt.xlabel('broj ponovaka')
    
    if ime_readova: 
        tekst = f'Distribucija duzina alela od ukupno {len(brojevi_ponovaka)} read-ova za {ime_readova[0][0]}'
        plt.title(tekst)
    else: 
        tekst = f'Distribucija duzina alela od ukupno {len(brojevi_ponovaka)} read-ova'
        plt.title(tekst)
    ime_slike = tekst+'.jpg'
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
    

def nacrtaj_sve_grafike_EDA(lista_readova, *ime_readova):
    
    if ime_readova:
        nacrtaj_scatterplot_broja_ponovaka(lista_readova,ime_readova )
        nacrtaj_scatterplot_duzina_readova(lista_readova,ime_readova)
        nacrtaj_scatter_ctg_vs_read_length(lista_readova,ime_readova)
        nacrtaj_counts_histogram(lista_readova,ime_readova)
        nacrtaj_density(lista_readova,ime_readova)
        nacrtaj_quality_histogram(lista_readova,ime_readova)
        return 
    
    # scatter plotovi svih brojeva ponovaka u on-target readovima
    nacrtaj_scatterplot_broja_ponovaka(lista_readova)

    # Scatterplot svih duzina read-ova
    nacrtaj_scatterplot_duzina_readova(lista_readova)

    # broj ponovaka vs duzina read-ova
    nacrtaj_scatter_ctg_vs_read_length(lista_readova)

    # histogram broja ponovaka svih on target readova
    nacrtaj_counts_histogram(lista_readova)

    # density plot with histogram za on target reads
    nacrtaj_density(lista_readova)

    # Histogram skorova kvaliteta svih on-target read-ova 
    nacrtaj_quality_histogram(lista_readova)



def nacrtaj_stacked_barplot_prekida(lista_readova, df_ucestalosti):
    
    brojevi_prekida = [len(x.gaps) for x in lista_readova]
    # grupe su brojevi prekida 
    group_names = list(set(brojevi_prekida))

    # categories su brojevi ponovaka. Treba nam da prebrojimo koliko reads po grupi ima koji broj ponovaka
    categories = list(df_ucestalosti['Broj ponovaka'])
    num_categories = len(categories)

    groups = []
    #group_values = [0]*len(brojevi ponovaka)
    for group in group_names:
        group_name = f'{group} prekida'
        group_values = []
        reads = [x.no_repeats for x in lista_readova if len(x.gaps)==group]
        element_counts = Counter(reads)
        for i in categories: 
        #    print(i)
             group_values.append(element_counts[i])
             
        groups.append({'label': group_name, 'values': group_values})     
        
    bottom_values = [0] * num_categories

    for group in groups:
        plt.bar(categories, group['values'], label=group['label'], bottom=bottom_values)
        bottom_values = [b + v for b, v in zip(bottom_values, group['values'])]

    # Adding labels and title
    plt.xlabel('Brojevi ponovaka')
    plt.ylabel('Broj readova')
    plt.title('Broj prekida u zavisnosti od broja ponovaka')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    plt.savefig('Broj prekida u zavisnosti od broja ponovaka.jpg', format='jpeg')
    plt.show()

def find_outliers_IQR(df):
    #import pandas
    q1=df.quantile(0.25)
    q3=df.quantile(0.75)
    IQR=q3-q1
    outliers = df[((df<(q1-1.5*IQR)) | (df>(q3+1.5*IQR)))]
    return outliers


def find_main_aleles(df):
    # uzima df u kojoj imamo kolonu '%' i kolonu 'Broj ponovaka'
    # ideja ove funkcije je da olaksa identifikaciju klastera alela u uzorku sa nepoznatim brojem pomesanih osoba/alelskih distribucija. 
    # pretpostavka je da ce aleli koji su mod distibucije biti outlieri u smislu ucestalosti sa kojom se javljaju, tj. bice ih dosta vise od ostalih
    
    udeli_cistih_ctg = df['%']
    outliers = find_outliers_IQR(udeli_cistih_ctg)
    
    #nadji broj ponovaka koji su outlieri po ucestalosti
    main_ponovci = []
    for outlier in outliers: 
        # nadji indeks gde je on u df_pure['%']
        a = df.index[df['%'] == outlier].tolist()
        broj_ponovaka_outlier = df.at[a[0], 'Broj ponovaka']
        main_ponovci.append(broj_ponovaka_outlier)
    return main_ponovci    


def nacrtaj_heatmap(reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs):
     
    a = [x for x in reads_list if x.no_repeats == broj_ponovaka_za_plotovanje]
    print(len(a))
    if len(a) == 0:
        return
        
    complement_dict = {'A': '1', 'T': '2', 'C': '3', 'G': '4'}
    
    vectors=[]
    for read in a:
        vector = list('0'*len(read.repeats[0]))
        #print(vector)
        prekidi = read.gaps
        
        for ii in range(len(prekidi)):
            prekid = prekidi[ii]  # za svaki prekid
            seq = prekid[0] # izvuci sekvencu prekida
            #print(seq)
            numeroid = ''.join(complement_dict.get(base, base) for base in seq)  # koristi dictionary odozgo i menja nukleotide brojevima
            vector[prekid[1]:prekid[2]] = list(numeroid) # zameni nule brojevima odozgo
            vector = [int(digit) for digit in vector]
            #print(vector, len(vector))
        vectors.append(vector) # ima vektore razlicitih velicina
        
    # Pad vectors with zeros to create a matrix
    matrix = np.array([vector + [-1] * (max(map(len, vectors)) - len(vector)) for vector in vectors])
    print(max(map(len, vectors)))
    # Create a DataFrame from the matrix
    df = pd.DataFrame(matrix)
    #print(df)
    df = df.sort_values(by=kolona_za_sortiranje, ascending = True)

    custom_colors = ['#000000', '#D2D1CD', '#F0E01C',  '#1CA0F0', '#EB3016','#23C64F']
    # #000000 - crna - rupa za popunjavanje prostora 
    # #D2D1CD - siva - nule (gde je CTG)
    # #F0E01C - zuta - adenin 
    # #1CA0F0 - plava - timin
    # #EB3016 - crvena - citozin
    # #23C64F - zelena - guanin
    
    #xticklabels = ["/", "CTG", "A", "T", "C", "G"]
    
    sns.set()
    heatmap = sns.heatmap(df, cmap = sns.color_palette(custom_colors))
    
    cbar = heatmap.collections[0].colorbar
    cbar.set_ticklabels(["/", "CTG", "A", "T", "C", "G"])
    if kwargs:
        plt.xlim(0, kwargs[0])
    ime_slike = f'Heatmap_{broj_ponovaka_za_plotovanje}_repeats.jpg'
    plt.title(ime_slike)
    plt.savefig(ime_slike, format = 'jpeg')
    plt.show()
    
def napravi_df_ucestalosti_prekida(lista_readova, broj_ponovaka_za_plotovanje):
    
    
    a = [x.gaps[0] for x in lista_readova if x.no_repeats == broj_ponovaka_za_plotovanje]
    a_dict = Counter(a)

    df_prekidi = pd.DataFrame()
    df_prekidi['prekid'] = list(a_dict.keys())
    df_prekidi['ucestalost'] = list(a_dict.values())
    df_prekidi = df_prekidi.sort_values(by='ucestalost', ascending =False)
    
    ime_tabele=f'Ucestalost_prekida_{broj_ponovaka_za_plotovanje}_ctg_ponovaka.csv'
    df_prekidi.to_csv(ime_tabele)
    return df_prekidi   
    

def izracunaj_sve(lista_readova, kategorija):

    """  
    kategorija mora biti string koji ce biti deo naslova nekih grafika i imena fajlova
    npr. kategorija = 'mut2_reads'
    Automatski se cuva sve 
    """
    # nacrtaj sve grafike za EDA
    nacrtaj_sve_grafike_EDA(lista_readova, kategorija) 

    # crtamo dataframe za ucestalosti brojeva ponovaka
    brojevi_ponovaka = [x.no_repeats for x in lista_readova]                            #   <======== !!!
    brojevi_ponovaka.sort()

    element_counts = Counter(brojevi_ponovaka)
    elements = list(element_counts.keys())
    counts = list(element_counts.values())


    df_category  = pd.DataFrame()
    df_category['Broj ponovaka'] = elements
    df_category['Ucestalost'] = counts 
    df_category['%'] = [x/len(brojevi_ponovaka)*100 for x in counts]  
    
    ime_tabele = 'Repeat_frequency_'+ kategorija +'.csv'
    df_category.to_csv(ime_tabele)                                       #   <======== !!!
    print(f"Ucestalost broja ponovaka u {kategorija}")
    print(df_category) 

    brojevi_ponovaka_pure = [x.no_repeats for x in lista_readova if x.gaps_bool == False]
    brojevi_ponovaka_pure.sort()

    element_counts = Counter(brojevi_ponovaka_pure)
    elements = list(element_counts.keys())
    counts = list(element_counts.values())

    df_category_pure  = pd.DataFrame()
    df_category_pure['Broj ponovaka'] = elements
    df_category_pure['Ucestalost'] = counts 
    df_category_pure['%'] = [x/len(brojevi_ponovaka_pure)*100 for x in counts]  
    
    ime_tabele = 'Repeat_frequency_'+ kategorija +'_pure.csv'
    df_category_pure.to_csv(ime_tabele)                                       #   <======== !!!
    print(f"Ucestalost broja cistih ponovaka u {kategorija}")
    print(df_category_pure) 
                        

    # plotovi udela zastupljenosti readova sa dati brojem ponovaka u kategoriji read-ova
    
    x1 = df_category['Broj ponovaka']
    y1 = df_category["%"]
    
    x2 = df_category_pure['Broj ponovaka']
    y2 = df_category_pure["%"]
    
    fig, axs = plt.subplots(1,2, figsize=(10,4))
    axs[0].scatter(x1, y1)
    axs[0].set_title(' All on-target reads')
    axs[0].set_ylabel('% of repeat number')
    
    axs[1].scatter(x2, y2)
    axs[1].set_title('Pure on-target reads')
    axs[1].set_ylabel('% of repeat number')
    
    ime_slike = f'Udeo CTG ponovaka u {kategorija}.jpg'
    plt.savefig(ime_slike,format='jpeg' )
    plt.show()

    # nalazenje outlier-a -> identifikacija wt clastera

    print(f'\nMain aleles in all reads {find_main_aleles(df_category)}')
    print(f'Main aleles in pure reads {find_main_aleles(df_category_pure)}')

    #Stacked bar plor prekida
    # ne moze biti df_pure_wt jer tu nema prekida za plotovanje
    nacrtaj_stacked_barplot_prekida(lista_readova, df_category)
    
    
    cat_sa_prekidima = [x for x in lista_readova if x.gaps_bool == True]
    
    
    if kategorija == 'wt_reads':
        start = 0; stop = 37
    elif kategorija == 'premut_reads':
        start = 37; stop = 50
    elif kategorija == 'mut1_reads':
        start = 50; stop = 100
    elif kategorija == 'mut2_reads':
        start = 101; stop = 330
    elif kategorija == 'mut3_reads':
        start = 330; stop = 830
    
    
    
    for ctg in range(start,stop+1):
        broj_ponovaka_za_plotovanje = ctg
        kolona_za_sortiranje = 0 # default
       
        # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
        nacrtaj_heatmap(cat_sa_prekidima, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
        napravi_df_ucestalosti_prekida(cat_sa_prekidima, broj_ponovaka_za_plotovanje) 


    return df_category, df_category_pure

#%%============================================================================ 
# Distribucije i procenti - TOTAL READS
#==============================================================================

# brojevi readova po kategorijama
broj_ukupno_reads = len(reads_list)
broj_on_target = len(on_target_reads)
broj_wt = len(wt_reads)
broj_premut = len(premut_reads)
broj_mut1 = len(mut1_reads)
broj_mut2 = len(mut2_reads)
broj_mut3 = len(mut3_reads)
broj_over = len(over_reads)

svi_cisti = len([x for x in reads_list if x.gaps_bool == False])
svi_on_target_cisti = len([x for x in on_target_reads if x.gaps_bool == False])
broj_cistih_wt = len([x for x in wt_reads if x.gaps_bool == False])
broj_cistih_premut = len([x for x in premut_reads if x.gaps_bool == False])
broj_cistih_mut1 = len([x for x in mut1_reads if x.gaps_bool == False])
broj_cistih_mut2 = len([x for x in mut2_reads if x.gaps_bool == False])
broj_cistih_mut3 = len([x for x in mut3_reads if x.gaps_bool == False])
broj_cistih_over = len([x for x in over_reads if x.gaps_bool == False])


# procenti i cuvanje u dataframe-u
df_brojevi_procenti = pd.DataFrame()
df_brojevi_procenti['Kategorija'] = ['Total reads', 'On-target reads', 'WT reads', 'Premut reads', 'Mut1 reads', 'Mut2 reads', 'Mut3 reads', 'Over 2.5kb reads']
df_brojevi_procenti['Reads number'] = [broj_ukupno_reads,broj_on_target, broj_wt,broj_premut,broj_mut1, broj_mut2, broj_mut3, broj_over]
df_brojevi_procenti['% in On Target reads'] = [x/broj_on_target*100 for x in df_brojevi_procenti['Reads number']]
df_brojevi_procenti['% in Total reads'] = [x/broj_ukupno_reads*100 for x in df_brojevi_procenti['Reads number']]

df_brojevi_procenti['Pure Reads number'] = [svi_cisti, svi_on_target_cisti, broj_cistih_wt,broj_cistih_premut ,broj_cistih_mut1, broj_cistih_mut2, broj_cistih_mut3, broj_cistih_over]
df_brojevi_procenti['% pure'] = [x/svi_cisti*100 for x in df_brojevi_procenti['Pure Reads number']]


# sacuvaj dataframe u fajlu 
df_brojevi_procenti.to_csv('Reads by category.csv')

# prikazi dataframe
print(df_brojevi_procenti)

# nacrtaj sve grafike za EDA
nacrtaj_sve_grafike_EDA(on_target_reads, 'on-target reads')



# za kompletnu automatizaciju odkomentovati ovo ispod: 
    
#izracunaj_sve(wt_reads, 'wt_reads')
#izracunaj_sve(premut_reads, 'premut_reads')
#izracunaj_sve(mut1_reads, 'mut1_reads')
#izracunaj_sve(mut2_reads, 'mut2_reads')
#izracunaj_sve(mut3_reads, 'mut3_reads')

#%%============================================================================ 
#         WT READS  
#==============================================================================

# nalaze se u listi wt_reads

# nacrtaj sve grafike za EDA
nacrtaj_sve_grafike_EDA(wt_reads, 'wt_reads')

# crtamo dataframe za ucestalosti brojeva ponovaka
brojevi_ponovaka = [x.no_repeats for x in wt_reads]
brojevi_ponovaka.sort()

element_counts = Counter(brojevi_ponovaka)
elements = list(element_counts.keys())
counts = list(element_counts.values())
                     
df_wt  = pd.DataFrame()
df_wt['Broj ponovaka'] = elements
df_wt['Ucestalost'] = counts 
df_wt['%'] = [x/len(brojevi_ponovaka)*100 for x in counts]  
#df_wt.to_csv('RepeatFrequency_wt.csv')                                         #   <======== !!!
print("Ucestalost broja ponovaka u svim wt read-ovima")
print(df_wt)

#  df_pure je od wt ponovaka                     

wt_brojevi_ponovaka_pure = [x.no_repeats for x in wt_reads if x.gaps_bool == False]
wt_brojevi_ponovaka_pure.sort()

element_counts = Counter(wt_brojevi_ponovaka_pure)
elements = list(element_counts.keys())
counts = list(element_counts.values())

df_pure_wt  = pd.DataFrame()
df_pure_wt['Broj ponovaka'] = elements
df_pure_wt['Ucestalost'] = counts 
df_pure_wt['%'] = [x/len(wt_brojevi_ponovaka_pure)*100 for x in counts]  
#df_pure_wt.to_csv('RepeatFrequency_wt_pure_CTG.csv')                           #   <======== !!!
print("\nUcestalost broja ponovaka u cistim wt read-ovima")
print(df_pure_wt)

# plotovi udela zastupljenosti readova sa dati brojem ponovaka u kategoriji read-ova

x1 = df_wt['Broj ponovaka']
y1 = df_wt["%"]

x2 = df_pure_wt['Broj ponovaka']
y2 = df_pure_wt["%"]

fig, axs = plt.subplots(1,2, figsize=(10,4))
axs[0].scatter(x1, y1)
axs[0].set_title(' All on-target reads')
axs[0].set_ylabel('% of repeat number')


axs[1].scatter(x2, y2)
axs[1].set_title('Pure on-target reads')
axs[1].set_ylabel('% of repeat number')
plt.savefig('Udeo CTG ponovaka u dobijenim read-ovima.jpg',format='jpeg' )
plt.show()

# nalazenje outlier-a -> identifikacija wt clastera

print(f'\nMain aleles in all wt_reads {find_main_aleles(df_wt)}')
print(f'Main aleles in pure wt_reads {find_main_aleles(df_pure_wt)}')

#Stacked bar plor prekida
# ne moze biti df_pure_wt jer tu nema prekida za plotovanje
nacrtaj_stacked_barplot_prekida(wt_reads, df_wt)

#%% WT grupa 4-5 CTG ponovaka

wt_do_5_ctg = [x for x in wt_reads if x.no_repeats <= 5 if x.gaps_bool == True]

kolona_za_sortiranje = 0
# reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje,                *kwargs <- xlimit

for ctg in range(6):
    broj_ponovaka_za_plotovanje = ctg
    kolona_za_sortiranje = 0 # default
   
    # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
    nacrtaj_heatmap(wt_do_5_ctg, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
    napravi_df_ucestalosti_prekida(wt_do_5_ctg, broj_ponovaka_za_plotovanje)  
    
#%% WT grupa 6-25 ponovaka

wt_5_25_ctg = [x for x in wt_reads if x.no_repeats > 5 if x.no_repeats<=25 if x.gaps_bool == True ]

kolona_za_sortiranje = 0
# reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje,                *kwargs <- xlimit

for ctg in range(6,26):
    broj_ponovaka_za_plotovanje = ctg
    kolona_za_sortiranje = 0 # default
   
    # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
    nacrtaj_heatmap(wt_5_25_ctg, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
    napravi_df_ucestalosti_prekida(wt_5_25_ctg, broj_ponovaka_za_plotovanje)  
    

#%% wt aleli 26-37

wt_26_37_ctg = [x for x in wt_reads if x.no_repeats > 26 if x.gaps_bool == True ]
kolona_za_sortiranje = 0
# reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje,                *kwargs <- xlimit

for ctg in range(26,38):
    broj_ponovaka_za_plotovanje = ctg
    kolona_za_sortiranje = 0 # default
   
    # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
    nacrtaj_heatmap(wt_26_37_ctg, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
    napravi_df_ucestalosti_prekida(wt_26_37_ctg, broj_ponovaka_za_plotovanje) 

#%%============================================================================ 
#         PREMUT READS
#==============================================================================


# nalaze se u listi premut_reads

# nacrtaj sve grafike za EDA
nacrtaj_sve_grafike_EDA(premut_reads, 'premut_reads')                           #   <======== !!!

# crtamo dataframe za ucestalosti brojeva ponovaka
brojevi_ponovaka = [x.no_repeats for x in premut_reads]                         #   <======== !!!
brojevi_ponovaka.sort()

element_counts = Counter(brojevi_ponovaka)
elements = list(element_counts.keys())
counts = list(element_counts.values())
                     
df_premut  = pd.DataFrame()
df_premut['Broj ponovaka'] = elements
df_premut['Ucestalost'] = counts 
df_premut['%'] = [x/len(brojevi_ponovaka)*100 for x in counts]  
#df_premut.to_csv('RepeatFrequency_premutacija.csv')                                         #   <======== !!!
print("Ucestalost broja ponovaka u svim premut read-ovima")
print(df_premut)

#  df_pure je od wt ponovaka                     

premut_brojevi_ponovaka_pure = [x.no_repeats for x in premut_reads if x.gaps_bool == False]
premut_brojevi_ponovaka_pure.sort()

element_counts = Counter(premut_brojevi_ponovaka_pure)
elements = list(element_counts.keys())
counts = list(element_counts.values())

df_pure_premut  = pd.DataFrame()
df_pure_premut['Broj ponovaka'] = elements
df_pure_premut['Ucestalost'] = counts 
df_pure_premut['%'] = [x/len(premut_brojevi_ponovaka_pure)*100 for x in counts]  
#df_pure_premut.to_csv('RepeatFrequency_premut_pure_CTG.csv')                           #   <======== !!!
print("\nUcestalost broja ponovaka u cistim premut read-ovima")
print(df_pure_premut)

# plotovi udela zastupljenosti readova sa dati brojem ponovaka u kategoriji read-ova

x1 = df_premut['Broj ponovaka']
y1 = df_premut["%"]

x2 = df_pure_premut['Broj ponovaka']
y2 = df_pure_premut["%"]

fig, axs = plt.subplots(1,2, figsize=(10,4))
axs[0].scatter(x1, y1)
axs[0].set_title(' All on-target reads')
axs[0].set_ylabel('% of repeat number')


axs[1].scatter(x2, y2)
axs[1].set_title('Pure on-target reads')
axs[1].set_ylabel('% of repeat number')
plt.savefig('Udeo CTG ponovaka u dobijenim premut read-ovima.jpg',format='jpeg' )
plt.show()

# nalazenje outlier-a -> identifikacija wt clastera

print(f'\nMain aleles in all premut_reads {find_main_aleles(df_premut)}')
print(f'Main aleles in pure premut_reads {find_main_aleles(df_pure_premut)}')

#Stacked bar plor prekida
# ne moze biti df_pure_wt jer tu nema prekida za plotovanje
nacrtaj_stacked_barplot_prekida(premut_reads, df_premut)

#%% Heatmaps

premut_sa_prekidima = [x for x in premut_reads if x.gaps_bool == True]

for ctg in range(10,49):
    broj_ponovaka_za_plotovanje = ctg
    kolona_za_sortiranje = 0 # default
   
    # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
    nacrtaj_heatmap(premut_sa_prekidima, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
    napravi_df_ucestalosti_prekida(premut_sa_prekidima, broj_ponovaka_za_plotovanje) 



#%%============================================================================ 
#         MUT1 READS
#==============================================================================


# nalaze se u listi premut_reads

# nacrtaj sve grafike za EDA
nacrtaj_sve_grafike_EDA(mut1_reads, 'mut1_reads')                                #   <======== !!!

# crtamo dataframe za ucestalosti brojeva ponovaka
brojevi_ponovaka = [x.no_repeats for x in mut1_reads]                            #   <======== !!!
brojevi_ponovaka.sort()

element_counts = Counter(brojevi_ponovaka)
elements = list(element_counts.keys())
counts = list(element_counts.values())
                     
df_mut1  = pd.DataFrame()
df_mut1['Broj ponovaka'] = elements
df_mut1['Ucestalost'] = counts 
df_mut1['%'] = [x/len(brojevi_ponovaka)*100 for x in counts]  
df_mut1.to_csv('RepeatFrequency_mut1.csv')                                       #   <======== !!!
print("Ucestalost broja ponovaka u svim mut1 read-ovima")
print(df_mut1)

#  df_pure je od wt ponovaka                     

mut1_brojevi_ponovaka_pure = [x.no_repeats for x in mut1_reads if x.gaps_bool == False]
mut1_brojevi_ponovaka_pure.sort()

element_counts = Counter(mut1_brojevi_ponovaka_pure)
elements = list(element_counts.keys())
counts = list(element_counts.values())

df_pure_mut1  = pd.DataFrame()
df_pure_mut1['Broj ponovaka'] = elements
df_pure_mut1['Ucestalost'] = counts 
df_pure_mut1['%'] = [x/len(mut1_brojevi_ponovaka_pure)*100 for x in counts]  
#df_pure_mut1.to_csv('RepeatFrequency_mut1_pure_CTG.csv')                   #   <======== !!!
print("\nUcestalost broja ponovaka u cistim mut1 read-ovima")
print(df_pure_mut1)

# plotovi udela zastupljenosti readova sa dati brojem ponovaka u kategoriji read-ova

x1 = df_mut1['Broj ponovaka']
y1 = df_mut1["%"]

x2 = df_pure_mut1['Broj ponovaka']
y2 = df_pure_mut1["%"]

fig, axs = plt.subplots(1,2, figsize=(10,4))
axs[0].scatter(x1, y1)
axs[0].set_title(' All on-target reads')
axs[0].set_ylabel('% of repeat number')


axs[1].scatter(x2, y2)
axs[1].set_title('Pure on-target reads')
axs[1].set_ylabel('% of repeat number')
plt.savefig('Udeo CTG ponovaka u dobijenim mut1 read-ovima.jpg',format='jpeg' )
plt.show()

# nalazenje outlier-a -> identifikacija wt clastera

print(f'\nMain aleles in all mut1_reads {find_main_aleles(df_mut1)}')
print(f'Main aleles in pure mut1_reads {find_main_aleles(df_pure_mut1)}')

#Stacked bar plor prekida
# ne moze biti df_pure_wt jer tu nema prekida za plotovanje
nacrtaj_stacked_barplot_prekida(mut1_reads, df_mut1)


mut1_sa_prekidima = [x for x in mut1_reads if x.gaps_bool == True]

for ctg in range(50,101):
    broj_ponovaka_za_plotovanje = ctg
    kolona_za_sortiranje = 0 # default
   
    # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
    nacrtaj_heatmap(mut1_sa_prekidima, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
    napravi_df_ucestalosti_prekida(mut1_sa_prekidima, broj_ponovaka_za_plotovanje) 




#%%============================================================================ 
#         MUT2 READS
#==============================================================================


# nalaze se u listi premut_reads

# nacrtaj sve grafike za EDA
nacrtaj_sve_grafike_EDA(mut2_reads, 'mut2_reads')                                #   <======== !!!

# crtamo dataframe za ucestalosti brojeva ponovaka
brojevi_ponovaka = [x.no_repeats for x in mut2_reads]                            #   <======== !!!
brojevi_ponovaka.sort()

element_counts = Counter(brojevi_ponovaka)
elements = list(element_counts.keys())
counts = list(element_counts.values())
                     
df_mut2  = pd.DataFrame()
df_mut2['Broj ponovaka'] = elements
df_mut2['Ucestalost'] = counts 
df_mut2['%'] = [x/len(brojevi_ponovaka)*100 for x in counts]  
df_mut2.to_csv('RepeatFrequency_mut2.csv')                                       #   <======== !!!
print("Ucestalost broja ponovaka u svim mut2 read-ovima")
print(df_mut2)

#  df_pure je od wt ponovaka                     

mut2_brojevi_ponovaka_pure = [x.no_repeats for x in mut2_reads if x.gaps_bool == False]
mut2_brojevi_ponovaka_pure.sort()

element_counts = Counter(mut2_brojevi_ponovaka_pure)
elements = list(element_counts.keys())
counts = list(element_counts.values())

df_pure_mut2  = pd.DataFrame()
df_pure_mut2['Broj ponovaka'] = elements
df_pure_mut2['Ucestalost'] = counts 
df_pure_mut2['%'] = [x/len(mut2_brojevi_ponovaka_pure)*100 for x in counts]  
#df_pure_mut2.to_csv('RepeatFrequency_mut2_pure_CTG.csv')                   #   <======== !!!
print("\nUcestalost broja ponovaka u cistim mut2 read-ovima")
print(df_pure_mut2)

# plotovi udela zastupljenosti readova sa dati brojem ponovaka u kategoriji read-ova

x1 = df_mut2['Broj ponovaka']
y1 = df_mut2["%"]

x2 = df_pure_mut2['Broj ponovaka']
y2 = df_pure_mut2["%"]

fig, axs = plt.subplots(1,2, figsize=(10,4))
axs[0].scatter(x1, y1)
axs[0].set_title(' All on-target reads')
axs[0].set_ylabel('% of repeat number')


axs[1].scatter(x2, y2)
axs[1].set_title('Pure on-target reads')
axs[1].set_ylabel('% of repeat number')
plt.savefig('Udeo CTG ponovaka u dobijenim mut2 read-ovima.jpg',format='jpeg' )
plt.show()

# nalazenje outlier-a -> identifikacija wt clastera

print(f'\nMain aleles in all mut2_reads {find_main_aleles(df_mut2)}')
print(f'Main aleles in pure mut2_reads {find_main_aleles(df_pure_mut2)}')

#Stacked bar plor prekida
# ne moze biti df_pure_wt jer tu nema prekida za plotovanje
nacrtaj_stacked_barplot_prekida(mut2_reads, df_mut2)


mut2_sa_prekidima = [x for x in mut2_reads if x.gaps_bool == True]

for ctg in range(101,331):
    broj_ponovaka_za_plotovanje = ctg
    kolona_za_sortiranje = 0 # default
   
    # reads_list, broj_ponovaka_za_plotovanje, kolona_za_sortiranje, *kwargs <- xlimit
    nacrtaj_heatmap(mut2_sa_prekidima, broj_ponovaka_za_plotovanje, kolona_za_sortiranje)
    napravi_df_ucestalosti_prekida(mut2_sa_prekidima, broj_ponovaka_za_plotovanje) 













