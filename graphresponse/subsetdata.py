"""Load in and subset the data by the different treatments"""
import math


def get_gen0sd(species='Simulans'):
    if species == 'Simulans':
        data = [line.rstrip('\n') for line in open('resources/Simulans2018Master.csv') if 'Gen0a' in line]
        datb = [line.rstrip('\n') for line in open('resources/Simulans2018Master.csv') if 'Gen0b' in line]
    if species == 'Melanogaster':
        data = [line.rstrip('\n') for line in open('resources/SelectionMelanogasterCombinedData.csv')
                if 'Gen0a' in line]
        datb = [line.rstrip('\n') for line in open('resources/SelectionMelanogasterCombinedData.csv')
                if 'Gen0b' in line]
    dindexa = [float(line.split(',')[201]) for line in data]
    ma = sum(dindexa)/len(dindexa)
    vara = sum((xi - ma) ** 2 for xi in dindexa)/len(dindexa)
    dindexb = [float(line.split(',')[201]) for line in datb]
    mb = sum(dindexb)/len(dindexb)
    varb = sum((xi - mb) ** 2 for xi in dindexb)/len(dindexb)
    # print(vara, varb)
    return math.sqrt(vara), math.sqrt(varb)


def create_dictionary(species='Simulans'):
    """ Create the data dictionary """
    if species == 'Simulans':
        replicateA = [line.rstrip('\n').replace('"', '') for line in open('resources/Simulans2018ResponseData.csv')
                      if line.replace('"', '').split(',')[0].endswith('a')]

        replicateB = [line.rstrip('\n').replace('"', '') for line in open('resources/Simulans2018ResponseData.csv')
                      if line.replace('"', '').split(',')[0].endswith('b')]
    if species == 'Melanogaster':
        replicateA = [line.rstrip('\n').replace('"', '') for line in open('resources/SelectionResponseMelanogaster2018.csv')
                      if line.replace('"', '').split(',')[0].endswith('a')]
        replicateB = [line.rstrip('\n').replace('"', '') for line in
                      open('resources/SelectionResponseMelanogaster2018.csv')
                      if line.replace('"', '').split(',')[0].endswith('b')]
    dictA = {'control': list(), 'up1': list(), 'up2': list(), 'down1': list(), 'down2': list()}
    dictB = {'control': list(), 'up1': list(), 'up2': list(), 'down1': list(), 'down2': list()}
    for key in dictA.keys():
        for line in replicateA:
            if key in line:
                dictA[key].append(line)
    for key in dictB.keys():
        for line in replicateB:
            if key in line:
                dictB[key].append(line)
    return dictA, dictB


def get_gens(sublist):
    g_list = list()
    for line in sublist:
        g_list.append(int(line.split(',')[6]))
    return g_list


def get_cums(sublist):
    s_list = list()
    for line in sublist:
        s_list.append(float(line.split(',')[11]))
    return s_list


def get_resp(sublist):
    r_list = list()
    for line in sublist:
        r_list.append(float(line.split(',')[13]))
    return r_list


def get_n(sublist):
    n_list = list()
    for line in sublist:
        n_list.append(float(line.split(',')[9]))
    return n_list


def get_var(sublist):
    v_list = list()
    for line in sublist:
        v_list.append(float(line.split(',')[8]))
    return v_list


def make_error(sublist, spec):
    sa, sb = get_gen0sd(species=spec)
    se_list = list()
    for line in sublist:
        # stndard error in selection index is sqrt(Vp)/sqrt(N)
        # but gets divided by sd at gen0, either for replicate A or B
        temp_val = math.sqrt(float(line.split(',')[8]))/math.sqrt(float(line.split(',')[9]))
        if line.split(',')[0].endswith('a'):
            se_list.append((1.96*temp_val)/sa)
        if line.split(',')[0].endswith('b'):
            se_list.append((1.96*temp_val)/sb)
    # 95% confidence interval is +/- these values
    return se_list


if __name__ == '__main__':
    pass
