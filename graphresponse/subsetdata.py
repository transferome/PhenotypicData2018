"""Load in and subset the data by the different treatments"""
import matplotlib.pyplot as plt


replicateA = [line.rstrip('\n').replace('"', '') for line in open('resources/Simulans2018ResponseData.csv')
              if line.replace('"', '').split(',')[0].endswith('a')]

replicateB = [line.rstrip('\n').replace('"', '') for line in open('resources/Simulans2018ResponseData.csv')
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
        r_list.append(float(line.split(',')[12]))
    return r_list


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
ax.set_title('D. simulans')
for key in dictA.keys():
    ax.scatter(get_gens(dictA[key]), get_resp(dictA[key]))
    ax.plot(get_gens(dictA[key]), get_resp(dictA[key]))
    ax.scatter(get_gens(dictB[key]), get_resp(dictB[key]))
    ax.plot(get_gens(dictB[key]), get_resp(dictB[key]))


fig.show()
test = get_cums(dictA['control'])
test1 = get_resp(dictA['control'])

