#!/bin/env python3

import numpy as np

#========#
#========#

def reading_law(filename):

    xval = []
    yval = []
    dataset = []
    data = []

    with open(filename) as f:
        for line in f:
            if line=='\n':
                continue
            elif line=='\t\n':
                continue
            elif line.startswith('#'):
                continue
            else:
                xi,yi=line.split()

                xi = float(xi)
                yi = float(yi)

                if xi == 0.0 and len(xval)>0:
                    dataset = [xval,yval]
                    data.append(dataset)
                    xval = []
                    yval = []

                xval.append(float(xi))
                yval.append(float(yi))

    #-- Last Dataset append --#
    dataset = [xval,yval]
    data.append(dataset)

    f.close()

    return data

#-------#

def reading_data_Xeval(lawdata,filename):

    xval = []
    yval = []
    dataset = []
    data = []
    index=0
    time0=0
    time_set=False

    with open(filename) as f:
        for line in f:
            if line == '\n':
                continue
            elif line == '\t\n':
                continue
            elif line.startswith('#'):
                continue
            else:
                xi,yi=line.split()

                xi = float(xi)
                yi = float(yi)

                # Flag for resetting the timer to zero
                if time_set == True and xi < time0:
                    time_set=False

                # Determines if you either have a dataset or not
                if func_inter(xi, lawdata[index]) >= 5:
                    law5_flag = True
                else:
                    law5_flag = False

                # Is read data in range of the flag or not
                if law5_flag==True:

                    if time_set==False:
                        time0    = xi
                        xi       = 0.0
                        time_set = True

                    if xi == 0.0 and len(xval)>0:
                        dataset = [xval,yval]
                        data.append(dataset)
                        xval = []
                        yval = []
                        index += 1

                    if xi == 0.0:
                        xval.append(float(xi))
                        yval.append(float(yi))
                    else:
                        xval.append(float(xi)-time0)
                        yval.append(float(yi))

    #-- Last Dataset append --#
    dataset = [xval,yval]
    data.append(dataset)

    f.close()

    return data

#========#
#========#

def read_data(filename):

    mydata1 = {}
    mydata2 = {}
    header_dict = {}
    flag_header = False
    idcounter = 0
    linecount = 0
    flag_id = False
    flag_output = False
    output = []


    # Determine EOF
    with open(filename) as f:
        f.readlines()
        numEOF = f.tell()
    f.close()

    with open(filename) as f:
        for line in f:
            if line=='\n' or line.startswith('='):

                if line.startswith('===T') and flag_header==False:
                    header_list = line.split()
                    header_list[0] = 'time'
                    for i,val in enumerate(header_list):
                        mydata1.update({val: []})
                        mydata2.update({val: []})
                        header_dict.update({val: []})
                    # adding an additional id for post processing
                    #mydata.update({'id':[]})
                    #mydata.update({'count':[]})
                    flag_header=True

                elif line.startswith('===EOF'): # and len(mydata1['time'])>10:
                    #--------------------------------------------------------------#
                    flag_id = False

                    output.append({})

                    for i, val in enumerate(header_list):
                        output[idcounter].update({val:[]})
                        for i in range(len(mydata1[val])):
                            output[idcounter][val].append(mydata1[val][i])

                        mydata1[val].clear()
                    idcounter += 1
                    #--------------------------------------------------------------#

                else:
                    continue

            else:
                tmp = [float(n) for n in line.split()]

                if (float(tmp[0]) == 0.0 and len(mydata1['time'])>10):

                    # --------------------------------------------------------------#
                    flag_id = False

                    output.append({})

                    for i, val in enumerate(header_list):
                        output[idcounter].update({val:[]})
                        for i in range(len(mydata1[val])):
                            output[idcounter][val].append(mydata1[val][i])

                        mydata1[val].clear()
                    idcounter += 1
                    #--------------------------------------------------------------#


                if float(tmp[0]) == 0.0 and flag_id == False:
                    flag_id = True
                    flag_output = True
                    for i,val in enumerate(header_list):
                        mydata1[val].append(tmp[i])
                        mydata2[val].append(tmp[i])
                elif float(tmp[0]) == 0.0 and flag_id == True:
                    continue
                else:
                    for i,val in enumerate(header_list):
                        mydata1[val].append(tmp[i])
                        mydata2[val].append(tmp[i])

                #elif float(tmp[0]) == 0.0 and flag_id == False:
                #    flag_id = True

                #for i,val in enumerate(header_list):

                    #print(i,val,counter)
                    #if val=='id':
                    #    if float(tmp[0]) == 0.0 and idcounter != 0:
                    #        idcounter += 1
                    #        linecount = 0
                    #    elif float(tmp[0]) == 0.0 and linecount != 0:
                    #        linecount += 1
                    #    mydata[val].append(idcounter)

                #    mydata1[val].append(tmp[i])
                #    mydata2[val].append(tmp[i])

                    #if float(tmp[0]) == 0.0  and flag_id == False:
                    #    flag_id = True
                    #    flag_header = False
                    #    mydata[val].append(tmp[i])
                    #elif float(tmp[0]) != 0.0:
                    #    # std appending of data
                    #    mydata[val].append(tmp[i])
                    #elif float(tmp[0]) == 0.0  and flag_id == True:
                    #    output.append(mydata)
                    #    # cleanup
                    #    for i,val in enumerate(header_list):
                    #        mydata[val].clear()
                    #    flag_id = False
                    #    mydata[val].append(tmp[i])



    return mydata2,output

#========#
#========#

def read_xy(filename):

    mydata = {}
    header_list = {'time':[],'pmass':[]}
    flag_header = False
    #counter = 0
    with open(filename) as f:
        for line in f:
            #counter += 1
            if line=='\n':
               continue
            else:
                tmp = [float(n) for n in line.split()]

                for i,val in enumerate(header_list):
                    #print(i,val,counter)
                    mydata[val].append(tmp[i])
    return mydata

#========#
#========#

def averages(threshold,crit,filter,data):

    tmp = []
    averages = {}
    medians = {}
    maxi = {}
    mini = {}

    for i,val in enumerate(filter):
        averages.update({val:0})
        medians.update({val:0})

    for i,val in enumerate(crit):
        maxi.update({val:0})
        mini.update({val:0})

    for j, val in enumerate(crit):
        maxi[val] = np.array(data[val]).max()
        mini[val] = np.array(data[val]).min()
        critmax = maxi[val] * threshold
        critmin = mini[val] * threshold

        if critmax > 0 and critmin >= 0:
            crit[val] = critmax
        elif critmax <= 0 and critmin < 0:
            crit[val] = critmin
        else:
            raise IOError('Min and Max problem')

    # --------------------
    for j, val in enumerate(filter):

        for i in range(len(data[val])):
            # test = filter[val]
            # test = crit[filter[val]]
            # checks effF_CO2 -> dXdt_CO2 ...
            if abs(data[filter[val]][i]) >= abs(crit[filter[val]]):
                # appends data[effF ... C1 ... ]
                dud = data[val][i]
                tmp.append(data[val][i])

        averages[val] = np.average(tmp)
        medians[val] = np.median(tmp)


        # in average -- every value is the same -- spikes included
        #averages[val] = np.average(tmp)

        tmp = []

    return averages,medians,maxi

#========#
#========#

def func_inter(timex,data):
    time_array = data[0]
    func_array = data[1]
    funcX = np.interp(timex,time_array,func_array)
    return float(funcX)

#========#
#========#

#####################################################
## START
#####################################################

threshold = 0.25

crit = {
    'dXdt_O2':  0,
    'dXdt_CO2': 0,
    'dXdt_H2O': 0,
}

#-------------------------------
filter = {
    'effF_O2' :'dXdt_O2',
    'effF_CO2':'dXdt_CO2',
    'effF_H2O':'dXdt_H2O',
    'C1_O2'   :'dXdt_O2',
    'C1_CO2'  :'dXdt_CO2',
    'C1_H2O'  :'dXdt_H2O',
}

#-------------------------------
data_names = [
    "./output/PEFR_inj1.out",
    "./output/PEFR_inj2.out",
    "./output/PEFR_inj3.out",
    "./output/PEFR_inj4.out",
    "./output/PEFR_inj5.out",
]

pmass_files = [
    "./data/inj1/inj1_p_mass.xy",
    "./data/inj2/inj2_p_mass.xy",
    "./data/inj3/inj3_p_mass.xy",
    "./data/inj4/inj4_p_mass.xy",
    "./data/inj5/inj5_p_mass.xy",
]

plaw_files = [
    "./data/inj1/inj1_p_law.xy",
    "./data/inj2/inj2_p_law.xy",
    "./data/inj3/inj3_p_law.xy",
    "./data/inj4/inj4_p_law.xy",
    "./data/inj5/inj5_p_law.xy",
]

##----------------------------
averages_collection = {}
for i,val in enumerate(data_names):
    averages_collection.update({val:0})

#----------------------------
#-- OUTPUT ------------------
#----------------------------
#----------------------------
for i,name in enumerate(data_names):
    #--------------------
    # one,two read_data()
    # one .. total set
    # two .. individual sets
    zeroData1 = read_data(data_names[i])[0]

    #--------------------

    averages_collection[name] = averages(threshold,crit,filter,zeroData1)

print('-----------------------------------')
print('-- average rate values ------------')
print('-----------------------------------')
for i, name in enumerate(data_names):
    print(name,averages_collection[name][0])

print('-----------------------------------')
print('-- media rate values ------------')
print('-----------------------------------')
for i, name in enumerate(data_names):
    print(name,averages_collection[name][1])

print('-------------------------------')
print('-- max rate values ------------')
print('-------------------------------')
for i,name in enumerate(data_names):
    print(name,averages_collection[name][2])

print('------------------------------------')
print('-- dmdt fluent/0D ratio ------------')
print('------------------------------------')


#lawdataset = []
for i,name in enumerate(plaw_files):
    # dataset
    # dataset[] .. datafiles (law1 law2 law3 law4 ..)
    # dataset[][] .. datasets in each datafile (law1: 15 datasets)
    # dataset[][][0] .. x in the dataset
    # dataset[][][1] .. y in the dataset
    # dataset[][][0][] .. x value
    lawdataset = reading_law(name)
    pmass      = reading_data_Xeval(lawdataset,pmass_files[i])
    zeroData2  = read_data(data_names[i])[1]
    av_dict    = { 'deltadmpdt':[] }

    # [ time , pmass , dmdt_fluent , dmdt_zeroD , dmdt_fluent/dmdt_zeroD ]
    #pmass_array = np.zeros([len(pmass) , len(pmass[k][0]) , 5])

    for k in range(len(pmass)):
        pmass_array = np.zeros([len(pmass[k][0]), 5])

        for j in range(len(pmass[k][0])):

            pmass_array[j][0] = float(pmass[k][0][j])
            pmass_array[j][1] = float(pmass[k][1][j])

            if j == 0:
                continue
            else:
                if (pmass_array[j][0]-pmass_array[j-1][0]) == 0.0:
                    pmass_array[j][2] = 0.0
                else:
                    # 2..fluent
                    # 3..zeroDim
                    pmass_array[j][2] = (pmass_array[j][1]-pmass_array[j-1][1])/(pmass_array[j][0]-pmass_array[j-1][0])

                    if zeroData2[k]['time'][-1] < pmass_array[j][0]:
                        pmass_array[j][3] = 0.0
                    else:
                        pmass_array[j][3] = func_inter(pmass_array[j][0],[zeroData2[k]['time'],zeroData2[k]['dmpdt']])

                    if pmass_array[j][3] == 0.0:
                        pmass_array[j][4] = 0.0
                    else:
                        pmass_array[j][4] = pmass_array[j][2]/pmass_array[j][3]

                    av_dict['deltadmpdt'].append(pmass_array[j][4])

#        pmass_array = np.zeros([len(pmass), len(pmass[k][0]), 5])

#        for j in range(len(pmass[k][0])):
#            print(k,j)
#            pmass_array[k][j][0] = float(pmass[k][0][j])
#            pmass_array[k][j][1] = float(pmass[k][1][j])
#
#            if j == 0:
#                continue
#            else:
#                if (pmass_array[k][j][0]-pmass_array[k][j-1][0]) == 0.0:
#                    pmass_array[k][j][2] = 0.0
#                else:
#                    # 2..fluent
#                    # 3..zeroDim
#                    pmass_array[k][j][2] = (pmass_array[k][j][1]-pmass_array[k][j-1][1])/(pmass_array[k][j][0]-pmass_array[k][j-1][0])
#
#                    if zeroData2[k]['time'][-1] < pmass_array[k][j][0]:
#                        pmass_array[k][j][3] = 0.0
#                    else:
#                        pmass_array[k][j][3] = func_inter(pmass_array[k][j][0],[zeroData2[k]['time'],zeroData2[k]['dmpdt']])
#
#                    if pmass_array[k][j][3] == 0.0:
#                        pmass_array[k][j][4] = 0.0
#                    else:
#                        pmass_array[k][j][4] = pmass_array[k][j][2]/pmass_array[k][j][3]
#
#                    av_dict['dmpdt'].append(pmass_array[k][j][4])

    # todo vorzeichen von deltadmpdt may has an influence, see averages!!!
    av_values = averages(0.02,{'deltadmpdt':0},{'deltadmpdt':'deltadmpdt'},av_dict)
    print( 'average',av_values[0],'median',av_values[1] )


#    print('yep')



    #lawdataset.append(reading_law(name))

#pmass = []
#for i,name in enumerate(pmass_files):
#    pmass.append(reading_data_Xeval(lawdataset[i],name))
#
#
#
##lawdataset = reading_law('inj3_p_law.xy')
#
## anzahl datensaetze datensatz(t,valuo)->2,values
##pmass = reading_data_Xeval(lawdataset,'inj3_p_mass.xy')
##zeroData1,zeroData2 = read_data('PEFR_inj3.out')
#
##tmp_list = []
##tmp_list_dataset = []
#
##for i,val in enumerate(zeroData):
##
##    if val == 'time'
##        if float(zeroData[val]) == 0.0 and len(tmp_list_dataset[0]) > 10:
##            tmp_list.append(tmp_list_dataset)
##            tmp_list_dataset = []
##        else:
#
## Array time,mass,dmdt(<- starts at 1 instead of 0)
## Anzahl datensaetze, anzahl daten in datensaetze, welche daten
#pmass_array = np.zeros([len(pmass)*len(pmass[0]),len(pmass[0][0][0]),3])
##zeroD_array = np.zeros([len(zeroData['time'])])
#
#
#
#for i in range (len(pmass_array)):
#    for j in range (len(pmass_array[i])):
#        pmass_array[i][j][0] = pmass[i][0][j]
#        pmass_array[i][j][1] = pmass[i][1][j]
#
#        if j == 0:
#            continue
#        else:
#            pmass_array[i][j][2] = (pmass_array[i][j][1]-pmass_array[i][j-1][1])/(pmass_array[i][j][0]-pmass_array[i][j-1][0])


print('-- finished --')



